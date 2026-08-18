#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""P0-3c：鹊桥二号（Queqiao-2）月球大椭圆冻结轨道运行动画（基于 poliastro）。

轨道根数参考 ``analysis/elevation_video/queqiao_relay.py::ORBIT_PARAMS``
（证据等级 A，来自孟占峰 2024 论文表 1 微分修正结果）：
- 半长轴 a = 10004.059 km，偏心率 e = 0.782356
- 轨道倾角 I = 119.203°（逆行冻结轨道），近月点幅角 ω = 90°
- Ω = 0°、平近点角 M = 215°（证据等级 B 的假设值，与 queqiao_relay.py 一致；
  由 M 解 Kepler 方程得真近点角 ν 后构造轨道）
- 周期：论文 24.84 h 为测控回归周期 T_R = 2π/(Ṁ''+ω̇'')（含地球三体摄动长期项）；
  本动画按月心二体传播，显示开普勒周期 24.94 h（poliastro 二体值）
- 24.84 vs 24.94 h 差异（约 6.1 分钟）主因 = 地球三体引力摄动（数值验证：
  J₂+地球三体全模型近点周期 24.840 h，与论文吻合到 0.01 分钟；仅 J₂ 贡献 <0.1 s
  可忽略）。二者是不同物理量（回归周期 vs 开普勒周期），非误差，不影响科学结论。

动画如何体现近月点/远月点的快慢：
1. 尾迹按**等时间间隔**采样 —— 近月点处点间距大（飞得快）、远月点处点间距小（飞得慢）
2. 尾迹按速度着色（coolwarm：红=快，蓝=慢）
3. 卫星头部带速度矢量箭头（长度正比于速率）
4. 速度–真近点角曲线（光标指示当前位置）+ 信息面板实时参数
   （高度 / 速度 / 距近月点 / 距远月点），并标注近/远月点速度对比

版面（--layout，手机横屏与竖屏都适配，关键数字均用大字号高亮卡片突出）：
- **portrait（默认，9:16 竖屏 1080×1920）**：标题 → 3D 轨道（放大置顶，--zoom 1.9，
  远月点弧段裁出画面以突出近月段大椭圆）→ 速度曲线 → 底部信息卡片四行。
- **landscape（16:9 横屏 1920×1080）**：左侧 3D 轨道大图，右上标题 + 右侧
  信息卡片列，底部通栏速度曲线。

用法::

    # 手机竖屏（默认）→ figures-and-logs/queqiao_orbit_animation.mp4
    python analysis/elevation_video/queqiao_orbit_animation.py --layout portrait

    # 手机横屏 → figures-and-logs/queqiao_orbit_animation_landscape.mp4
    python analysis/elevation_video/queqiao_orbit_animation.py --layout landscape

    # 只渲染首帧 PNG 检查版面
    python analysis/elevation_video/queqiao_orbit_animation.py --layout portrait --preview-only
    python analysis/elevation_video/queqiao_orbit_animation.py --layout landscape --preview-only

    # 自定义
    python analysis/elevation_video/queqiao_orbit_animation.py \
        --output-dir figures-and-logs --orbits 2 \
        --frames-per-orbit 360 --fps 30 --tail 90

要求：conda 环境 ``spaceport``（Python 3.10，poliastro 0.17.0 + astropy 5.3.4）::

    conda run -n spaceport python analysis/elevation_video/queqiao_orbit_animation.py

"""

from __future__ import annotations

import argparse
import math
import os
import sys

import numpy as np
from astropy import units as u
from astropy.time import Time

import matplotlib

matplotlib.use("Agg")
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.animation import FFMpegWriter, FuncAnimation
from matplotlib.patches import FancyBboxPatch

# poliastro 0.17（Python 3.10 的 spaceport 环境）
from poliastro.bodies import Moon
from poliastro.ephem import Ephem
from poliastro.twobody import Orbit

# 允许以脚本方式从仓库任意位置运行：将仓库根加入 sys.path
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from analysis.elevation_video.queqiao_relay import (  # noqa: E402
    M_ASSUMED_DEG,
    ORBIT_PARAMS,
    RAAN_ASSUMED_DEG,
)
from analysis.elevation_video.utils import MOON_R_KM  # noqa: E402

DEFAULT_OUTPUT_DIR = os.path.join(_REPO_ROOT, "logs", "elevation_video")

# ---------------------------------------------------------------------------
# 版面（portrait = 手机竖屏 9:16，landscape = 手机横屏 16:9）
# ---------------------------------------------------------------------------
# box_aspect 与 lim_scale 成比例（x:y:z 数据跨度 ∝ 屏幕盒子比），保证月球不变形：
#   portrait  xlim=±1.0lim  ylim=±1.15lim  zlim=±0.85lim（2.0:2.3:1.7 = 1:1.15:0.85）
#   landscape xlim=±1.25lim ylim=±1.0lim   zlim=±0.75lim（2.5:2:1.5 = 1.25:1:0.75）
LAYOUTS: dict = {
    "portrait": {
        "figsize": (9.0, 16.0),           # 1080×1920 @120dpi
        "kind": "stacked",                # 竖屏：标题→3D→曲线→信息 四行
        "height_ratios": [0.14, 2.40, 0.75, 1.60],
        "hspace": 0.16,
        "left": 0.03, "right": 0.97, "top": 0.975, "bottom": 0.025,
        "box_aspect": (1.0, 1.15, 0.85),
        "lim_scale": (1.0, 1.15, 0.85),
        "view_init": (18, 145),
        "default_zoom": 1.9,
        "sat_size": 190, "tail_size": 24, "vec_lw": 3.0,
        "marker_size": 90, "sat_lw": 1.8,
        "legend": False,                  # 图内已有近/远月点标注，去图例更简洁
        "normal_line": False,             # 去掉轨道面法向参考线
        "label_3d_fs": 18,
        "orbit_lw": 1.5,
        "card_value_fs": 32,
        "card_dist_fs": 24,
        "card_label_fs": 18,
    },
    "landscape": {
        "figsize": (16.0, 9.0),           # 1920×1080 @120dpi
        "kind": "split",                  # 横屏：左 3D 大图 + 右信息卡片 + 底通栏曲线
        "width_ratios": [1.5, 1.0],
        "height_ratios": [0.4, 1.2, 0.8],
        "hspace": 0.1, "wspace": 0.1,
        "left": 0.1, "right": 0.9, "top": 0.9, "bottom": 0.1,
        "box_aspect": (1.25, 1.0, 0.75),
        "lim_scale": (1.25, 1.0, 0.75),
        "view_init": (20, 150),
        "default_zoom": 2.5,
        "sat_size": 140, "tail_size": 18, "vec_lw": 2.6,
        "marker_size": 80, "sat_lw": 1.7,
        "legend": False,
        "normal_line": False,
        "label_3d_fs": 12,
        "orbit_lw": 1.5,
        "card_value_fs": 27,
        "card_dist_fs": 20,
        "card_label_fs": 16,
        "ticks_fs": 14,
    },
}

# ---------------------------------------------------------------------------
# 中文字体（macOS 常见 CJK 字体，按优先级 fallback）
# ---------------------------------------------------------------------------
# 简体中文字体优先（PingFang HK 等繁体字体会缺简体字形，如「远」）
_CJK_FONTS = [
    # "Alibaba PuHuiTi 3.0",
    # "Hiragino Sans GB",
    "Songti SC",
    "SimHei",
    # "PingFang SC",
    # "STHeiti",
    # "PingFang HK",
    "Heiti TC",
    "Arial Unicode MS",
    "Noto Sans CJK SC",
]


def setup_cjk_font() -> str | None:
    """配置 matplotlib 中文字体，返回选中的字体名（找不到返回 None）。"""
    avail = {f.name for f in fm.fontManager.ttflist}
    chosen = next((c for c in _CJK_FONTS if c in avail), None)
    if chosen is None:
        print("警告：未找到中文字体，中文可能显示为方块", file=sys.stderr)
        return None
    plt.rcParams["font.sans-serif"] = [chosen, "DejaVu Sans"]
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["axes.unicode_minus"] = False
    return chosen


# ---------------------------------------------------------------------------
# 轨道构造（根数取自 queqiao_relay.py ORBIT_PARAMS + 假设的 Ω/M）
# ---------------------------------------------------------------------------
def solve_kepler(M_rad: float, e: float, tol: float = 1e-12) -> float:
    """牛顿迭代解 Kepler 方程 E - e·sinE = M，返回偏近点角 E（弧度）。"""
    E = M_rad
    for _ in range(40):
        dE = (E - e * math.sin(E) - M_rad) / (1.0 - e * math.cos(E))
        E -= dE
        if abs(dE) < tol:
            break
    return E


def build_orbit(epoch_iso: str):
    """按 ORBIT_PARAMS 构造 poliastro 月心轨道（平近点角 M → 真近点角 ν）。"""
    a_km = ORBIT_PARAMS["semi_major_axis_km"]
    ecc = ORBIT_PARAMS["eccentricity"]
    inc_deg = ORBIT_PARAMS["inclination_deg"]
    argp_deg = ORBIT_PARAMS["periapsis_argument_deg"]
    raan_deg = RAAN_ASSUMED_DEG
    M_deg = M_ASSUMED_DEG

    # 由平近点角解出真近点角（与 queqiao_relay.build_orbital_state 一致）
    M_rad = math.radians(M_deg)
    E = solve_kepler(M_rad, ecc)
    nu = 2.0 * math.atan2(
        math.sqrt(1.0 + ecc) * math.sin(E / 2.0),
        math.sqrt(1.0 - ecc) * math.cos(E / 2.0),
    )
    nu_deg = math.degrees(nu % (2.0 * math.pi))

    orb = Orbit.from_classical(
        Moon,
        a_km * u.km,
        ecc * u.one,
        inc_deg * u.deg,
        raan_deg * u.deg,
        argp_deg * u.deg,
        nu_deg * u.deg,
        epoch=Time(epoch_iso, scale="utc"),
    )
    return orb


def sample_trajectory(orb, frames_per_orbit: int, n_orbits: int):
    """等时间间隔采样轨道（单周期数组 tile 到多圈），返回轨道平面坐标下的 r/v/ν/速度等。"""
    period_h = orb.period.to_value(u.h)
    epochs = orb.epoch + np.linspace(0, 1, frames_per_orbit, endpoint=False) * orb.period
    ephem = Ephem.from_orbit(orb, epochs)
    r_cycle, v_cycle = ephem.rv()
    r_cycle = r_cycle.to_value(u.km)
    v_cycle = v_cycle.to_value(u.km / u.s)

    r_full = np.tile(r_cycle, (n_orbits, 1))
    v_full = np.tile(v_cycle, (n_orbits, 1))

    # 轨道基：x = 近月点方向，y = 轨道面内垂直，z = 轨道面法向
    h = np.cross(r_cycle[0], v_cycle[0])
    h /= np.linalg.norm(h)
    i_peri = int(np.argmin(np.linalg.norm(r_cycle, axis=1)))
    i_apo = int(np.argmax(np.linalg.norm(r_cycle, axis=1)))
    u_p = r_cycle[i_peri] / np.linalg.norm(r_cycle[i_peri])
    u_q = np.cross(h, u_p)
    u_q /= np.linalg.norm(u_q)
    rot = np.vstack([u_p, u_q, h])  # 行向量为新基

    r_orb = r_full @ rot.T
    v_orb = v_full @ rot.T

    rnorm = np.linalg.norm(r_full, axis=1)
    speeds = np.linalg.norm(v_full, axis=1)
    alts = rnorm - MOON_R_KM

    # 真近点角 ν ∈ [0, 2π)（0° = 近月点）
    r_unit = r_full / rnorm[:, None]
    cos_nu = np.clip(r_unit @ u_p, -1.0, 1.0)
    sin_nu = np.clip(r_unit @ u_q, -1.0, 1.0)
    nu_deg = np.degrees(np.arctan2(sin_nu, cos_nu) % (2.0 * np.pi))

    return {
        "period_h": period_h,
        "r_orb": r_orb,
        "v_orb": v_orb,
        "rnorm": rnorm,
        "speeds": speeds,
        "alts": alts,
        "nu_deg": nu_deg,
        "i_peri": i_peri,
        "i_apo": i_apo,
        "v_peri": float(np.linalg.norm(v_cycle[i_peri])),
        "v_apo": float(np.linalg.norm(v_cycle[i_apo])),
        "h_peri_km": float(np.linalg.norm(r_cycle[i_peri]) - MOON_R_KM),
        "h_apo_km": float(np.linalg.norm(r_cycle[i_apo]) - MOON_R_KM),
    }


# ---------------------------------------------------------------------------
# 绘图
# ---------------------------------------------------------------------------
def draw_moon(ax, R: float) -> None:
    """绘制月球球体 + 经纬网格。"""
    n_th, n_ph = 56, 28
    theta = np.linspace(0, 2 * np.pi, n_th)
    phi = np.linspace(0, np.pi, n_ph)
    xs = R * np.outer(np.cos(theta), np.sin(phi))
    ys = R * np.outer(np.sin(theta), np.sin(phi))
    zs = R * np.outer(np.ones(n_th), np.cos(phi))
    ax.plot_surface(
        xs, ys, zs,
        rstride=2, cstride=2, color="#aab6c8", alpha=0.92,
        edgecolor="#7e8ba0", linewidth=0.15, antialiased=True,
    )


def setup_3d_axes(ax, lim: float, L: dict) -> None:
    """按版面配置 3D 盒体（box_aspect ∝ lim_scale → 月球保持正球）。"""
    ls = L["lim_scale"]
    ax.set_box_aspect(L["box_aspect"])
    ax.set_xlim(-lim * ls[0], lim * ls[0])
    ax.set_ylim(-lim * ls[1], lim * ls[1])
    ax.set_zlim(-lim * ls[2], lim * ls[2])
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_ticks([])
        axis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.view_init(*L["view_init"])


# ---------------------------------------------------------------------------
# 信息面板（两种手机版面：竖屏 stacked / 横屏 split）
# ---------------------------------------------------------------------------
def build_info_portrait(ax, d: dict, L: dict) -> dict:
    """手机竖屏信息面板：大字号 + 高亮色块卡片，突出关键数字。

    布局（transAxes，y 自下而上）：
      顶部副标题（轨道根数）→ 卡片行1（当前高度/速度，32pt 高亮）
      → 卡片行2（距近/远月点，24pt）→ 近/远月点对比 → 圈次/时间/注。
    """
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    def card(x, y, w, h):
        ax.add_patch(FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.018",
            transform=ax.transAxes, facecolor="#eef2f7",
            edgecolor="#c3d0de", linewidth=1.2, zorder=0,
        ))

    def txt(x, y, s, size, color, weight="normal", ha="left", va="center",
            family=None):
        return ax.text(x, y, s, transform=ax.transAxes, fontsize=size,
                       color=color, weight=weight, ha=ha, va=va,
                       family=family, zorder=3)

    # 顶部副标题（轨道根数，简洁一行）
    txt(0.5, 0.952,
        f"a={ORBIT_PARAMS['semi_major_axis_km']:.1f} km   "
        f"e={ORBIT_PARAMS['eccentricity']:.4f}   "
        f"i={ORBIT_PARAMS['inclination_deg']:.1f}°   "
        f"二体周期≈{d['period_h']:.2f} h",
        13, "#455a64", ha="center")
    ax.plot([0.05, 0.95], [0.905, 0.905], transform=ax.transAxes,
            color="#c3d0de", lw=1.0, zorder=1)

    # 卡片行1：当前高度 / 当前速度（最关键数字，最大最亮）
    card(0.02, 0.72, 0.46, 0.185)
    card(0.52, 0.72, 0.46, 0.185)
    txt(0.055, 0.868, "当前高度", 16, "#546e7a")
    t_h = txt(0, 0.772, "", L["card_value_fs"], "#0d47a1", "bold")
    txt(0.555, 0.868, "当前速度", 16, "#546e7a")
    t_v = txt(0.5, 0.772, "", L["card_value_fs"], "#0d47a1", "bold")

    # 卡片行2：距下一近月点 / 距下一远月点
    card(0.02, 0.505, 0.46, 0.185)
    card(0.52, 0.505, 0.46, 0.185)
    txt(0.055, 0.653, "距下一近月点", 16, "#546e7a")
    t_peri = txt(0.055, 0.56, "", L["card_dist_fs"], "#d62728", "bold")
    txt(0.555, 0.653, "距下一远月点", 16, "#546e7a")
    t_apo = txt(0.555, 0.56, "", L["card_dist_fs"], "#1f77b4", "bold")

    # 近/远月点速度对比（静态）
    txt(0.04, 0.41, f"近月点: h={d['h_peri_km']:.0f} km,   v={d['v_peri']:.3f} km/s（快）",
        14.5, "#d62728", "bold")
    txt(0.04, 0.36, f"远月点: h={d['h_apo_km']:.0f} km,    v={d['v_apo']:.3f} km/s（慢）",
        14.5, "#1f77b4", "bold")
    ax.plot([0.02, 0.98], [0.325, 0.325], transform=ax.transAxes, color="#c3d0de", lw=1.0, zorder=1)

    # 圈次 / 时间 / 注
    t_orbit = txt(0.04, 0.262, "", 13, "#37474f")
    t_time = txt(0.04, 0.198, "", 13, "#37474f", family="DejaVu Sans Mono")
    

    return {"t_time": t_time, "t_orbit": t_orbit, "t_h": t_h, "t_v": t_v,
            "t_peri": t_peri, "t_apo": t_apo}


def build_info_landscape(ax, d: dict, L: dict) -> dict:
    """手机横屏信息面板：右侧卡片式，关键数字大字号高亮。

    布局（transAxes）：
      2×2 卡片（高度/速度 27pt 深蓝，距近/远月点 20pt 红/蓝）
      → 近/远月点对比 → 圈次/时间/注。
    """
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    def card(x, y, w, h):
        ax.add_patch(FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.018",
            transform=ax.transAxes, facecolor="#eef2f7",
            edgecolor="#c3d0de", linewidth=1.2, zorder=0,
        ))

    def txt(x, y, s, size, color, weight="normal", ha="left", va="center",
            family=None):
        return ax.text(x, y, s, transform=ax.transAxes, fontsize=size,
                       color=color, weight=weight, ha=ha, va=va,
                       family=family, zorder=3)

    # 2×2 卡片
    card(0.02, 0.74, 0.46, 0.22)
    card(0.52, 0.74, 0.46, 0.22)
    txt(0.055, 0.905, "当前高度", L['card_label_fs'], "#546e7a")
    t_h = txt(0.055, 0.795, "", L["card_value_fs"], "#0d47a1", "bold")
    txt(0.555, 0.905, "当前速度", L['card_label_fs'], "#546e7a")
    t_v = txt(0.555, 0.795, "", L["card_value_fs"], "#0d47a1", "bold")

    card(0.02, 0.48, 0.46, 0.22)
    card(0.52, 0.48, 0.46, 0.22)
    txt(0.055, 0.655, "距下一近月点", L['card_label_fs'], "#546e7a")
    t_peri = txt(0.055, 0.545, "", L["card_dist_fs"], "#d62728", "bold")
    txt(0.555, 0.655, "距下一远月点", L['card_label_fs'], "#546e7a")
    t_apo = txt(0.555, 0.545, "", L["card_dist_fs"], "#1f77b4", "bold")

    # 近/远月点速度对比（静态）
    txt(0.04, 0.415, f"近月点　h={d['h_peri_km']:.0f} km　v={d['v_peri']:.3f} km/s（快）",
        12.5, "#d62728", "bold")
    txt(0.04, 0.355, f"远月点　h={d['h_apo_km']:.0f} km　v={d['v_apo']:.3f} km/s（慢）",
        12.5, "#1f77b4", "bold")
    ax.plot([0.02, 0.98], [0.285, 0.285], transform=ax.transAxes,
            color="#c3d0de", lw=1.0, zorder=1)

    # 圈次 / 时间 / 注
    t_orbit = txt(0.04, 0.215, "", 11.5, "#37474f")
    t_time = txt(0.04, 0.14, "", 11.5, "#37474f", family="DejaVu Sans Mono")

    return {"t_time": t_time, "t_orbit": t_orbit, "t_h": t_h, "t_v": t_v,
            "t_peri": t_peri, "t_apo": t_apo}


def main() -> None:
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    font = setup_cjk_font()
    L = LAYOUTS[args.layout]
    is_portrait = L["kind"] == "stacked"
    zoom = args.zoom if args.zoom is not None else L["default_zoom"]

    # ---- 轨道与采样 -------------------------------------------------------
    print(f"[queqiao_orbit_animation] 构造轨道 (epoch={args.epoch}, layout={args.layout}) ...")
    orb = build_orbit(args.epoch)
    d = sample_trajectory(orb, args.frames_per_orbit, args.orbits)

    N = len(d["r_orb"])
    N_per = args.frames_per_orbit
    dt_h = d["period_h"] / N_per
    epoch = Time(args.epoch, scale="utc")

    print(
        f"[queqiao_orbit_animation] 二体开普勒周期={d['period_h']:.3f} h "
        f"(论文回归周期 24.84 h，差异由地球三体摄动长期项引起), "
        f"近月点 h={d['h_peri_km']:.1f} km v={d['v_peri']:.3f} km/s, "
        f"远月点 h={d['h_apo_km']:.1f} km v={d['v_apo']:.3f} km/s"
    )

    cmap = matplotlib.colormaps["coolwarm"]
    norm = mcolors.Normalize(vmin=d["v_apo"], vmax=d["v_peri"])
    # --zoom 内容放大：lim 缩小 → 轨道与月球在画面中放大
    lim = d["rnorm"].max() * 1.08 / zoom

    # ---- 画布布局 ---------------------------------------------------------
    # portrait：9:16 竖屏 1080×1920，标题→3D 轨道→速度曲线→信息面板 四行；
    # landscape：16:9 横屏 1920×1080，左 3D 大图 + 右上标题 + 右信息卡片 + 底通栏曲线。
    fig = plt.figure(figsize=L["figsize"], dpi=args.dpi)
    if is_portrait:
        gs = fig.add_gridspec(
            4, 1, height_ratios=L["height_ratios"], hspace=L["hspace"],
            left=L["left"], right=L["right"], top=L["top"], bottom=L["bottom"],
        )
        ax_title = fig.add_subplot(gs[0, 0])
        ax_title.axis("off")
        ax_title.text(0.5, 0.5, "鹊桥二号 · 月球大椭圆冻结轨道",
                      ha="center", va="center", fontsize=22, weight="bold",
                      color="#1a237e")
        ax3d = fig.add_subplot(gs[1, 0], projection="3d")
        ax_spd = fig.add_subplot(gs[2, 0])
        ax_info = fig.add_subplot(gs[3, 0])
        ax_info.axis("off")
    else:
        gs = fig.add_gridspec(
            3, 2, width_ratios=L["width_ratios"],
            height_ratios=L["height_ratios"],
            hspace=L["hspace"], wspace=L["wspace"],
            left=L["left"], right=L["right"], top=L["top"], bottom=L["bottom"],
        )
        ax3d = fig.add_subplot(gs[0:2, 0], projection="3d")
        ax_title = fig.add_subplot(gs[0, 1])
        ax_title.axis("off")
        ax_title.text(0.04, 0.60, "鹊桥二号 · 月球大椭圆冻结轨道",
                      transform=ax_title.transAxes, fontsize=22,
                      weight="bold", color="#1a237e", va="center")
        ax_title.text(
            0.04, 0.18,
            f"a={ORBIT_PARAMS['semi_major_axis_km']:.1f} km   "
            f"e={ORBIT_PARAMS['eccentricity']:.4f}   "
            f"i={ORBIT_PARAMS['inclination_deg']:.1f}°   "
            f"二体周期≈{d['period_h']:.2f} h",
            transform=ax_title.transAxes, fontsize=13, color="#455a64",
            va="center")
        ax_info = fig.add_subplot(gs[1, 1])
        ax_info.axis("off")
        ax_spd = fig.add_subplot(gs[2, :])

    # ---- 3D 主图（静态元素） ----------------------------------------------
    draw_moon(ax3d, MOON_R_KM)
    setup_3d_axes(ax3d, lim, L)

    # 轨道背景线（单周期闭环）
    cyc = np.concatenate([d["r_orb"][:N_per], d["r_orb"][:1]])
    (orb_line,) = ax3d.plot(
        cyc[:, 0], cyc[:, 1], cyc[:, 2],
        color="#3e4a59", lw=L["orbit_lw"], alpha=0.85, label="轨道",
    )

    # 轨道面法向参考线（两种手机版面均去掉，更简洁）
    if L["normal_line"]:
        ax3d.plot(
            [0, 0], [0, 0],
            [-MOON_R_KM * 1.25, MOON_R_KM * 1.25],
            color="#5d6d7e", lw=1.0, ls=":", alpha=0.8,
        )
        ax3d.text(0, 0, MOON_R_KM * 1.6, "轨道面法向", color="#5d6d7e",
                  fontsize=9, ha="center")

    # 近月点 / 远月点标记
    # matplotlib ≥3.4 的 Axes3D 默认 computed_zorder=True：会把 collection（月球
    # plot_surface、scatter、尾迹、卫星）的 zorder 按相机深度**覆盖**，因此原先的
    # scatter(zorder=8) 在近月点位于月球背面视角时会被月球按深度盖住（图标+label
    # 均不可见）。而 Line3D（ax.plot）与 Text3D 的 zorder 不会被覆盖，故改用
    # ax.plot 画标记并给一个高于月球球体的 zorder，保证近/远月点始终显示在月球之前。
    p = d["r_orb"][d["i_peri"]]
    a = d["r_orb"][d["i_apo"]]
    fs3d = L["label_3d_fs"]
    _APSIS_Z = 50  # 高于月球球体（computed 排序约 3~4）与卫星/尾迹（computed）
    # scatter s（pt² 面积）→ plot ms（直径 pt）：ms = 2·√(s/π)
    ms3d = 2.0 * math.sqrt(L["marker_size"] / math.pi)
    ax3d.plot([p[0]], [p[1]], [p[2]], marker="^", ms=ms3d,
              color="#d62728", lw=0, label="近月点", zorder=_APSIS_Z)
    ax3d.plot([a[0]], [a[1]], [a[2]], marker="v", ms=ms3d,
              color="#1f77b4", lw=0, label="远月点", zorder=_APSIS_Z)
    ax3d.text(p[0] + 300, p[1], p[2] + 500, "近月点", color="#d62728",
              fontsize=fs3d, ha="center", weight="bold", zorder=_APSIS_Z)
    ax3d.text(a[0], a[1], a[2] + 1000, "远月点", color="#1f77b4",
              fontsize=fs3d, ha="center", weight="bold", zorder=_APSIS_Z)

    # 动态元素：尾迹 / 卫星 / 速度矢量
    sc_tail = ax3d.scatter([], [], [], s=L["tail_size"], c=[], cmap=cmap,
                           norm=norm, depthshade=False, zorder=5)
    sc_sat = ax3d.scatter([], [], [], s=L["sat_size"], color="#ffd54f",
                          edgecolors="#37474f", linewidths=L["sat_lw"],
                          depthshade=False, zorder=10)
    (line_vec,) = ax3d.plot([], [], [], color="#ff7043", lw=L["vec_lw"],
                            alpha=0.95)
    if L["legend"]:
        ax3d.legend(loc="center left", fontsize=15, framealpha=0.85)

    # ---- 速度–相位曲线（静态曲线 + 动态光标） ------------------------------
    nu_cyc = d["nu_deg"][:N_per]
    v_cyc = d["speeds"][:N_per]
    ax_spd.plot(nu_cyc, v_cyc, color="#546e7a",
                lw=2.2 if is_portrait else 2.0)
    ax_spd.axvline(0.0, color="#d62728", ls="--",
                   lw=1.1 if is_portrait else 1.0, alpha=0.75)
    ax_spd.axvline(180.0, color="#1f77b4", ls="--",
                   lw=1.1 if is_portrait else 1.0, alpha=0.75)
    ax_spd.text(4, d["v_peri"] * 1.05, "近月点(快)", color="#d62728",
                fontsize=8.5 if is_portrait else 12)
    ax_spd.text(184, d["v_apo"] + (d["v_peri"] - d["v_apo"]) * 0.35,
                "远月点(慢)", color="#1f77b4",
                fontsize=8.5 if is_portrait else 12)
    ax_spd.set_xlim(0, 360)
    ax_spd.set_ylim(0, d["v_peri"] * 1.15)
    if is_portrait:
        ax_spd.set_xlabel("真近点角 ν（0° = 近月点）", fontsize=9)
        ax_spd.set_ylabel("速度 (km/s)", fontsize=9)
        ax_spd.tick_params(labelsize=8)
    else:
        ax_spd.set_xlabel("真近点角 ν（0° = 近月点）", fontsize=L['ticks_fs'])
        ax_spd.set_ylabel("速度 (km/s)", fontsize=L['ticks_fs'])
        ax_spd.tick_params(labelsize=L['ticks_fs'])
    ax_spd.set_title("速度 vs 轨道相位", fontsize=18, weight="bold")
    ax_spd.grid(alpha=0.3)
    (cursor_pt,) = ax_spd.plot([], [], "o", color="#ffd54f",
                               ms=10 if is_portrait else 9,
                               mec="#37474f", mew=1.4, zorder=6)
    (cursor_vl,) = ax_spd.plot([], [], color="#ffd54f",
                               lw=1.4 if is_portrait else 1.3, alpha=0.6)

    # ---- 信息面板 ----------------------------------------------------------
    panel = (build_info_portrait(ax_info, d, L) if is_portrait
             else build_info_landscape(ax_info, d, L))

    # ---- 动画更新 ----------------------------------------------------------
    def update(i: int):
        x0, y0, z0 = d["r_orb"][i]
        # 尾迹（等时间间隔 → 间距反映快慢）
        s = max(0, i - args.tail + 1)
        e = i + 1
        sc_tail._offsets3d = (d["r_orb"][s:e, 0], d["r_orb"][s:e, 1],
                              d["r_orb"][s:e, 2])
        sc_tail.set_array(d["speeds"][s:e])
        # 卫星当前位置
        sc_sat._offsets3d = ([x0], [y0], [z0])
        # 速度矢量（长度 ∝ 速率）
        vi = d["v_orb"][i]
        vn = vi / np.linalg.norm(vi)
        alen = max(1200.0, np.linalg.norm(vi) * 2200.0)
        line_vec.set_data_3d(
            [x0, x0 + vn[0] * alen],
            [y0, y0 + vn[1] * alen],
            [z0, z0 + vn[2] * alen],
        )
        # 时间 / 圈次
        t_now = epoch + i * dt_h * u.h
        orbit_idx = i // N_per + 1
        frac = (i % N_per) / N_per
        panel["t_time"].set_text(f"UTC {t_now.isot}")
        panel["t_orbit"].set_text(
            f"第 {orbit_idx}/{args.orbits} 圈 · 圈内相位 {frac * 100:5.1f}%"
        )
        # 高度 / 速度 / 距近远月点（只填数字，靠大字号高亮）
        # 距下一近/远月点 = 距下一次经过的**剩余时间（倒数）**：从周期倒计到 0，
        # 而非距上次经过的已过时间（0→周期递增）。
        d_peri = (N_per - (i - d["i_peri"]) % N_per) % N_per * dt_h
        d_apo = (N_per - (i - d["i_apo"]) % N_per) % N_per * dt_h
        panel["t_h"].set_text(f"{d['alts'][i]:7.1f} km")
        panel["t_v"].set_text(f"{d['speeds'][i]:7.3f} km/s")
        panel["t_peri"].set_text(f"{d_peri:6.2f} h")
        panel["t_apo"].set_text(f"{d_apo:6.2f} h")
        # 速度曲线光标
        nu_i = d["nu_deg"][i]
        cursor_pt.set_data([nu_i], [d["speeds"][i]])
        cursor_vl.set_data([nu_i, nu_i], [0, d["v_peri"] * 1.15])

        if i % 60 == 0:
            print(f"[queqiao_orbit_animation] 帧 {i}/{N} ...")
        return (sc_tail, sc_sat, line_vec, cursor_pt, cursor_vl,
                panel["t_time"], panel["t_orbit"], panel["t_h"], panel["t_v"],
                panel["t_peri"], panel["t_apo"])

    # ---- 预览或渲染 --------------------------------------------------------
    suffix = "" if is_portrait else "_landscape"
    if args.preview_only:
        update(0)
        png_path = os.path.join(args.output_dir,
                                f"queqiao_orbit_preview{suffix}.png")
        fig.savefig(png_path, dpi=args.dpi)
        print(f"[queqiao_orbit_animation] 预览帧 → {png_path}")
        plt.close(fig)
        return

    mp4_path = os.path.join(args.output_dir,
                            f"queqiao_orbit_animation{suffix}.mp4")
    print(f"[queqiao_orbit_animation] 渲染 {N} 帧 @ {args.fps} fps → {mp4_path}")
    anim = FuncAnimation(fig, update, frames=N, interval=1000.0 / args.fps,
                         blit=False)
    # yuv420p 要求宽高均为偶数：figsize@dpi 可能产生奇数像素（如 1705px），
    # 用 scale 滤镜强制取偶，否则 libx264 编码失败（ffmpeg 退出码 187）。
    writer = FFMpegWriter(
        fps=args.fps, codec="libx264", bitrate=-1,
        extra_args=[
            "-pix_fmt", "yuv420p",
            "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",
        ],
    )
    anim.save(mp4_path, writer=writer, dpi=args.dpi)
    plt.close(fig)
    print(f"[queqiao_orbit_animation] DONE → {mp4_path}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR,
                   help="输出目录，默认 figures-and-logs")
    p.add_argument("--epoch", default="2026-08-01 00:00:00",
                   help="轨道历元（UTC），默认 2026-08-01 00:00:00")
    p.add_argument("--layout", choices=["portrait", "landscape"],
                   default="landscape",
                   help="手机画幅：landscape=横屏 16:9（默认）; portrait=竖屏 9:16")
    p.add_argument("--orbits", type=int, default=2, help="动画展示轨道圈数，默认 2")
    p.add_argument("--frames-per-orbit", type=int, default=720,
                   help="每圈采样帧数，默认 720")
    p.add_argument("--tail", type=int, default=30,
                   help="尾迹帧数（等时间间隔，反映快慢），默认 30")
    p.add_argument("--fps", type=int, default=60, help="视频帧率，默认 60")
    p.add_argument("--dpi", type=int, default=120,
                   help="渲染 dpi，默认 120（portrait→1080×1920；landscape→1920×1080）")
    p.add_argument("--zoom", type=float, default=None,
                   help="3D 轨道/月球放大倍率，默认 portrait=1.9 / landscape=1.5；"
                        "远月点弧段会出画，看完整轨道用 1.0")
    p.add_argument("--preview-only", action="store_true",
                   help="只渲染首帧 PNG（用于检查版面），不生成视频")
    return p.parse_args()


if __name__ == "__main__":
    main()
