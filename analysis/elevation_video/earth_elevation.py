#!/usr/bin/env python3
"""P0-2：地球在月球南极的高度角范围。
科学口径：参考面是月球参考椭球（IAU_MOON RADII），不是坑底局部斜面。
算法复用 ``illumination._spice_sun_vector`` 的 ENU 投影几何，仅将目标天体
"SUN" 替换为 "EARTH"。

产出：
- ``earth_elevation_timeseries.json`` — 当前时段时间序列 + 理论极值
- ``earth_elevation_timeseries.png`` — 高度角时间序列图（标注 min/max 线）
- ``earth_elevation_histogram.png`` — 高度角分布直方图（当前范围 vs 理论极值）
- ``earth_sky_view_*.png`` — 2D 天空视角单帧（可选，--sky-view-time）
- ``earth_sky_view_animation.mp4`` — 2D 天空视角延时动画（可选，--sky-animation，
  16:9 手机横屏，每帧 2-4 小时，显示地球在任意纬度/任意时刻天空中的位置）
- ``earth_pov_view.png`` — POV 第一人称地平线视角单帧（可选，--pov-view-time）
- ``earth_pov_animation.mp4`` — POV 地平线视角 4K 延时动画（可选，--pov-animation，
  站在月面平坦地面看地球升起/落下）

天空视角（等距鱼眼全景，适合移动端横屏查看）：
- 黑色圆盘 = 宇宙 / 天空（天顶在圆心，地平线在圆环）
- 圆盘外灰色 = 月面（替换 DEM 地形渲染）
- 蓝色球体 = 地球（按 SPICE 计算的高度角/方位角/视径定位，带淡蓝轨迹）
- 顶部/底部四角 = 简洁演示文字（UTC 时间、高度角、视直径、方位角、坐标）
- 画面四周留白 ≥10%-15% 安全区（--sky-safe-frac 可调）

POV 第一人称地平线视角（针孔透视，站在月面平坦地面）：
- 黑色天空 + 星空，灰色月面（近处更深 + 巨石透视），蓝色地球
- 相机朝向地球方位角（--pov-az 可手动指定，默认自动指向）
- 地球在地平线以下时被月面自然遮挡——可见完整的升起/落下
  （月球潮汐锁定 → 地球仅在极区/边缘区接近地平线，周期 ~27.3 天）
- 右下角显示动画时段内的升起/落下时刻（↑/↓，UTC）
- 四角 HUD 距画布边 ≥10% 安全区（--pov-tilt 调俯仰、--pov-fov 调视场）

用法::

    # 原有功能：高度角时间序列 + 直方图 + JSON
    python analysis/elevation_video/earth_elevation.py \\
        --lat -89.67 --lon 129.78 --h 0.0 \\
        --start 2026-08-01 --end 2027-12-31 --interval-h 1.0 \\
        --theoretical-extent-years 18.6 --output-dir figures-and-logs

    # 2D 天空视角单帧（任意纬度/时刻）
    python analysis/elevation_video/earth_elevation.py \\
        --lat -89.67 --lon 129.78 --sky-view-time 2026-08-15T12:00:00 \\
        --sky-dpi 120 --output-dir figures-and-logs

    # 2D 天空 4K 延时动画（一帧 3 小时，3840×2160 手机横屏）
    python analysis/elevation_video/earth_elevation.py \\
        --lat -89.67 --lon 129.78 \\
        --sky-animation --sky-start 2026-08-01 --sky-end 2026-08-31 \\
        --sky-interval-h 3.0 --sky-fps 24 --sky-dpi 240 \\
        --output-dir figures-and-logs

    # POV 第一人称单帧（站在月面看地球）
    python analysis/elevation_video/earth_elevation.py \\
        --lat -89.67 --lon 129.78 --pov-view-time 2026-08-05T12:00:00 \\
        --pov-dpi 120 --output-dir figures-and-logs

    # POV 4K 延时动画（地球升起/落下，一帧 3 小时）
    python analysis/elevation_video/earth_elevation.py \\
        --lat -89.67 --lon 129.78 \\
        --pov-animation --pov-start 2026-08-01 --pov-end 2026-08-31 \\
        --pov-interval-h 3.0 --pov-fps 24 --pov-dpi 240 \\
        --output-dir figures-and-logs
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from datetime import datetime, timedelta, timezone

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter, FuncAnimation
from matplotlib.patches import FancyBboxPatch

# 允许以脚本方式从仓库任意位置运行：将仓库根加入 sys.path
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from analysis.elevation_video.utils import (  # noqa: E402
    SHACKLETON_LAT_DEG,
    SHACKLETON_LON_DEG,
    cleanup_spice,
    enu_elevation_azimuth,
    load_spice_kernels,
    mjd_to_utc_iso,
    moon_radii_m,
    utc_iso_to_mjd,
    make_time_grid,
    setup_cjk_font,
)

# 地球赤道半径（m，视大小计算用，IAU 2009 推荐值）
EARTH_R_EQ_M: float = 6378137.0

DEFAULT_OUTPUT_DIR = os.path.join(_REPO_ROOT, "logs", "elevation_video")

# 理论极值时间范围（完整天平动周期 ~18.6 年）
THEORETICAL_START = "2008-01-01T00:00:00"
THEORETICAL_END = "2026-12-31T23:00:00"


def earth_elevation_angle(
    lon_deg: float,
    lat_deg: float,
    h_m: float,
    t_mjd: float,
    R_eq: float,
    R_pol: float,
) -> float:
    """计算地球在观测点的当地水平坐标系下的高度角（度）。

    利用 ``utils.enu_elevation_angle`` 复用共享 ENU 投影逻辑，
    仅替换目标天体为 EARTH。

    :param lon_deg, lat_deg: 观测点经纬度（度）
    :param h_m: 观测点高程（米，参考椭球面以上）
    :param t_mjd: 时间（MJD）
    :param R_eq: 月球赤道半径（米，从 PCK 读取）
    :param R_pol: 月球极半径（米，从 PCK 读取）
    """
    import spiceypy

    jd = t_mjd + 2400000.5
    et = spiceypy.unitim(jd, "JDTDB", "ET")
    state, _lt = spiceypy.spkezr("EARTH", et, "IAU_MOON", "LT+S", "MOON")
    earth_pos_m = np.asarray(state[:3], dtype=float) * 1000.0  # km -> m
    # 复用共享 ENU 高度角（参考椭球面法线为 Up）
    from analysis.elevation_video.utils import enu_elevation_angle

    return enu_elevation_angle(earth_pos_m, lon_deg, lat_deg, h_m, R_eq, R_pol)


def earth_az_el(
    lon_deg: float,
    lat_deg: float,
    h_m: float,
    t_mjd: float,
    R_eq: float,
    R_pol: float,
) -> tuple[float, float]:
    """地球在当地水平坐标系下的 (方位角, 高度角)，单位度。

    与 :func:`earth_elevation_angle` 使用相同的 SPICE 观测几何（EARTH 在
    IAU_MOON 月固系），额外返回方位角：0°=北、90°=东、180°=南、270°=西
    （地平线图约定，与 utils.enu_elevation_azimuth 一致）。

    :return: (azimuth_deg, elevation_deg)
    """
    import spiceypy

    jd = t_mjd + 2400000.5
    et = spiceypy.unitim(jd, "JDTDB", "ET")
    state, _lt = spiceypy.spkezr("EARTH", et, "IAU_MOON", "LT+S", "MOON")
    earth_pos_m = np.asarray(state[:3], dtype=float) * 1000.0  # km -> m
    el, az = enu_elevation_azimuth(earth_pos_m, lon_deg, lat_deg, h_m, R_eq, R_pol)
    return az, el


def earth_apparent_radius_deg(
    lon_deg: float,
    lat_deg: float,
    h_m: float,
    t_mjd: float,
    R_eq: float,
    R_pol: float,
) -> float:
    """地球视半径（度）——用于 2D 天空图中地球球体的绘制大小。

    计算观测点到地球中心的距离，再按地球赤道半径得到视半径
    （sin⁻¹(R_E / d)，球面近似，地平线附近会有轻微折射差异可忽略）。
    月球表面典型值：视半径 ≈0.9°–1.0°（视直径 ≈1.9°–2.0°，
    地月距离 ~38–40 万 km）。
    """
    import spiceypy

    jd = t_mjd + 2400000.5
    et = spiceypy.unitim(jd, "JDTDB", "ET")
    state, _lt = spiceypy.spkezr("EARTH", et, "IAU_MOON", "LT+S", "MOON")
    earth_pos_m = np.asarray(state[:3], dtype=float) * 1000.0  # km -> m

    lon_rad = math.radians(lon_deg)
    lat_rad = math.radians(lat_deg)
    radius = R_eq + h_m
    obs = np.array(
        [
            radius * math.cos(lat_rad) * math.cos(lon_rad),
            radius * math.cos(lat_rad) * math.sin(lon_rad),
            radius * math.sin(lat_rad),
        ]
    )
    dist_m = float(np.linalg.norm(earth_pos_m - obs))
    return math.degrees(math.asin(min(1.0, EARTH_R_EQ_M / dist_m)))


# ---------------------------------------------------------------------------
# 2D 等距鱼眼天空视角（任意纬度 / 任意时刻）
# ---------------------------------------------------------------------------
# 等距投影：径向距离 ∝ 天顶角（高度角 h → 半径 r = R·(90°-h)/90°）。
# 方位角 az（北=0°、东=90°）→ 极角 θ = radians(90°-az)，北在顶部。

# 画面配色（简洁演示风：黑色宇宙 / 灰色月面 / 蓝色地球）
_SKY_COLORS = {
    "space": "#05070d",          # 宇宙深黑
    "space_ring": "#111827",     # 天空圆盘边缘
    "moon": "#8a8f98",           # 月面灰
    "moon_light": "#b8bcc4",     # 月面纹理亮带
    "horizon": "#d8dce2",        # 地平线环
    "grid": "#5a6472",           # 高度圈 / 方位线
    "grid_label": "#c8d0da",     # 网格标注
    "earth": "#3b82f6",          # 地球蓝
    "earth_rim": "#93c5fd",      # 地球亮边
    "earth_track": "#60a5fa",    # 地球轨迹
    "text": "#f1f5f9",           # 演示文字
    "text_dim": "#a8b3c2",       # 次级文字
    "card": "#16212e",           # 信息卡底色
    "card_edge": "#2e4057",      # 信息卡描边
}

# 灰色月面（圆盘外）上的简单陨石坑纹理：(x, y, 半径)，画面单位
_MOON_CRATERS: tuple[tuple[float, float, float], ...] = (
    (1.8, 6.8, 0.55), (14.0, 7.0, 0.45), (1.5, 2.0, 0.62), (14.2, 2.2, 0.50),
    (2.6, 4.5, 0.32), (13.4, 4.5, 0.38), (5.9, 8.2, 0.30), (10.1, 8.4, 0.28),
    (5.6, 0.7, 0.26), (10.4, 0.65, 0.30),
)


def _fisheye_xy(az_deg: float, el_deg: float, R: float):
    """等距鱼眼投影：(方位角, 高度角) → 圆盘内 (x, y)（圆心=天顶）。"""
    theta = math.radians(90.0 - az_deg)
    r = R * (90.0 - el_deg) / 90.0
    return r * math.cos(theta), r * math.sin(theta)


def _card(ax, x, y, w, h, title, value, sub=None, fs_title=13, fs_value=17):
    """半透明信息卡（角落，transAxes 坐标）。"""
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.008",
        transform=ax.transAxes, facecolor=_SKY_COLORS["card"],
        edgecolor=_SKY_COLORS["card_edge"], linewidth=1.1,
        alpha=0.82, zorder=8,
    ))
    ax.text(x + 0.028, y + h - 0.028, title, transform=ax.transAxes,
            fontsize=fs_title, color=_SKY_COLORS["text_dim"],
            ha="left", va="top", zorder=9)
    ax.text(x + 0.028, y + 0.028, value, transform=ax.transAxes,
            fontsize=fs_value, color=_SKY_COLORS["text"],
            ha="left", va="bottom", weight="bold", zorder=9)
    if sub:
        ax.text(x + w - 0.028, y + 0.028, sub, transform=ax.transAxes,
                fontsize=fs_title, color=_SKY_COLORS["text_dim"],
                ha="right", va="bottom", zorder=9)


def draw_sky_view(
    ax,
    t_mjd: float,
    lat_deg: float,
    lon_deg: float,
    h_m: float,
    R_eq: float,
    R_pol: float,
    sky_scale: float = 1.0,
    trail_hours: float = 48.0,
    show_track: bool = True,
):
    """在给定 axes 上绘制一帧 2D 等距鱼眼天空视角。

    :param ax: matplotlib Axes（已配置为 16:9 等比例圆盘画面）
    :param t_mjd: 当前时刻（MJD）
    :param lat_deg, lon_deg, h_m: 观测点
    :param R_eq, R_pol: 月球椭球半径
    :param sky_scale: 天空圆盘半径倍率（默认 1.0 = 占画面高度约 45%）
    :param trail_hours: 地球轨迹回看时长（小时），默认 48
    :param show_track: 是否绘制地球历史轨迹
    """
    az_now, el_now = earth_az_el(lon_deg, lat_deg, h_m, t_mjd, R_eq, R_pol)
    app_r = earth_apparent_radius_deg(lon_deg, lat_deg, h_m, t_mjd, R_eq, R_pol)

    # ---- 画布（16:9）----
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 9)
    ax.set_aspect("equal")
    ax.axis("off")

    # ---- 灰色月面（圆盘外）----
    ax.set_facecolor(_SKY_COLORS["moon"])
    fig = ax.figure
    fig.patch.set_facecolor(_SKY_COLORS["moon"])

    # 圆盘圆心/半径（等距鱼眼全景，四周留 ≥10%-15% 安全区）
    # 画布 16×9：圆盘直径 ≈ 短边 78%（R=3.5），上下边距各 ≈11%（1.0/9），
    # 左右边距更宽（4.5/9 ≈ 50%），动画内容不会贴边。
    cx, cy, R = 8.0, 4.5, 3.5 * sky_scale
    margin = min(16.0 - (cx + R), cx - R, 9.0 - (cy + R), cy - R)
    safe_frac = margin / 9.0
    assert safe_frac >= 0.10, f"sky disk margin {safe_frac:.1%} < 10% safe area"

    # ---- 月面陨石坑纹理（灰色月面细节，被圆盘覆盖的部分自动隐藏）----
    for crx, cry, crr in _MOON_CRATERS:
        if (crx - cx) ** 2 + (cry - cy) ** 2 > (R + crr) ** 2:
            ax.add_patch(plt.Circle(
                (crx, cry), crr,
                facecolor=_SKY_COLORS["moon_light"], alpha=0.45,
                edgecolor=_SKY_COLORS["moon"], linewidth=1.0, zorder=1,
            ))

    # ---- 黑色宇宙（天空圆盘）----
    sky = plt.Circle((cx, cy), R, facecolor=_SKY_COLORS["space"],
                     edgecolor=_SKY_COLORS["space_ring"], linewidth=2.5, zorder=2)
    ax.add_patch(sky)

    # ---- 高度圈（30°/60°）与方位线（N/E/S/W）----
    for elv, ls, lw in ((60.0, (0, (4, 4)), 1.1), (30.0, (0, (4, 4)), 1.1)):
        th = np.linspace(0, 2 * np.pi, 361)
        rr = R * (90.0 - elv) / 90.0
        ax.plot(cx + rr * np.cos(th), cy + rr * np.sin(th),
                color=_SKY_COLORS["grid"], ls=ls, lw=lw, zorder=3)
    # 方位线（北/东/南/西）
    for azd in (0, 90, 180, 270):
        x, y = _fisheye_xy(azd, 0.0, R)
        ax.plot([cx, cx + x], [cy, cy + y], color=_SKY_COLORS["grid"],
                lw=0.9, ls=":", alpha=0.7, zorder=3)

    # ---- 高度圈刻度（30°/60°）----
    for elv in (30, 60):
        x, y = _fisheye_xy(0, elv, R)  # 北向刻度
        ax.text(cx + x, cy + y + 0.09, f"{elv}°", color=_SKY_COLORS["grid_label"],
                fontsize=8.5, ha="center", va="bottom", zorder=4)
    # 方位角标签
    for azd, lab, dx, dy in ((0, "北", 0, 0.35), (90, "东", 0.35, 0),
                             (180, "南", 0, -0.42), (270, "西", -0.42, 0)):
        x, y = _fisheye_xy(azd, 0.0, R)
        ax.text(cx + x + dx, cy + y + dy, lab, color=_SKY_COLORS["grid_label"],
                fontsize=9.5, ha="center", va="center", zorder=4)

    # ---- 地平线环 ----
    th = np.linspace(0, 2 * np.pi, 361)
    ax.plot(cx + R * np.cos(th), cy + R * np.sin(th),
            color=_SKY_COLORS["horizon"], lw=2.2, zorder=4)
    ax.text(cx, cy - R - 0.38, "地平线 0°", color=_SKY_COLORS["horizon"],
            fontsize=9, ha="center", va="top", zorder=4)

    # ---- 地球历史轨迹（淡蓝虚线）----
    if show_track and trail_hours > 0:
        trail_mjd = t_mjd - np.linspace(0, trail_hours / 24.0, 60)[::-1]
        tr_x, tr_y = [], []
        for tm in trail_mjd:
            az_t, el_t = earth_az_el(lon_deg, lat_deg, h_m, tm, R_eq, R_pol)
            if el_t >= -89.0:
                x, y = _fisheye_xy(az_t, el_t, R)
                tr_x.append(cx + x)
                tr_y.append(cy + y)
        if tr_x:
            ax.plot(tr_x, tr_y, color=_SKY_COLORS["earth_track"], lw=1.8,
                    ls="--", alpha=0.45, zorder=5)

    # ---- 地球（蓝色球体，按高度角/方位角定位）----
    if el_now >= -89.0:
        ex, ey = _fisheye_xy(az_now, el_now, R)
        earth_diam = 2.0 * app_r * R / 90.0
        # 视径下限（0.3 画面单位），保证低分辨率下可见
        earth_diam = max(earth_diam, 0.42)
        earth = plt.Circle((cx + ex, cy + ey), earth_diam / 2.0,
                           facecolor=_SKY_COLORS["earth"],
                           edgecolor=_SKY_COLORS["earth_rim"], linewidth=1.4,
                           zorder=7)
        ax.add_patch(earth)
        # 高度角引线（当前时刻）
        hx, hy = _fisheye_xy(az_now, 0.0, R)
        ax.plot([cx + hx, cx + ex], [cy + hy, cy + ey],
                color=_SKY_COLORS["earth_rim"], lw=1.0, alpha=0.65, zorder=6)

    return {"az": az_now, "el": el_now, "apparent_radius_deg": app_r}


def render_sky_view_frame(
    fig,
    ax,
    t_mjd: float,
    lat_deg: float,
    lon_deg: float,
    h_m: float,
    R_eq: float,
    R_pol: float,
    trail_hours: float = 48.0,
    show_track: bool = True,
) -> dict:
    """渲染一帧天空视角（清空 axes 后调用 draw_sky_view + 文字信息卡）。

    :return: {"az", "el", "apparent_radius_deg"}（供动画更新复用）
    """
    ax.clear()
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 9)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_facecolor(_SKY_COLORS["moon"])
    fig.patch.set_facecolor(_SKY_COLORS["moon"])

    info = draw_sky_view(ax, t_mjd, lat_deg, lon_deg, h_m, R_eq, R_pol,
                         trail_hours=trail_hours, show_track=show_track)
    az, el = info["az"], info["el"]
    app_r = info["apparent_radius_deg"]

    # ---- 简洁演示文字（四角信息卡 + 顶部标题）----
    # 时间（UTC）
    dt = datetime(1858, 11, 17, tzinfo=timezone.utc) + timedelta(days=float(t_mjd))
    time_str = dt.strftime("%Y-%m-%d %H:%M UTC")
    _card(ax, 0.10, 0.6, 0.30, 0.16, "时间 · UTC",
          time_str, sub=f"Δt={trail_hours:.0f} h 轨迹", fs_value=14)
    # 高度角
    _card(ax, 0.62, 0.6, 0.30, 0.16, "地球高度角",
          f"{el:+.2f}°", sub=f"视直径 {2*app_r:.2f}°", fs_value=18)
    # 方位角
    _card(ax, 0.10, 0.4, 0.30, 0.16, "地球方位角",
          f"{az:5.1f}°", sub="北=0°", fs_value=18)
    # 观测点
    _card(ax, 0.62, 0.4, 0.30, 0.16, "观测点",
          f"{lat_deg:.2f}°N  {lon_deg:.2f}°E", sub="月球表面", fs_value=14)
    # 顶部标题（居中）
    ax.text(0.5, 0.955, "月球天空 · 地球视角", transform=ax.transAxes,
            fontsize=17, color=_SKY_COLORS["text"], ha="center", va="top",
            weight="bold", zorder=9)
    ax.text(0.5, 0.905, "等距鱼眼全景 · 天顶在圆心 · 地平线在圆环",
            transform=ax.transAxes, fontsize=10.5, color=_SKY_COLORS["text_dim"],
            ha="center", va="top", zorder=9)
    return info


def make_sky_view_figure(figsize=(16.0, 9.0), dpi=120):
    """创建天空视角 Figure（16:9 手机横屏，深色画布）。"""
    fig = plt.figure(figsize=figsize, dpi=dpi)
    fig.patch.set_facecolor(_SKY_COLORS["moon"])
    ax = fig.add_axes([0, 0, 1, 1])
    return fig, ax


def render_sky_view_single(
    t_mjd: float,
    lat_deg: float,
    lon_deg: float,
    h_m: float,
    R_eq: float,
    R_pol: float,
    output_path: str,
    dpi: int = 120,
    show_track: bool = True,
) -> str:
    """渲染单帧 2D 天空视角 PNG（任意纬度/任意时刻）。"""
    setup_cjk_font()
    fig, ax = make_sky_view_figure((16.0, 9.0), dpi)
    render_sky_view_frame(fig, ax, t_mjd, lat_deg, lon_deg, h_m, R_eq, R_pol,
                          show_track=show_track)
    fig.savefig(output_path, dpi=dpi, facecolor=fig.get_facecolor())
    plt.close(fig)
    return output_path


def render_sky_view_animation(
    t_mjd_grid: np.ndarray,
    lat_deg: float,
    lon_deg: float,
    h_m: float,
    R_eq: float,
    R_pol: float,
    output_path: str,
    fps: int = 24,
    dpi: int = 240,
    trail_hours: float = 48.0,
    show_track: bool = True,
    preview_only: bool = False,
) -> str:
    """渲染 2D 天空视角 4K 延时动画（16:9 手机横屏，MP4）。

    每帧时间步长由 ``t_mjd_grid`` 相邻点决定（``--sky-interval-h``，默认 3 小时，
    即“一帧 2-4 小时”的延时节奏）。使用 FuncAnimation + FFMpegWriter（与
    queqiao_orbit_animation.py 相同的 yuv420p + 偶数宽高处理，libx264）。
    ``dpi=240 @ 16×9`` → 3840×2160（4K UHD）。

    :param t_mjd_grid: 动画帧时刻（MJD，升序均匀网格）
    :param preview_only: 只渲染首帧 PNG（版面检查），不编码视频
    :return: 输出的 PNG（preview）或 MP4 路径
    """
    setup_cjk_font()
    fig, ax = make_sky_view_figure((16.0, 9.0), dpi)

    def update(i: int):
        render_sky_view_frame(
            fig, ax, float(t_mjd_grid[i]), lat_deg, lon_deg, h_m, R_eq, R_pol,
            trail_hours=trail_hours, show_track=show_track,
        )
        return ()

    N = len(t_mjd_grid)
    if preview_only:
        update(0)
        png_path = os.path.splitext(output_path)[0] + "_preview.png"
        fig.savefig(png_path, dpi=dpi, facecolor=fig.get_facecolor())
        plt.close(fig)
        print(f"[earth_elevation] 动画首帧预览 → {png_path}")
        return png_path

    print(f"[earth_elevation] 渲染天空视角动画 {N} 帧 @ {fps} fps → {output_path}")
    anim = FuncAnimation(fig, update, frames=N, interval=1000.0 / fps, blit=False)
    # yuv420p 要求宽高均为偶数：用 scale 滤镜强制取偶，避免 libx264 编码失败
    writer = FFMpegWriter(
        fps=fps, codec="libx264", bitrate=-1,
        extra_args=[
            "-pix_fmt", "yuv420p",
            "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",
        ],
    )
    anim.save(output_path, writer=writer, dpi=dpi)
    plt.close(fig)
    print(f"[earth_elevation] DONE → {output_path}")
    return output_path


# ---------------------------------------------------------------------------
# POV 第一人称地平线视角（站在月面平坦地面，看向地球，可见升起/落下）
# ---------------------------------------------------------------------------
# 相机模型：针孔透视。yaw=view_az（相机朝向方位角），pitch=tilt（俯仰角）。
# 屏幕：x∈[0,1] 归一化横坐标（绘制时乘画面宽高比 A=16/9），y∈[0,1] 纵坐标
# （0=画面底，1=画面顶，地平线 el=0 在 y=0.5+tilt 修正处）。
# 月面为平坦地面（观测者眼高 eye_h），地平线是屏幕上的水平直线。

_POV_COLORS = {
    "space": "#04060b",          # 宇宙深黑
    "star": "#c9d4e0",           # 星空
    "ground": "#5c6169",         # 月面灰（近处更深）
    "ground_near": "#464b53",    # 近处月面
    "horizon_glow": "#8b929b",   # 地平线微光带
    "horizon": "#d8dce2",        # 地平线
    "boulder": "#7b8088",        # 月面巨石
    "earth": "#3b82f6",          # 地球蓝
    "earth_rim": "#93c5fd",      # 地球亮边
    "earth_night": "#1e3a8a",    # 地球夜面
    "earth_track": "#60a5fa",    # 地球轨迹
    "guide": "#93c5fd",          # 高度引线
    "text": "#f1f5f9",
    "text_dim": "#a8b3c2",
    "card": "#16212e",
    "card_edge": "#2e4057",
}
POV_ASPECT = 16.0 / 9.0         # 手机横屏


def _pov_camera_basis(view_az_deg: float, tilt_deg: float):
    """相机基 (F, R, U)（ENU 东/北/上坐标）。yaw=view_az, pitch=tilt。

    F=视向单位向量，R=视向右手单位向量，U=画面上方向量（正交右手系）。
    """
    yaw = math.radians(view_az_deg)
    pit = math.radians(tilt_deg)
    F = np.array([
        math.sin(yaw) * math.cos(pit),
        math.cos(yaw) * math.cos(pit),
        math.sin(pit),
    ])
    R = np.array([math.cos(yaw), -math.sin(yaw), 0.0])
    U = np.array([
        -math.sin(yaw) * math.sin(pit),
        -math.cos(yaw) * math.sin(pit),
        math.cos(pit),
    ])
    return F, R, U


def _pov_dir(az_deg: float, el_deg: float) -> np.ndarray:
    """方位角/高度角 → 单位方向向量（ENU）。"""
    az = math.radians(az_deg)
    el = math.radians(el_deg)
    return np.array([
        math.cos(el) * math.sin(az),
        math.cos(el) * math.cos(az),
        math.sin(el),
    ])


def _pov_project(az_deg, el_deg, F, R, U, fov_x_deg, fov_y_deg):
    """(方位角, 高度角) → 屏幕归一化坐标 (x, y)；视场外返回 None。

    x∈[0,1]（画图时乘 POV_ASPECT），y∈[0,1]（0=底 1=顶）。
    """
    D = _pov_dir(az_deg, el_deg)
    dx = float(D @ R)
    dy = float(D @ U)
    dz = float(D @ F)
    if dz <= 1e-6:
        return None
    kx = 0.5 / math.tan(math.radians(fov_x_deg) / 2.0)
    ky = 0.5 / math.tan(math.radians(fov_y_deg) / 2.0)
    x = 0.5 + kx * dx / dz
    # U 朝上：dy>0（目标在相机上方）→ y>0.5（数据坐标更高，屏幕上更靠上）
    y = 0.5 + ky * dy / dz
    if not (-0.25 <= x <= 1.25 and -0.25 <= y <= 1.25):
        return None
    return x, y


def _pov_project_ground(e, n_, eye_h, F, R, U, fov_x_deg, fov_y_deg):
    """月面地面点 (e, n)（米，ENU 水平）→ (x, y, 距离)；画布外返回 None。"""
    D = np.array([e, n_, -eye_h], dtype=float)
    dist = float(np.linalg.norm(D))
    D /= dist
    dx = float(D @ R)
    dy = float(D @ U)
    dz = float(D @ F)
    if dz <= 1e-6:
        return None
    kx = 0.5 / math.tan(math.radians(fov_x_deg) / 2.0)
    ky = 0.5 / math.tan(math.radians(fov_y_deg) / 2.0)
    x = 0.5 + kx * dx / dz
    # 地面点在眼高以下（dy<0）→ y<0.5 → 屏幕上更低（在地平线下方）
    y = 0.5 + ky * dy / dz
    if not (-0.1 <= x <= 1.1 and -0.1 <= y <= 1.1):
        return None
    return x, y, dist


def _pov_fov_y(fov_x_deg: float, aspect: float = POV_ASPECT) -> float:
    """由水平视场角与画面宽高比推导垂直视场角（度）。"""
    return 2.0 * math.degrees(math.atan(math.tan(math.radians(fov_x_deg) / 2.0) / aspect))


def _pov_stars(rng, view_az, tilt_deg, fov_x_deg, fov_y_deg, n=220):
    """生成天空星点（仅地平线以上，屏幕坐标 + 点大小）。"""
    F, R, U = _pov_camera_basis(view_az, tilt_deg)
    max_el = tilt_deg + fov_y_deg / 2.0
    stars = []
    for _ in range(n * 4):
        az = view_az + rng.uniform(-fov_x_deg * 0.62, fov_x_deg * 0.62)
        el = rng.uniform(0.0, max_el * 0.96)
        p = _pov_project(az, el, F, R, U, fov_x_deg, fov_y_deg)
        if p is not None and p[1] < 1.0:
            stars.append((p[0], p[1], rng.uniform(0.7, 1.9)))
        if len(stars) >= n:
            break
    return stars


def _pov_boulders(rng, view_az, fov_x_deg, n=42, eye_h=1.7):
    """生成月面巨石（ENU 位置 + 半径 + 距离），透视深度感。"""
    azs = view_az + rng.uniform(-fov_x_deg * 0.68, fov_x_deg * 0.68, n)
    dists = 10.0 ** rng.uniform(1.8, 3.8, n)   # ~63 m … ~6300 m
    radii = rng.uniform(0.4, 2.4, n)
    items = []
    for a, d, r in zip(azs, dists, radii):
        ar = math.radians(a)
        items.append((d * math.sin(ar), d * math.cos(ar), r, d))
    return items


def draw_pov_view(
    ax,
    t_mjd: float,
    lat_deg: float,
    lon_deg: float,
    h_m: float,
    R_eq: float,
    R_pol: float,
    view_az: float,
    tilt_deg: float,
    fov_x: float = 100.0,
    eye_h: float = 1.7,
    trail_hours: float = 48.0,
    show_track: bool = True,
    star_seed: int = 42,
    boulder_seed: int = 7,
):
    """在给定 axes 上绘制一帧 POV 地平线视角。

    布局（16:9 数据坐标，x∈[0,POV_ASPECT] y∈[0,1]，equal aspect）：
    - 上半/地平线以上 = 黑色宇宙 + 星空
    - 地平线以下 = 灰色月面（近处更深）+ 巨石透视
    - 蓝色地球 = 按 SPICE 高度角/方位角/视径定位（+ 历史轨迹 + 引线）

    :param view_az: 相机朝向方位角（度，北=0）
    :param tilt_deg: 相机俯仰角（度，>0 向上看），决定地平线在画面中的位置
    :param fov_x: 水平视场角（度）
    :param eye_h: 眼高（米）
    """
    az_now, el_now = earth_az_el(lon_deg, lat_deg, h_m, t_mjd, R_eq, R_pol)
    app_r = earth_apparent_radius_deg(lon_deg, lat_deg, h_m, t_mjd, R_eq, R_pol)

    # ---- 画布（equal aspect，x 乘宽高比保持 16:9）----
    A = POV_ASPECT
    ax.set_xlim(0, A)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_facecolor(_POV_COLORS["space"])
    fig = ax.figure
    fig.patch.set_facecolor(_POV_COLORS["space"])

    fov_y = _pov_fov_y(fov_x, A)
    F, R, U = _pov_camera_basis(view_az, tilt_deg)
    ky = 0.5 / math.tan(math.radians(fov_y) / 2.0)
    # 地平线 = 相机朝向方位角处 el=0 的投影 y（与 _pov_project 自洽）
    _hp = _pov_project(view_az, 0.0, F, R, U, fov_x, fov_y)
    y_horizon = _hp[1] if _hp is not None else 0.5

    # ---- 星空（地平线以上，相机固定 → 星点静态）----
    rng_star = np.random.default_rng(star_seed)
    for sx, sy, ss in _pov_stars(rng_star, view_az, tilt_deg, fov_x, fov_y):
        ax.plot(A * sx, sy, ".", color=_POV_COLORS["star"], ms=ss, alpha=0.85,
                zorder=1)

    # ---- 灰色月面（地平线以下）+ 近处加深 + 地平线微光 ----
    # zorder 4 > 地球 zorder 3：地平线以下的地球部分被月面遮挡（升起/落下效果）
    yh = max(-0.05, min(1.05, y_horizon))
    ax.fill([0, A, A, 0], [yh, yh, -0.2, -0.2],
            color=_POV_COLORS["ground"], zorder=4)
    if yh > 0.0:
        ax.fill([0, A, A, 0], [min(0.0, yh), min(0.0, yh), 0.22, 0.22],
                color=_POV_COLORS["ground_near"], zorder=4)
        # 地平线上方微光带（真空辉光示意）
        ax.fill([0, A, A, 0], [yh, yh, min(1.0, yh + 0.03), min(1.0, yh + 0.03)],
                color=_POV_COLORS["horizon_glow"], alpha=0.30, zorder=4)
    # 地平线
    ax.plot([0, A], [yh, yh], color=_POV_COLORS["horizon"], lw=2.4, zorder=6)

    # ---- 月面巨石（透视：近大远小、近低远高贴地平线）----
    rng_bo = np.random.default_rng(boulder_seed)
    for be, bn, br, bd in _pov_boulders(rng_bo, view_az, fov_x, eye_h=eye_h):
        pb = _pov_project_ground(be, bn, eye_h, F, R, U, fov_x, fov_y)
        if pb is None or pb[1] > yh - 1e-4:
            continue
        bx, by, bdist = pb
        b_r = max(ky * br / bdist, 0.0015)   # 屏幕半径（y 单位）
        ax.add_patch(plt.Circle((A * bx, by), b_r,
                                facecolor=_POV_COLORS["boulder"],
                                edgecolor=_POV_COLORS["ground"],
                                linewidth=0.5, zorder=3))

    # ---- 地球历史轨迹（淡蓝虚线）----
    # 先画（zorder 2 < 月面 zorder 4），地平线以下部分被月面自然遮挡
    if show_track and trail_hours > 0:
        trail_mjd = t_mjd - np.linspace(0, trail_hours / 24.0, 28)[::-1]
        tr_x, tr_y = [], []
        for tm in trail_mjd:
            az_t, el_t = earth_az_el(lon_deg, lat_deg, h_m, tm, R_eq, R_pol)
            p = _pov_project(az_t, el_t, F, R, U, fov_x, fov_y)
            if p is not None:
                tr_x.append(A * p[0])
                tr_y.append(p[1])
        if tr_x:
            ax.plot(tr_x, tr_y, color=_POV_COLORS["earth_track"], lw=1.8,
                    ls="--", alpha=0.5, zorder=2)

    # ---- 地球（蓝色球体 / 地平线以下用虚线球体）----
    # 地平线以上：实心蓝色球体（zorder 3 < 月面 zorder 4，升起/落下时
    #   只有露出地平线的部分可见，真实 POV 效果）。
    # 地平线以下：半透明虚线球体（zorder 9 > 月面 zorder 4，叠加在地面上，
    #   表示地球此刻的实际位置——已在地平线之下）。
    p_az = _pov_project(az_now, 0.0, F, R, U, fov_x, fov_y)
    p_earth = _pov_project(az_now, el_now, F, R, U, fov_x, fov_y)
    earth_drawn = False
    if p_earth is not None and -0.15 <= p_earth[1] <= 1.15:
        ex, ey = A * p_earth[0], p_earth[1]
        earth_r = max(ky * math.tan(math.radians(app_r)), 0.011)
        if el_now < 0.0:
            # 虚线球体：浅色半透明填充 + 蓝色虚线圆（叠加在地面上方）
            ax.add_patch(plt.Circle((ex, ey), earth_r,
                                    facecolor=_POV_COLORS["earth"],
                                    edgecolor=_POV_COLORS["earth_rim"],
                                    linewidth=1.6, linestyle="--",
                                    alpha=0.42, zorder=9))
            # 内圈虚线（球体感）
            ax.add_patch(plt.Circle((ex, ey), earth_r * 0.62,
                                    facecolor="none",
                                    edgecolor=_POV_COLORS["earth_rim"],
                                    linewidth=0.9, linestyle="--",
                                    alpha=0.32, zorder=9))
        else:
            ax.add_patch(plt.Circle((ex, ey), earth_r,
                                    facecolor=_POV_COLORS["earth"],
                                    edgecolor=_POV_COLORS["earth_rim"],
                                    linewidth=1.6, zorder=3))
            # 夜面（假影：暗蓝新月，立体感）
            off = 0.45 * earth_r
            ax.add_patch(plt.Circle((ex - off, ey - off), earth_r * 0.86,
                                    facecolor=_POV_COLORS["earth_night"],
                                    alpha=0.85, zorder=3))
            # 中央亮斑（大气辉光）
            ax.add_patch(plt.Circle((ex, ey), earth_r * 0.32,
                                    facecolor="#bfdbfe", alpha=0.5, zorder=3))
        earth_drawn = True

    # ---- 高度引线（地球方位角处的垂线，从地平线到地球，画在月面上）----
    if p_az is not None and p_earth is not None and p_earth[1] >= yh - 0.02:
        ax.plot([A * p_az[0], A * p_earth[0]], [yh, p_earth[1]],
                color=_POV_COLORS["guide"], lw=1.0, alpha=0.55,
                ls=":", zorder=7)

    return {
        "az": az_now, "el": el_now, "apparent_radius_deg": app_r,
        "y_horizon": yh, "fov_y": fov_y, "earth_drawn": earth_drawn,
    }


def render_pov_frame(
    fig,
    ax,
    t_mjd: float,
    lat_deg: float,
    lon_deg: float,
    h_m: float,
    R_eq: float,
    R_pol: float,
    view_az: float,
    tilt_deg: float,
    fov_x: float = 100.0,
    eye_h: float = 1.7,
    trail_hours: float = 48.0,
    show_track: bool = True,
    rise_set_text: str | None = None,
) -> dict:
    """渲染一帧 POV 地平线视角（清空 axes 后绘制 + 信息卡）。

    :param rise_set_text: 升起/落下时间摘要（如 "↑ 08-07 04:12\n↓ 08-20 16:40"）
    :return: {"az", "el", "apparent_radius_deg", "y_horizon", "fov_y", "earth_drawn"}
    """
    ax.clear()
    ax.set_xlim(0, POV_ASPECT)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_facecolor(_POV_COLORS["space"])
    fig.patch.set_facecolor(_POV_COLORS["space"])

    info = draw_pov_view(
        ax, t_mjd, lat_deg, lon_deg, h_m, R_eq, R_pol,
        view_az, tilt_deg, fov_x=fov_x, eye_h=eye_h,
        trail_hours=trail_hours, show_track=show_track,
    )
    az, el = info["az"], info["el"]
    app_r = info["apparent_radius_deg"]

    dt = datetime(1858, 11, 17, tzinfo=timezone.utc) + timedelta(days=float(t_mjd))
    time_str = dt.strftime("%Y-%m-%d %H:%M UTC")

    # ---- HUD 信息卡（四角，四周留 ≥10% 安全区，大字号紧凑卡片）----
    # 布局（transAxes，16:9 横屏）：
    #   左上 = 时间 UTC，右上 = 地球高度角，左下 = 地球方位角，右下 = 升起/落下（或视角）
    # 安全区：卡片距画布边 ≥10%（x≥0.10、x+w≤0.90、y≥0.11、y+h≤0.89）
    _card(ax, 0.2, 0.6, 0.25, 0.13, "时间 · UTC", time_str,
          fs_title=15, fs_value=18)
    _card(ax, 0.55, 0.6, 0.25, 0.13, "地球高度角",
          f"{el:+.2f}°", sub=f"视直径 {2*app_r:.2f}°",
          fs_title=15, fs_value=22)
    _card(ax, 0.2, 0.2, 0.25, 0.13, "地球方位角",
          f"{az:5.1f}°", sub=_fmt_latlon(lat_deg, lon_deg),
          fs_title=15, fs_value=21)
    if rise_set_text:
        _card(ax, 0.55, 0.2, 0.25, 0.13, "升起 / 落下 · UTC",
              rise_set_text, fs_title=15, fs_value=13)
    else:
        _card(ax, 0.55, 0.2, 0.25, 0.13, "视角",
              f"朝向 {view_az:5.1f}°", sub=f"俯仰 {tilt_deg:+.1f}°",
              fs_title=15, fs_value=16)
    # 顶部标题（居中，距顶 ≥10%）
    ax.text(0.5, 0.85, "月球表面 · 地球升起与落下", transform=ax.transAxes,
            fontsize=18, color=_POV_COLORS["text"], ha="center", va="top",
            weight="bold", zorder=11)
    ax.text(0.5, 0.80, f"第一人称 · 视场 {fov_x:.0f}° · 地平线即 0° 高度角",
            transform=ax.transAxes, fontsize=11, color=_POV_COLORS["text_dim"],
            ha="center", va="top", zorder=11)
    return info


def _fmt_latlon(lat_deg: float, lon_deg: float) -> str:
    """纬度/经度 → 简洁 N/S E/W 字符串。"""
    ns = "N" if lat_deg >= 0 else "S"
    ew = "E" if lon_deg >= 0 else "W"
    return f"{abs(lat_deg):.2f}°{ns} {abs(lon_deg):.2f}°{ew}"


def make_pov_figure(figsize=(16.0, 9.0), dpi=120):
    """创建 POV Figure（16:9 手机横屏，黑色画布）。"""
    fig = plt.figure(figsize=figsize, dpi=dpi)
    fig.patch.set_facecolor(_POV_COLORS["space"])
    ax = fig.add_axes([0, 0, 1, 1])
    return fig, ax


def render_pov_single(
    t_mjd: float,
    lat_deg: float,
    lon_deg: float,
    h_m: float,
    R_eq: float,
    R_pol: float,
    output_path: str,
    view_az: float,
    tilt_deg: float,
    fov_x: float = 100.0,
    eye_h: float = 1.7,
    dpi: int = 120,
    show_track: bool = True,
) -> str:
    """渲染单帧 POV 地平线视角 PNG（站在月面，看向地球方位角）。

    ``view_az`` 默认应传当前时刻地球方位角（也可手动指定朝向）。
    """
    setup_cjk_font()
    fig, ax = make_pov_figure((16.0, 9.0), dpi)
    render_pov_frame(fig, ax, t_mjd, lat_deg, lon_deg, h_m, R_eq, R_pol,
                     view_az, tilt_deg, fov_x=fov_x, eye_h=eye_h,
                     show_track=show_track)
    fig.savefig(output_path, dpi=dpi, facecolor=fig.get_facecolor())
    plt.close(fig)
    return output_path


def compute_rise_set(
    lat_deg: float,
    lon_deg: float,
    h_m: float,
    R_eq: float,
    R_pol: float,
    start_mjd: float,
    end_mjd: float,
    search_step_h: float = 0.25,
) -> list[dict]:
    """在 [start, end] 内搜索地球高度角过零（升起/落下）时刻。

    由于月球潮汐锁定，地球高度角仅在很小的角度范围内摆动
    （南极约 -7°~+6°，受天平动调制），过零即“升起/落下”。
    用粗扫 + 线性插值定位过零时刻，返回 [{t_mjd, type}] 列表
    （type="rise" 升、"set" 落）。

    :param search_step_h: 粗扫步长（小时），默认 15 分钟
    """
    t = make_time_grid(start_mjd, end_mjd, search_step_h)
    prev = None
    events: list[dict] = []
    for i, tm in enumerate(t):
        _, el = earth_az_el(lon_deg, lat_deg, h_m, tm, R_eq, R_pol)
        if prev is not None:
            if prev[1] < 0.0 <= el:
                # 升（负 → 正），线性插值
                f = -prev[1] / (el - prev[1])
                events.append({"t_mjd": prev[0] + f * (tm - prev[0]),
                               "type": "rise"})
            elif prev[1] > 0.0 >= el:
                f = prev[1] / (prev[1] - el)
                events.append({"t_mjd": prev[0] + f * (tm - prev[0]),
                               "type": "set"})
        prev = (tm, el)
    return events


def _fmt_rise_set(events: list[dict]) -> str:
    """升起/落下事件 → 简洁两行文字（如 '↑ 08-07 04:12\\n↓ 08-20 16:40'）。"""
    lines = []
    for ev in events[:2]:
        dt = datetime(1858, 11, 17, tzinfo=timezone.utc) + timedelta(
            days=float(ev["t_mjd"]))
        arrow = "↑" if ev["type"] == "rise" else "↓"
        lines.append(f"{arrow} {dt.strftime('%m-%d %H:%M')}")
    if not lines:
        return "近期无过零"
    return "\n".join(lines)


def render_pov_animation(
    t_mjd_grid: np.ndarray,
    lat_deg: float,
    lon_deg: float,
    h_m: float,
    R_eq: float,
    R_pol: float,
    output_path: str,
    view_az: float,
    tilt_deg: float,
    fov_x: float = 100.0,
    eye_h: float = 1.7,
    fps: int = 24,
    dpi: int = 240,
    trail_hours: float = 48.0,
    show_track: bool = True,
    preview_only: bool = False,
) -> str:
    """渲染 POV 地平线视角 4K 延时动画（16:9 手机横屏，MP4）。

    相机固定朝向地球方位角，地平线在画面中固定；地球随天平动在
    地平线附近缓慢升降（南极可见完整的升起/落下）。每帧时间步长由
    ``t_mjd_grid`` 决定（默认 3 小时，一帧 2-4 小时）。
    首帧前计算升起/落下时刻，显示在右下角信息卡。

    :param view_az: 相机朝向方位角（建议传动画时段内地球平均方位角）
    :param preview_only: 只渲染首帧 PNG，不编码视频
    """
    setup_cjk_font()
    fig, ax = make_pov_figure((16.0, 9.0), dpi)

    # 升起/落下摘要（动画时段内）
    rise_set = compute_rise_set(lat_deg, lon_deg, h_m, R_eq, R_pol,
                                float(t_mjd_grid[0]), float(t_mjd_grid[-1]))
    rs_text = _fmt_rise_set(rise_set)

    def update(i: int):
        render_pov_frame(
            fig, ax, float(t_mjd_grid[i]), lat_deg, lon_deg, h_m, R_eq, R_pol,
            view_az, tilt_deg, fov_x=fov_x, eye_h=eye_h,
            trail_hours=trail_hours, show_track=show_track,
            rise_set_text=rs_text,
        )
        return ()

    N = len(t_mjd_grid)
    if preview_only:
        update(0)
        png_path = os.path.splitext(output_path)[0] + "_preview.png"
        fig.savefig(png_path, dpi=dpi, facecolor=fig.get_facecolor())
        plt.close(fig)
        print(f"[earth_elevation] POV 动画首帧预览 → {png_path}")
        return png_path

    print(f"[earth_elevation] 渲染 POV 动画 {N} 帧 @ {fps} fps → {output_path}")
    anim = FuncAnimation(fig, update, frames=N, interval=1000.0 / fps, blit=False)
    writer = FFMpegWriter(
        fps=fps, codec="libx264", bitrate=-1,
        extra_args=[
            "-pix_fmt", "yuv420p",
            "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",
        ],
    )
    anim.save(output_path, writer=writer, dpi=dpi)
    plt.close(fig)
    print(f"[earth_elevation] DONE → {output_path}")
    return output_path


def compute_timeseries(
    lon_deg: float,
    lat_deg: float,
    h_m: float,
    start_mjd: float,
    end_mjd: float,
    interval_h: float,
    R_eq: float,
    R_pol: float,
) -> tuple[np.ndarray, np.ndarray]:
    """计算高度角时间序列。SPICE 非线程安全，单线程串行。"""
    t_mjd = make_time_grid(start_mjd, end_mjd, interval_h)
    elev = np.empty_like(t_mjd)
    for i, tm in enumerate(t_mjd):
        elev[i] = earth_elevation_angle(lon_deg, lat_deg, h_m, tm, R_eq, R_pol)
    return t_mjd, elev


def plot_timeseries(
    t_mjd: np.ndarray,
    elev: np.ndarray,
    current_range: tuple[float, float],
    output_dir: str,
    start_utc: str,
    end_utc: str,
) -> str:
    """时间序列图：横轴=日期，纵轴=地球高度角（度），标注 min/max 线。"""
    fig, ax = plt.subplots(figsize=(13, 5.5))
    # MJD → datetime
    dates = [datetime(1858, 11, 17, tzinfo=timezone.utc) + __import__(
        "datetime").timedelta(days=float(t)) for t in t_mjd]

    ax.plot(dates, elev, lw=0.9, color="#1f77b4", label="Earth elevation")
    lo, hi = current_range
    ax.axhline(hi, color="#d62728", ls="--", lw=1.2, label=f"max {hi:.2f}°")
    ax.axhline(lo, color="#2ca02c", ls="--", lw=1.2, label=f"min {lo:.2f}°")
    ax.set_xlabel("UTC date")
    ax.set_ylabel("Earth elevation (deg, ellipsoid reference)")
    ax.set_title(
        f"Earth elevation at Shackleton floor "
        f"({start_utc[:10]} → {end_utc[:10]})"
    )
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    fig.autofmt_xdate()
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    path = os.path.join(output_dir, "earth_elevation_timeseries.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_histogram(
    current_elev: np.ndarray,
    theoretical_elev: np.ndarray | None,
    current_range: tuple[float, float],
    theoretical_range: tuple[float, float] | None,
    output_dir: str,
) -> str:
    """直方图：x=高度角，y=频次，标注当前范围 vs 理论极值范围。"""
    fig, ax = plt.subplots(figsize=(11, 5.5))
    if theoretical_elev is not None:
        ax.hist(
            theoretical_elev, bins=180, color="#ff7f0e", alpha=0.55,
            label=f"Theoretical 18.6-yr (n={len(theoretical_elev)})",
        )
    ax.hist(
        current_elev, bins=90, color="#1f77b4", alpha=0.7,
        label=f"Current period (n={len(current_elev)})",
    )
    lo, hi = current_range
    ax.axvline(hi, color="#d62728", ls="--", lw=1.2,
               label=f"current max {hi:.2f}°")
    ax.axvline(lo, color="#2ca02c", ls="--", lw=1.2,
               label=f"current min {lo:.2f}°")
    if theoretical_range is not None:
        tlo, thi = theoretical_range
        ax.axvspan(tlo, thi, color="#ff7f0e", alpha=0.15,
                   label=f"theoretical range [{tlo:.2f}, {thi:.2f}]°")
    ax.set_xlabel("Earth elevation (deg)")
    ax.set_ylabel("Frequency")
    ax.set_title("Earth elevation distribution at Shackleton floor")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    path = os.path.join(output_dir, "earth_elevation_histogram.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def _series_to_records(t_mjd: np.ndarray, elev: np.ndarray) -> list[dict]:
    """时间序列 → [{t_utc, elev_deg}, ...]，供 JSON 输出。"""
    return [
        {"t_utc": mjd_to_utc_iso(float(t)), "elev_deg": round(float(e), 4)}
        for t, e in zip(t_mjd, elev)
    ]


def build_json(
    lon_deg: float,
    lat_deg: float,
    h_m: float,
    start_utc: str,
    end_utc: str,
    interval_h: float,
    t_cur: np.ndarray,
    elev_cur: np.ndarray,
    t_theo: np.ndarray,
    elev_theo: np.ndarray,
) -> dict:
    """组装 JSON 记录。"""
    kernels = [
        "latest_leapseconds.tls",
        "moon_pa_de421_1900-2050.bpc",
        "pck00011_n0066.tpc",
        "de442s.bsp",
    ]
    return {
        "site": {"lat_deg": lat_deg, "lon_deg": lon_deg, "h_m": h_m},
        "time_range": {"start_utc": start_utc, "end_utc": end_utc},
        "interval_h": interval_h,
        "spice_kernels": kernels,
        "current_period": {
            "min_elev_deg": round(float(elev_cur.min()), 4),
            "max_elev_deg": round(float(elev_cur.max()), 4),
            "mean_elev_deg": round(float(elev_cur.mean()), 4),
            "n_samples": int(len(elev_cur)),
            "timeseries": _series_to_records(t_cur, elev_cur),
        },
        "theoretical_extreme": {
            "period_years": 18.6,
            "start_utc": THEORETICAL_START,
            "end_utc": THEORETICAL_END,
            "min_elev_deg": round(float(elev_theo.min()), 4),
            "max_elev_deg": round(float(elev_theo.max()), 4),
            "n_samples": int(len(elev_theo)),
            "note": "完整天平动周期极值，供参考",
        },
        "evidence_level": "A",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--lat", type=float, default=SHACKLETON_LAT_DEG)
    p.add_argument("--lon", type=float, default=SHACKLETON_LON_DEG)
    p.add_argument("--h", type=float, default=0.0,
                   help="观测点高程（米，参考椭球面以上）")
    p.add_argument("--start", default="2026-08-01",
                   help="当前时段起始（UTC ISO，默认 2026-08-01）")
    p.add_argument("--end", default="2027-12-31",
                   help="当前时段结束（UTC ISO，默认 2027-12-31）")
    p.add_argument("--interval-h", type=float, default=1.0,
                   help="当前时段采样间隔（小时），默认 1.0")
    p.add_argument("--theoretical-extent-years", type=float, default=18.6,
                   help="理论极值时间范围（年），默认 18.6")
    p.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR,
                   help="输出目录，默认 figures-and-logs")

    # ---- 2D 天空视角（单帧 PNG / 4K 延时动画 MP4）----
    p.add_argument("--sky-view-time", default=None,
                   help="2D 天空视角单帧时刻（UTC ISO，如 2026-08-15T12:00:00）；"
                        "指定后输出 PNG 并跳过原时间序列分析")
    p.add_argument("--sky-animation", action="store_true",
                   help="输出 2D 天空视角 4K 延时动画（16:9 手机横屏）")
    p.add_argument("--sky-start", default="2026-08-01",
                   help="天空动画起始（UTC ISO，默认 2026-08-01）")
    p.add_argument("--sky-end", default="2026-08-31",
                   help="天空动画结束（UTC ISO，默认 2026-08-31）")
    p.add_argument("--sky-interval-h", type=float, default=3.0,
                   help="天空动画每帧时间步长（小时），默认 3.0（一帧 2-4 小时）")
    p.add_argument("--sky-fps", type=int, default=24,
                   help="天空动画帧率，默认 24")
    p.add_argument("--sky-dpi", type=int, default=240,
                   help="天空动画渲染 dpi，默认 240（16×9 → 3840×2160 4K）")
    p.add_argument("--sky-trail-h", type=float, default=48.0,
                   help="天空动画地球轨迹回看时长（小时），默认 48")
    p.add_argument("--sky-no-track", action="store_true",
                   help="天空图不绘制地球历史轨迹")
    p.add_argument("--sky-preview-only", action="store_true",
                   help="只渲染天空动画首帧 PNG（版面检查），不编码视频")

    # ---- POV 第一人称地平线视角（站在月面，看地球升起/落下）----
    p.add_argument("--pov-view-time", default=None,
                   help="POV 单帧时刻（UTC ISO）；指定后输出 PNG 并跳过原分析")
    p.add_argument("--pov-animation", action="store_true",
                   help="输出 POV 4K 延时动画（16:9 手机横屏，地球升起/落下）")
    p.add_argument("--pov-az", type=float, default=None,
                   help="相机朝向方位角（度，北=0）；默认自动指向地球方位")
    p.add_argument("--pov-tilt", type=float, default=4.0,
                   help="相机俯仰角（度，>0 向上看，默认 4.0）")
    p.add_argument("--pov-fov", type=float, default=100.0,
                   help="水平视场角（度，默认 100）")
    p.add_argument("--pov-eye-h", type=float, default=1.7,
                   help="眼高（米，默认 1.7）")
    p.add_argument("--pov-start", default="2026-08-01",
                   help="POV 动画起始（UTC ISO，默认 2026-08-01）")
    p.add_argument("--pov-end", default="2026-08-31",
                   help="POV 动画结束（UTC ISO，默认 2026-08-31）")
    p.add_argument("--pov-interval-h", type=float, default=3.0,
                   help="POV 动画每帧时间步长（小时），默认 3.0（一帧 2-4 小时）")
    p.add_argument("--pov-fps", type=int, default=24,
                   help="POV 动画帧率，默认 24")
    p.add_argument("--pov-dpi", type=int, default=240,
                   help="POV 动画渲染 dpi，默认 240（16×9 → 3840×2160 4K）")
    p.add_argument("--pov-trail-h", type=float, default=48.0,
                   help="POV 地球轨迹回看时长（小时），默认 48")
    p.add_argument("--pov-no-track", action="store_true",
                   help="POV 图不绘制地球历史轨迹")
    p.add_argument("--pov-preview-only", action="store_true",
                   help="只渲染 POV 动画首帧 PNG（版面检查），不编码视频")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    load_spice_kernels()
    try:
        R_eq, R_pol = moon_radii_m()
        print(f"[earth_elevation] Moon radii: R_eq={R_eq:.0f} m, R_pol={R_pol:.0f} m")

        # ---- 2D 天空视角单帧（任意纬度 / 任意时刻）----
        if args.sky_view_time is not None:
            t_mjd = utc_iso_to_mjd(args.sky_view_time)
            out = os.path.join(args.output_dir, "earth_sky_view.png")
            p = render_sky_view_single(
                t_mjd, args.lat, args.lon, args.h, R_eq, R_pol,
                out, dpi=args.sky_dpi, show_track=not args.sky_no_track,
            )
            print(f"[earth_elevation] 天空视角单帧 → {p}")
            print("[earth_elevation] DONE")
            return

        # ---- 2D 天空视角 4K 延时动画（手机横屏）----
        if args.sky_animation:
            sky_start = utc_iso_to_mjd(args.sky_start)
            sky_end = utc_iso_to_mjd(args.sky_end)
            t_grid = make_time_grid(sky_start, sky_end, args.sky_interval_h)
            out = os.path.join(args.output_dir, "earth_sky_view_animation.mp4")
            p = render_sky_view_animation(
                t_grid, args.lat, args.lon, args.h, R_eq, R_pol,
                out, fps=args.sky_fps, dpi=args.sky_dpi,
                trail_hours=args.sky_trail_h, show_track=not args.sky_no_track,
                preview_only=args.sky_preview_only,
            )
            print(f"[earth_elevation] 天空视角动画 → {p}")
            print("[earth_elevation] DONE")
            return

        # ---- POV 第一人称地平线视角单帧（站在月面，看地球）----
        if args.pov_view_time is not None:
            t_mjd = utc_iso_to_mjd(args.pov_view_time)
            az_now, el_now = earth_az_el(args.lon, args.lat, args.h,
                                         t_mjd, R_eq, R_pol)
            view_az = args.pov_az if args.pov_az is not None else az_now
            out = os.path.join(args.output_dir, "earth_pov_view.png")
            p = render_pov_single(
                t_mjd, args.lat, args.lon, args.h, R_eq, R_pol,
                out, view_az=view_az, tilt_deg=args.pov_tilt,
                fov_x=args.pov_fov, eye_h=args.pov_eye_h,
                dpi=args.pov_dpi, show_track=not args.pov_no_track,
            )
            print(f"[earth_elevation] POV 单帧 → {p}")
            print("[earth_elevation] DONE")
            return

        # ---- POV 4K 延时动画（地球升起/落下）----
        if args.pov_animation:
            pov_start = utc_iso_to_mjd(args.pov_start)
            pov_end = utc_iso_to_mjd(args.pov_end)
            t_grid = make_time_grid(pov_start, pov_end, args.pov_interval_h)
            # 相机默认指向动画时段内地球平均方位角
            if args.pov_az is not None:
                view_az = args.pov_az
            else:
                mid = t_grid[len(t_grid) // 2]
                az_mid, _ = earth_az_el(args.lon, args.lat, args.h,
                                        float(mid), R_eq, R_pol)
                view_az = az_mid
            out = os.path.join(args.output_dir, "earth_pov_animation.mp4")
            p = render_pov_animation(
                t_grid, args.lat, args.lon, args.h, R_eq, R_pol,
                out, view_az=view_az, tilt_deg=args.pov_tilt,
                fov_x=args.pov_fov, eye_h=args.pov_eye_h,
                fps=args.pov_fps, dpi=args.pov_dpi,
                trail_hours=args.pov_trail_h, show_track=not args.pov_no_track,
                preview_only=args.pov_preview_only,
            )
            print(f"[earth_elevation] POV 动画 → {p}")
            print("[earth_elevation] DONE")
            return

        # ---- 原有功能：高度角时间序列 + 直方图 + JSON ----
        start_mjd = utc_iso_to_mjd(args.start)
        end_mjd = utc_iso_to_mjd(args.end)
        print(
            f"[earth_elevation] current period {args.start} → {args.end} "
            f"(MJD {start_mjd:.2f} → {end_mjd:.2f}), interval {args.interval_h} h"
        )
        t_cur, elev_cur = compute_timeseries(
            args.lon, args.lat, args.h,
            start_mjd, end_mjd, args.interval_h, R_eq, R_pol,
        )
        print(
            f"[earth_elevation] current: min={elev_cur.min():.3f}°, "
            f"max={elev_cur.max():.3f}°, mean={elev_cur.mean():.3f}°, "
            f"n={len(elev_cur)}"
        )

        # 理论极值（天平动周期）
        theo_start = utc_iso_to_mjd(THEORETICAL_START)
        theo_end = utc_iso_to_mjd(THEORETICAL_END)
        theo_interval = 6.0  # 6 小时间隔
        t_theo, elev_theo = compute_timeseries(
            args.lon, args.lat, args.h,
            theo_start, theo_end, theo_interval, R_eq, R_pol,
        )
        print(
            f"[earth_elevation] theoretical {args.theoretical_extent_years}yr: "
            f"min={elev_theo.min():.3f}°, max={elev_theo.max():.3f}°, n={len(elev_theo)}"
        )

        data = build_json(
            args.lon, args.lat, args.h,
            args.start, args.end, args.interval_h,
            t_cur, elev_cur, t_theo, elev_theo,
        )
        json_path = os.path.join(args.output_dir, "earth_elevation_timeseries.json")
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        print(f"[earth_elevation] wrote {json_path}")

        png1 = plot_timeseries(
            t_cur, elev_cur,
            (float(elev_cur.min()), float(elev_cur.max())),
            args.output_dir, args.start, args.end,
        )
        png2 = plot_histogram(
            elev_cur, elev_theo,
            (float(elev_cur.min()), float(elev_cur.max())),
            (float(elev_theo.min()), float(elev_theo.max())),
            args.output_dir,
        )
        print(f"[earth_elevation] wrote {png1}")
        print(f"[earth_elevation] wrote {png2}")
        print("[earth_elevation] DONE")
    finally:
        cleanup_spice()


if __name__ == "__main__":
    main()
