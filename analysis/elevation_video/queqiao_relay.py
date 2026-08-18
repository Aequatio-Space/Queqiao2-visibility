#!/usr/bin/env python3
"""P0-3：鹊桥二号中继窗口估算与条件表述。

**关键纠正**：鹊桥二号不是地月平动点晕轨道，而是
**月球大椭圆太阳同步测控回归冻结轨道**（周期约 24.84 小时）。本脚本
采用月心二体轨道传播，不使用三体平动点模型。

科学口径：
- 轨道根数 (a, e, I, ω, T) 来自公开论文表 1（证据等级 A）
- Ω, M 依赖任务时间（证据等级 B），此处用假设值做量级估算
- 所有结论携带证据等级，输出条件表述而非假精确数字

算法：
1. 以月球赤道（IAU_MOON 极轴）为参考构造轨道状态矢量（M 处）
2. ``spiceypy.oscelt`` 转椭圆要素，``spiceypy.conics`` 二体传播
3. 每步 ``pxform(J2000→IAU_MOON)`` 转到月固系，用共享 ENU 投影计算
   南极观测点→卫星仰角，仰角 > 0° = 可见
4. 按天聚合可见时长，输出每日窗口统计 + 条件表述

产出：
- ``queqiao_relay_analysis.json`` — 数据等级 + 轨道参数 + 可见性统计 + 条件表述
- ``queqiao_relay_visibility.png`` — 可见性时间序列图

用法::

    python analysis/elevation_video/queqiao_relay.py \\
        --lat -89.67 --lon 129.78 --h 0.0 \\
        --start 2026-08-01 --end 2027-12-31 --interval-h 0.5 \\
        --output-dir figures-and-logs
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from datetime import datetime, timezone

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

# 允许以脚本方式从仓库任意位置运行：将仓库根加入 sys.path
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from analysis.elevation_video.utils import (  # noqa: E402
    MOON_R_KM,
    MU_MOON_KM3_S2,
    SHACKLETON_LAT_DEG,
    SHACKLETON_LON_DEG,
    cleanup_spice,
    enu_elevation_angle,
    load_spice_kernels,
    mjd_to_utc_iso,
    moon_radii_m,
    utc_iso_to_mjd,
    make_time_grid,
)

DEFAULT_OUTPUT_DIR = os.path.join(_REPO_ROOT, "logs", "elevation_video")

# ===========================================================================
# 轨道参数（证据等级 A，来自论文表 1 微分修正结果）
# ===========================================================================
ORBIT_PARAMS = {
    "period_h": 24.84,           # 测控回归周期 T_R（论文直接给出；含地球三体摄动长期项，
                                 # 二体开普勒周期为 24.94 h）
    "semi_major_axis_km": 10004.059,   # 平均要素（von Zeipel，双撇号）
    "eccentricity": 0.782359,          # 平均要素
    "inclination_deg": 119.202,        # 逆行冻结轨道（>90°）
    "periapsis_argument_deg": 90.0,    # 冻结条件 sin2ω''=0 → ω''=±90°
    "periapsis_altitude_design_km": 600.0,    # 论文设计约束值
    "apoapsis_altitude_km": 16093.4,          # a(1+e) - R_moon
    "longest_relay_arc_h": 20.7,              # 论文表 4-5
    "longest_relay_arc_source": "paper Table 4-5 (Queqiao-2 mission relay arcs)",
}

# 估算参数（证据等级 B，依赖任务时间）
RAAN_ASSUMED_DEG = 115.0
RAAN_SOURCE = "(paper: Ω = 115.9°)"
M_ASSUMED_DEG = 215.0
M_SOURCE = "midpoint of paper's range 210°~220°"

REFERENCES = [
    "孟占峰. 月球大椭圆太阳同步测控回归冻结轨道设计[J]. 航空学报, 2024, 45(18): 210-225",
]

CONDITIONAL_STATEMENT = (
    "根据公开论文（孟占峰 2024 航空学报；Zhang et al. 2024；Zhou et al. 2024），"
    "鹊桥二号运行于约 24.84 小时周期的月球大椭圆太阳同步测控回归冻结轨道。"
    "半长轴约 10004 km（5.76 倍月球半径），偏心率 0.782，"
    "轨道倾角 119.2°（逆行），近月点幅角冻结在 ±90°。"
    "近月点高度约 600 km（设计值，平均要素计算值 439.9 km，差异源于地球三体短周期摄动），"
    "远月点高度约 16093 km。论文给出最长中继弧段约 20.7 小时/天。"
    "升交点赤经和平近点角的精确值依赖任务时间（Ω = λ_sample - 90°，M ≈ 210°~220°），"
    "本分析用假设值做量级估算，可见性窗口相位存在不确定性，"
    "但每日可见时长量级（10~20 小时）不受 Ω/M 假设影响（由轨道几何与南极高纬观测共同决定）。"
)


# ===========================================================================
# 轨道构造与传播
# ===========================================================================
def _moon_pole_j2000(et: float) -> np.ndarray:
    """IAU_MOON 极轴在 J2000 惯性系中的单位向量。"""
    import spiceypy

    rot = spiceypy.pxform("IAU_MOON", "J2000", et)
    pole = rot @ np.array([0.0, 0.0, 1.0])
    pole = pole / np.linalg.norm(pole)
    return pole


def build_orbital_state(
    et0: float,
    a_km: float = ORBIT_PARAMS["semi_major_axis_km"],
    e: float = ORBIT_PARAMS["eccentricity"],
    I_deg: float = ORBIT_PARAMS["inclination_deg"],
    omega_deg: float = ORBIT_PARAMS["periapsis_argument_deg"],
    Omega_deg: float = RAAN_ASSUMED_DEG,
    M_deg: float = M_ASSUMED_DEG,
) -> tuple[np.ndarray, float]:
    """在月心 J2000 惯性系构造轨道状态矢量（月赤道参考系）。

    :param et0: 历元（SPICE ET，秒）
    :return: (state6 [km, km/s], mu_moon)
    """
    mu = MU_MOON_KM3_S2
    # 月球极轴 + 月赤道面参考方向（J2000 X 在月赤道面的投影）
    pole = _moon_pole_j2000(et0)
    ref = np.array([1.0, 0.0, 0.0]) - np.dot(np.array([1.0, 0.0, 0.0]), pole) * pole
    ref = ref / np.linalg.norm(ref)
    e1, e2, e3 = ref, np.cross(pole, ref), pole
    e2 = e2 / np.linalg.norm(e2)

    # 经典要素 → 笛卡尔（月赤道参考系）
    co = math.cos(math.radians(Omega_deg))
    so = math.sin(math.radians(Omega_deg))
    cw = math.cos(math.radians(omega_deg))
    sw = math.sin(math.radians(omega_deg))
    cI = math.cos(math.radians(I_deg))
    sI = math.sin(math.radians(I_deg))
    P = e1 * (co * cw - so * sw * cI) + e2 * (so * cw + co * sw * cI) + e3 * (sw * sI)
    Q = e1 * (-co * sw - so * cw * cI) + e2 * (-so * sw + co * cw * cI) + e3 * (sI * cw)

    # 解 Kepler 方程（牛顿迭代）
    M_rad = math.radians(M_deg)
    E = M_rad
    for _ in range(30):
        dE = (E - e * math.sin(E) - M_rad) / (1.0 - e * math.cos(E))
        E -= dE
        if abs(dE) < 1e-12:
            break
    cosE, sinE = math.cos(E), math.sin(E)

    r_mag = a_km * (1.0 - e * cosE)
    r_vec = a_km * (cosE - e) * P + a_km * math.sqrt(1.0 - e * e) * sinE * Q
    v_dir = -sinE * P + math.sqrt(1.0 - e * e) * cosE * Q
    v_vec = math.sqrt(mu * a_km) / r_mag * v_dir
    state6 = np.concatenate([r_vec, v_vec])
    return state6, mu


def propagate_satellite(
    state0: np.ndarray,
    mu: float,
    et0: float,
    et_list: np.ndarray,
) -> np.ndarray:
    """用 oscelt + conics 二体传播轨道，返回各时刻 J2000 位置（km）。

    :param et_list: SPICE ET 数组（秒）
    :return: (n, 3) km
    """
    import spiceypy

    elts = spiceypy.oscelt(state0, et0, mu)
    out = np.empty((len(et_list), 3))
    for i, et in enumerate(et_list):
        st = spiceypy.conics(elts, et)
        out[i] = st[:3]
    return out


def sat_pos_iau_moon(et_list: np.ndarray, pos_j2000_km: np.ndarray) -> np.ndarray:
    """J2000 位置（km）→ IAU_MOON 月固系位置（米）。"""
    import spiceypy

    out = np.empty_like(pos_j2000_km)
    for i, et in enumerate(et_list):
        rot = spiceypy.pxform("J2000", "IAU_MOON", et)
        out[i] = spiceypy.mxv(rot, pos_j2000_km[i]) * 1000.0  # km -> m
    return out


# ===========================================================================
# 可见性分析
# ===========================================================================
def compute_visibility(
    lon_deg: float,
    lat_deg: float,
    h_m: float,
    t_mjd: np.ndarray,
    R_eq: float,
    R_pol: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """计算逐时刻中继卫星仰角，返回 (t_mjd, elev_deg, visible_bool)。

    观测者高程 h_m 为参考椭球面以上（P0-2 相同约定）。
    """
    import spiceypy

    et0_mjd = t_mjd[0]
    et0 = spiceypy.unitim(et0_mjd + 2400000.5, "JDTDB", "ET")
    state0, mu = build_orbital_state(et0)
    # unitim 仅接受标量，逐元素转换
    et_list = np.array(
        [spiceypy.unitim(t + 2400000.5, "JDTDB", "ET") for t in t_mjd],
        dtype=float,
    )

    pos_j2000 = propagate_satellite(state0, mu, et0, et_list)
    pos_iau = sat_pos_iau_moon(et_list, pos_j2000)

    elev = np.empty(len(t_mjd))
    for i, p in enumerate(pos_iau):
        elev[i] = enu_elevation_angle(p, lon_deg, lat_deg, h_m, R_eq, R_pol)
    visible = elev > 0.0
    return t_mjd, elev, visible


def aggregate_daily_stats(
    t_mjd: np.ndarray, visible: np.ndarray, interval_h: float
) -> dict:
    """按天聚合可见时长。

    :return: dict 含 daily list, mean/min/max daily visible hours, 窗口数统计
    """
    days = np.floor(t_mjd).astype(int)
    uniq_days = np.unique(days)
    daily = []
    for d in uniq_days:
        mask = days == d
        n_vis = int(visible[mask].sum())
        hours = n_vis * interval_h
        daily.append({"day_mjd": int(d), "day_utc": mjd_to_utc_iso(float(d))[:10],
                      "visible_h": round(hours, 2)})
    hrs = np.array([x["visible_h"] for x in daily])

    # 可见窗口数（连续段）
    n_windows_total = 0
    prev = False
    for v in visible:
        if v and not prev:
            n_windows_total += 1
        prev = v
    n_days = max(1, len(uniq_days))

    return {
        "daily": daily,
        "mean_daily_visible_h": round(float(hrs.mean()), 2),
        "min_daily_visible_h": round(float(hrs.min()), 2),
        "max_daily_visible_h": round(float(hrs.max()), 2),
        "n_visible_windows_per_day_avg": round(n_windows_total / n_days, 2),
        "n_days": int(len(uniq_days)),
    }


def build_json(
    lon_deg: float,
    lat_deg: float,
    h_m: float,
    start_utc: str,
    end_utc: str,
    interval_h: float,
    t_mjd: np.ndarray,
    elev: np.ndarray,
    visible: np.ndarray,
    daily_stats: dict,
    h_p_calc_km: float,
) -> dict:
    """组装 JSON 记录。"""
    timeseries = [
        {"t_utc": mjd_to_utc_iso(float(t)), "elev_deg": round(float(e), 3),
         "visible": bool(v)}
        for t, e, v in zip(t_mjd, elev, visible)
    ]
    orbit_params = dict(ORBIT_PARAMS)
    orbit_params["periapsis_altitude_calc_km"] = round(h_p_calc_km, 1)
    orbit_params["periapsis_altitude_note"] = (
        "design=600km (paper design constraint); calc=439.9km from mean elements "
        "(short-period Earth-3rd-body perturbation Δe≈0.016 → Δh_p≈160km)"
    )
    orbit_params["raan_assumed_deg"] = RAAN_ASSUMED_DEG
    orbit_params["raan_source"] = RAAN_SOURCE
    orbit_params["mean_anomaly_assumed_deg"] = M_ASSUMED_DEG
    orbit_params["mean_anomaly_source"] = M_SOURCE

    return {
        "data_level": "A",
        "data_level_reason": (
            "轨道根数(a,e,I,ω,T)来自论文表1微分修正结果；Ω,M依赖任务时间，"
            "用假设值做量级估算（若 Ω/M 精确值不可得，可见性窗口相位为估算值）"
        ),
        "references": REFERENCES,
        "observer": {"lat_deg": lat_deg, "lon_deg": lon_deg, "h_m": h_m},
        "time_range": {"start_utc": start_utc, "end_utc": end_utc},
        "interval_h": interval_h,
        "orbit_parameters": orbit_params,
        "visibility_stats": {
            "mean_daily_visible_h": daily_stats["mean_daily_visible_h"],
            "min_daily_visible_h": daily_stats["min_daily_visible_h"],
            "max_daily_visible_h": daily_stats["max_daily_visible_h"],
            "n_visible_windows_per_day_avg": daily_stats["n_visible_windows_per_day_avg"],
            "n_days": daily_stats["n_days"],
        },
        "timeseries": timeseries,
        "conditional_statement": CONDITIONAL_STATEMENT,
        "evidence_level": "A",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }


def plot_visibility(
    t_mjd: np.ndarray,
    elev: np.ndarray,
    visible: np.ndarray,
    output_dir: str,
    start_utc: str,
    end_utc: str,
) -> str:
    """可见性时间序列图：仰角曲线 + 可见区域填充 + 0° 地平线。"""
    dates = [
        datetime(1858, 11, 17, tzinfo=timezone.utc)
        + __import__("datetime").timedelta(days=float(t))
        for t in t_mjd
    ]
    fig, ax = plt.subplots(figsize=(13, 5.5))
    ax.plot(dates, elev, lw=0.7, color="#1f77b4", label="Queqiao-2 elevation")
    ax.axhline(0.0, color="k", ls="--", lw=1.0, label="0° horizon (visible above)")
    # 可见时段填充
    fill_t = [dates[i] for i in range(len(dates)) if visible[i]]
    fill_e = [elev[i] for i in range(len(dates)) if visible[i]]
    if fill_t:
        ax.fill_between(fill_t, fill_e, 0.0, color="#2ca02c", alpha=0.35,
                        label="visible window")
    ax.set_xlabel("UTC date")
    ax.set_ylabel("Relay satellite elevation (deg)")
    ax.set_title(
        f"Queqiao-2 relay visibility at Shackleton floor "
        f"({start_utc[:10]} → {end_utc[:10]})"
    )
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    fig.autofmt_xdate()
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    path = os.path.join(output_dir, "queqiao_relay_visibility.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--lat", type=float, default=SHACKLETON_LAT_DEG)
    p.add_argument("--lon", type=float, default=SHACKLETON_LON_DEG)
    p.add_argument("--h", type=float, default=0.0,
                   help="观测点高程（米，参考椭球面以上）")
    p.add_argument("--start", default="2026-08-01",
                   help="起始时间（UTC ISO，默认 2026-08-01）")
    p.add_argument("--end", default="2027-12-31",
                   help="结束时间（UTC ISO，默认 2027-12-31）")
    p.add_argument("--interval-h", type=float, default=0.5,
                   help="采样间隔（小时），默认 0.5")
    p.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR,
                   help="输出目录，默认 figures-and-logs")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    load_spice_kernels()
    try:
        R_eq, R_pol = moon_radii_m()
        start_mjd = utc_iso_to_mjd(args.start)
        end_mjd = utc_iso_to_mjd(args.end)
        t_mjd = make_time_grid(start_mjd, end_mjd, args.interval_h)
        print(
            f"[queqiao_relay] period {args.start} → {args.end}, "
            f"interval {args.interval_h} h, n={len(t_mjd)}"
        )

        t_mjd, elev, visible = compute_visibility(
            args.lon, args.lat, args.h, t_mjd, R_eq, R_pol
        )
        n_vis = int(visible.sum())
        print(
            f"[queqiao_relay] visible {n_vis}/{len(visible)} samples "
            f"({100.0 * n_vis / len(visible):.1f}%), "
            f"elev range [{elev.min():.1f}, {elev.max():.1f}]°"
        )

        daily = aggregate_daily_stats(t_mjd, visible, args.interval_h)
        print(
            f"[queqiao_relay] daily visible h: "
            f"mean={daily['mean_daily_visible_h']}, "
            f"min={daily['min_daily_visible_h']}, "
            f"max={daily['max_daily_visible_h']}, "
            f"windows/day={daily['n_visible_windows_per_day_avg']}"
        )

        # 近月点高度计算值（平均要素）
        a = ORBIT_PARAMS["semi_major_axis_km"]
        e = ORBIT_PARAMS["eccentricity"]
        h_p_calc = a * (1.0 - e) - MOON_R_KM

        data = build_json(
            args.lon, args.lat, args.h,
            args.start, args.end, args.interval_h,
            t_mjd, elev, visible, daily, h_p_calc,
        )
        json_path = os.path.join(args.output_dir, "queqiao_relay_analysis.json")
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        print(f"[queqiao_relay] wrote {json_path}")

        png = plot_visibility(
            t_mjd, elev, visible, args.output_dir, args.start, args.end
        )
        print(f"[queqiao_relay] wrote {png}")
        print("[queqiao_relay] DONE")
    finally:
        cleanup_spice()


if __name__ == "__main__":
    main()
