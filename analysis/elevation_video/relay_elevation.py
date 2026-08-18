#!/usr/bin/env python3
"""P0-3b：鹊桥二号中继卫星在沙克尔顿三观测点的高度角 + 方位角时间序列。

扩展目标：为 ``horizon_mask.py``（P0-1）的方位角剖面图提供中继卫星
逐方位角高度角带，与地球高度角带（P0-2）并列对比。

科学口径：
- 鹊桥二号不是地月平动点晕轨道，而是**月球大椭圆太阳同步测控回归冻结轨道**
  （周期约 24.84 小时）。轨道传播复用 ``queqiao_relay.py`` 的月心二体模型。
- 轨道根数 (a, e, I, ω, T) 来自公开论文表 1（证据等级 A）；
  Ω, M 依赖任务时间（证据等级 B），此处沿用 ``queqiao_relay.py`` 的假设值。
- 所有结论携带证据等级，输出条件表述而非假精确数字。

算法：
1. 以月球赤道（IAU_MOON 极轴）为参考构造轨道状态矢量（M 处）
2. ``spiceypy.oscelt`` 转椭圆要素，``spiceypy.conics`` 二体传播到每个采样时刻
3. 每步 ``pxform(J2000→IAU_MOON)`` 转到月固系，对三个观测点（坑底中心 /
   偏北 1 km / 偏南 1 km，坐标与 ``horizon_mask.py`` 的 SITES 一致）用共享
   ENU 投影计算 (仰角, 方位角)
4. 输出每站点的时间序列 JSON + 一张三站点叠加图

产出：
- ``relay_elevation_{site}.json`` — 三站点时间序列（t_utc, elev_deg, az_deg,
  visible, visible_terrain, horizon_elev_deg）+ 每站点 min/max/mean 统计，
  以及两个可见比例：``visible_fraction``（纯几何，水平线之上即算）与
  ``visible_terrain_fraction``（考虑地形遮挡：高度角须大于该方位角地形
  地平线，由 ``--terrain`` 启用以 P0-1 ``horizon_mask_{suffix}.json`` 判定）
- ``relay_elevation_timeseries.png`` — 三站点仰角时间序列叠加图

用法::

    python analysis/elevation_video/relay_elevation.py \\
        --lat -89.67 --lon 129.78 --h 0.0 \\
        --start 2026-08-01 --end 2027-12-31 --interval-h 0.5 \\
        --output-dir figures-and-logs
"""

from __future__ import annotations

import argparse
import json
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

from analysis.elevation_video.queqiao_relay import (  # noqa: E402
    CONDITIONAL_STATEMENT,
    ORBIT_PARAMS,
    RAAN_ASSUMED_DEG,
    M_ASSUMED_DEG,
    RAAN_SOURCE,
    M_SOURCE,
    REFERENCES,
    aggregate_daily_stats,
    build_orbital_state,
    propagate_satellite,
)
from analysis.elevation_video.utils import (  # noqa: E402
    DEM_PATH,
    SHACKLETON_LAT_DEG,
    SHACKLETON_LON_DEG,
    cleanup_spice,
    enu_elevation_azimuth,
    load_spice_kernels,
    mjd_to_utc_iso,
    moon_radii_m,
    utc_iso_to_mjd,
    make_time_grid,
)

DEFAULT_OUTPUT_DIR = os.path.join(_REPO_ROOT, "logs", "elevation_video")

# 三个观测点：与 horizon_mask.py 的 SITES 完全一致（极射投影 +y 北 / -y 南平移 1 km）
# (site_name, y_offset_m, 文件名后缀)
SITES: tuple[tuple[str, float, str], ...] = (
    ("shackleton_floor_center", 0.0, "center"),
    ("shackleton_floor_north1km", 1000.0, "north1km"),
    ("shackleton_floor_south1km", -1000.0, "south1km"),
)

SITE_LABELS = {
    "shackleton_floor_center": "Crater floor center",
    "shackleton_floor_north1km": "1 km north",
    "shackleton_floor_south1km": "1 km south",
}


def site_lonlat(
    lon_center: float, lat_center: float, y_offset: float
) -> tuple[float, float]:
    """中心经纬度 + 极射投影 +y 平移 → 平移后站点经纬度（与 horizon_mask.py 一致）。"""
    from analysis.elevation_video.utils import (
        lonlat_to_polar_stereo_xy,
        polar_stereo_xy_to_lonlat,
    )

    x0, y0 = lonlat_to_polar_stereo_xy(lon_center, lat_center)
    return polar_stereo_xy_to_lonlat(x0, y0 + y_offset)


def load_terrain_horizon(
    output_dir: str,
    suffix: str,
    y_offset_m: float,
    dem_path: str | None = None,
    lon_center: float | None = None,
    lat_center: float | None = None,
    max_radius_m: float = 50000.0,
    n_az: int = 360,
    n_samples: int = 300,
    roi_size: int = 4096,
) -> tuple[np.ndarray | None, dict]:
    """加载（或回退计算）单观测点的地形地平线仰角表（度，360 方位）。

    :param output_dir: 产物目录，优先读取 ``horizon_mask_{suffix}.json``（P0-1
        已生成时直接复用，保证与地平线分析口径完全一致）
    :param suffix: 站点文件后缀（center / north1km / south1km）
    :param y_offset_m: 站点沿极射投影 +y 的平移量（米），计算回退时使用
    :param dem_path: DEM GeoTIFF 路径；为 None 且 JSON 不存在时返回 None
    :return: (horizon_elev_deg[n_az] 或 None, meta) — meta 含来源与站点坐标；
        数据不可用时返回 (None, {"available": False})
    """
    from analysis.elevation_video.utils import (
        lonlat_to_polar_stereo_xy,
        polar_stereo_xy_to_lonlat,
    )

    # 优先复用 P0-1 产物（口径一致：曲率修正 + 50 km 射线）
    json_path = os.path.join(output_dir, f"horizon_mask_{suffix}.json")
    if os.path.exists(json_path):
        try:
            with open(json_path, encoding="utf-8") as f:
                rec = json.load(f)
            azs = np.asarray(rec["azimuths_deg"], dtype=float)
            hor = np.asarray(rec["horizon_elev_deg"], dtype=float)
            if len(azs) == len(hor) and len(azs) > 1:
                meta = {
                    "available": True,
                    "source": f"horizon_mask_{suffix}.json",
                    "dem_source": rec.get("dem_source"),
                    "max_radius_m": rec.get("max_radius_m"),
                    "n_az": int(len(azs)),
                    "site_lat_deg": rec.get("lat_deg"),
                    "site_lon_deg": rec.get("lon_deg"),
                    "elev_m": rec.get("elev_m"),
                }
                return hor, meta
        except (json.JSONDecodeError, KeyError, TypeError, ValueError):
            pass  # 产物损坏则回退计算

    # 回退：直接从 DEM 计算（P0-1 尚未运行 / 产物缺失）
    if dem_path is None or lon_center is None or lat_center is None:
        return None, {"available": False, "source": "none"}

    from analysis.elevation_video.horizon_mask import compute_horizon_mask

    x0, y0 = lonlat_to_polar_stereo_xy(lon_center, lat_center)
    result = compute_horizon_mask(
        dem_path, x0, y0 + y_offset_m,
        max_radius_m=max_radius_m, n_az=n_az, n_samples=n_samples, roi_size=roi_size,
    )
    azs = np.asarray(result["azimuths_deg"], dtype=float)
    hor = np.asarray(result["horizon_elev_deg"], dtype=float)
    lon_site, lat_site = polar_stereo_xy_to_lonlat(x0, y0 + y_offset_m)
    meta = {
        "available": True,
        "source": f"computed from DEM {os.path.basename(dem_path)}",
        "dem_source": os.path.basename(dem_path),
        "max_radius_m": result["max_radius_m"],
        "n_az": int(len(azs)),
        "site_lat_deg": lat_site,
        "site_lon_deg": lon_site,
        "elev_m": result["elev_m"],
    }
    return hor, meta


def horizon_at_azimuth(
    horizon: np.ndarray | None, az_deg: np.ndarray, n_az: int = 360
) -> np.ndarray:
    """把逐方位角地平线仰角表插值到样本方位角（线性插值，周期闭合）。

    方位角约定与 horizon_mask 一致：北=0° 顺时针。horizon 为 None 时
    返回全 -inf（等效纯几何可见性，不含地形）。
    """
    if horizon is None:
        return np.full_like(np.asarray(az_deg, dtype=float), -np.inf)
    az = np.asarray(az_deg, dtype=float)
    h = np.asarray(horizon, dtype=float)
    if h.size == 0:
        return np.full_like(az, -np.inf)
    # 周期延拓：把表按方位角排好并线性插值
    if h.size == n_az:
        idx = np.clip(np.floor(az % 360.0).astype(int), 0, n_az - 1)
        idx_next = (idx + 1) % n_az
        frac = (az % 360.0) - np.floor(az % 360.0)
        return h[idx] * (1.0 - frac) + h[idx_next] * frac
    # 非 360 表：np.interp 周期化
    az_wrap = np.concatenate([np.asarray(range(n_az)) - n_az, np.asarray(range(n_az)),
                              np.asarray(range(n_az)) + n_az])
    h_wrap = np.concatenate([h, h, h])
    return np.interp(az % 360.0, az_wrap + n_az, h_wrap)


def compute_relay_geometry(
    lon_center: float,
    lat_center: float,
    h_m: float,
    start_mjd: float,
    end_mjd: float,
    interval_h: float,
    R_eq: float,
    R_pol: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """计算中继卫星在三个观测点的高度角 + 方位角时间序列。

    轨道只传播一次（与观测点无关），对三个站点分别做 ENU 投影。

    :return: (t_mjd, elev_deg (3, n), az_deg (3, n), site_latlon (3, 2))
    """
    import spiceypy

    t_mjd = make_time_grid(start_mjd, end_mjd, interval_h)

    # 轨道传播一次（月心 J2000 → IAU_MOON，米）
    et0 = spiceypy.unitim(float(t_mjd[0]) + 2400000.5, "JDTDB", "ET")
    state0, mu = build_orbital_state(et0)
    et_list = np.array(
        [spiceypy.unitim(float(t) + 2400000.5, "JDTDB", "ET") for t in t_mjd],
        dtype=float,
    )
    pos_j2000 = propagate_satellite(state0, mu, et0, et_list)
    pos_iau = np.empty_like(pos_j2000)
    for i, et in enumerate(et_list):
        rot = spiceypy.pxform("J2000", "IAU_MOON", et)
        pos_iau[i] = spiceypy.mxv(rot, pos_j2000[i]) * 1000.0  # km -> m

    # 三个站点经纬度
    site_latlon = np.array(
        [site_lonlat(lon_center, lat_center, off) for _, off, _ in SITES],
        dtype=float,
    )

    elev = np.empty((len(SITES), len(t_mjd)))
    az = np.empty((len(SITES), len(t_mjd)))
    for s, (lon, lat) in enumerate(site_latlon):
        for i, p in enumerate(pos_iau):
            e, a = enu_elevation_azimuth(p, float(lon), float(lat), h_m, R_eq, R_pol)
            elev[s, i] = e
            az[s, i] = a
    return t_mjd, elev, az, site_latlon


def build_json(
    site_index: int,
    lon_center: float,
    lat_center: float,
    h_m: float,
    start_utc: str,
    end_utc: str,
    interval_h: float,
    t_mjd: np.ndarray,
    elev: np.ndarray,
    az: np.ndarray,
    site_latlon: np.ndarray,
    horizon: np.ndarray | None = None,
    horizon_meta: dict | None = None,
) -> dict:
    """组装单个观测点的 JSON 记录（与 ``horizon_mask_{suffix}.json`` 一一对应）。

    :param site_index: SITES 下标，决定该记录对应的观测点
    :param horizon: 该站点逐方位角地形地平线仰角表（度）；None 表示无地形数据
    :param horizon_meta: 地形地平线来源元信息（source / dem_source / n_az 等）
    """
    site_name, _, suffix = SITES[site_index]
    lon, lat = float(site_latlon[site_index, 0]), float(site_latlon[site_index, 1])

    # 纯几何可见性（水平线之上即算，不含地形）
    elev_s = elev[site_index]
    az_s = az[site_index]
    visible_geom = elev_s > 0.0

    # 地形可见性：高度角须大于该方位角的地形地平线（有效遮挡 = max(horizon, 0)）
    horizon_avail = horizon is not None and horizon_meta is not None and horizon_meta.get("available")
    if horizon_avail:
        hor_at_az = horizon_at_azimuth(horizon, az_s)
        effective_obstacle = np.maximum(hor_at_az, 0.0)
        visible_terrain = elev_s > effective_obstacle
    else:
        hor_at_az = np.full_like(elev_s, np.nan)
        visible_terrain = visible_geom  # 无地形数据时退化为纯几何

    timeseries = [
        {
            "t_utc": mjd_to_utc_iso(float(t)),
            "elev_deg": round(float(e), 3),
            "az_deg": round(float(a), 3),
            "visible": bool(vg),
            "visible_terrain": bool(vt),
            "horizon_elev_deg": (
                round(float(h), 3) if np.isfinite(h) else None
            ),
        }
        for t, e, a, vg, vt, h in zip(
            t_mjd, elev_s, az_s, visible_geom, visible_terrain, hor_at_az
        )
    ]

    orbit_params = dict(ORBIT_PARAMS)
    orbit_params["raan_assumed_deg"] = RAAN_ASSUMED_DEG
    orbit_params["raan_source"] = RAAN_SOURCE
    orbit_params["mean_anomaly_assumed_deg"] = M_ASSUMED_DEG
    orbit_params["mean_anomaly_source"] = M_SOURCE

    terrain = {
        "available": bool(horizon_avail),
        "source": (horizon_meta or {}).get("source"),
        "dem_source": (horizon_meta or {}).get("dem_source"),
        "max_radius_m": (horizon_meta or {}).get("max_radius_m"),
        "n_az": (horizon_meta or {}).get("n_az"),
    }

    # 带地形的每日可见时长统计（与 queqiao_relay.aggregate_daily_stats 同口径；
    # 无地形数据时 visible_terrain 退化为纯几何 visible_geom）
    daily_stats = aggregate_daily_stats(t_mjd, visible_terrain, interval_h)

    return {
        "data_level": "A",
        "data_level_reason": (
            "轨道根数(a,e,I,ω,T)来自论文表1微分修正结果；Ω,M依赖任务时间，"
            "用假设值做量级估算（若 Ω/M 精确值不可得，可见性窗口相位为估算值）"
        ),
        "references": REFERENCES,
        "site": {
            "site_name": site_name,
            "lat_deg": round(lat, 6),
            "lon_deg": round(lon, 6),
            "h_m": h_m,
        },
        "center_site": {"lat_deg": lat_center, "lon_deg": lon_center, "h_m": h_m},
        "time_range": {"start_utc": start_utc, "end_utc": end_utc},
        "interval_h": interval_h,
        "orbit_parameters": orbit_params,
        "terrain": terrain,
        "site_stats": {
            "min_elev_deg": round(float(elev_s.min()), 3),
            "max_elev_deg": round(float(elev_s.max()), 3),
            "mean_elev_deg": round(float(elev_s.mean()), 3),
            "visible_fraction": round(float(visible_geom.mean()), 4),
            "visible_terrain_fraction": round(float(visible_terrain.mean()), 4),
            "longest_relay_arc_h": daily_stats["max_daily_visible_h"],
            "mean_daily_visible_h": daily_stats["mean_daily_visible_h"],
            "n_samples": int(len(t_mjd)),
        },
        "timeseries": timeseries,
        "conditional_statement": CONDITIONAL_STATEMENT,
        "evidence_level": "A",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }


def plot_timeseries(
    t_mjd: np.ndarray,
    elev: np.ndarray,
    output_dir: str,
    start_utc: str,
    end_utc: str,
) -> str:
    """三站点仰角时间序列叠加图（横轴=日期，纵轴=仰角，含 0° 地平线）。"""
    from datetime import timedelta

    dates = [
        datetime(1858, 11, 17, tzinfo=timezone.utc) + timedelta(days=float(t))
        for t in t_mjd
    ]
    colors = ["#1f77b4", "#d62728", "#2ca02c"]
    fig, ax = plt.subplots(figsize=(13, 5.5))
    for s, (site_name, _, _) in enumerate(SITES):
        label = SITE_LABELS.get(site_name, site_name)
        ax.plot(dates, elev[s], lw=0.6, color=colors[s], label=label)
    ax.axhline(0.0, color="k", ls="--", lw=1.0, label="0° horizon (visible above)")
    ax.set_xlabel("UTC date")
    ax.set_ylabel("Relay satellite elevation (deg)")
    ax.set_title(
        f"Queqiao-2 relay elevation at Shackleton sites "
        f"({start_utc[:10]} → {end_utc[:10]})"
    )
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    fig.autofmt_xdate()
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    path = os.path.join(output_dir, "relay_elevation_timeseries.png")
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
    # ── 地形遮挡（P0-1 horizon_mask 交叉分析）──
    p.add_argument(
        "--terrain", action="store_true", default=False,
        help="启用地形遮挡判定：高度角须大于该方位角的地形地平线（读取 "
             "horizon_mask_{suffix}.json，缺失时回退用 DEM 计算）",
    )
    p.add_argument(
        "--dem", default=None,
        help="DEM GeoTIFF 路径（回退计算时用；默认与 utils.DEM_PATH 相同）",
    )
    p.add_argument("--max-radius-m", type=float, default=50000.0,
                   help="回退计算时最大射线距离（米），默认 50000")
    p.add_argument("--n-az", type=int, default=360,
                   help="回退计算时方位角采样数，默认 360")
    p.add_argument("--n-samples", type=int, default=300,
                   help="回退计算时每方位角射线采样点数，默认 300")
    p.add_argument("--roi-size", type=int, default=4096,
                   help="回退计算时 ROI 降采样目标像素数，默认 4096")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    load_spice_kernels()
    try:
        R_eq, R_pol = moon_radii_m()
        start_mjd = utc_iso_to_mjd(args.start)
        end_mjd = utc_iso_to_mjd(args.end)
        print(
            f"[relay_elevation] period {args.start} → {args.end}, "
            f"interval {args.interval_h} h"
        )

        t_mjd, elev, az, site_latlon = compute_relay_geometry(
            args.lon, args.lat, args.h,
            start_mjd, end_mjd, args.interval_h, R_eq, R_pol,
        )

        # ── 地形地平线（逐站点，P0-1 口径）──
        dem_path = args.dem if args.dem is not None else DEM_PATH
        terrain_horizons: list[np.ndarray | None] = []
        terrain_metas: list[dict] = []
        if args.terrain:
            for site_name, y_offset, suffix in SITES:
                hor, meta = load_terrain_horizon(
                    args.output_dir, suffix, y_offset,
                    dem_path=dem_path if os.path.exists(dem_path) else None,
                    lon_center=args.lon, lat_center=args.lat,
                    max_radius_m=args.max_radius_m, n_az=args.n_az,
                    n_samples=args.n_samples, roi_size=args.roi_size,
                )
                terrain_horizons.append(hor)
                terrain_metas.append(meta)
                if hor is not None:
                    print(
                        f"[relay_elevation] terrain horizon for {site_name}: "
                        f"loaded from {meta['source']} "
                        f"(max {float(hor.max()):.2f}°)"
                    )
                else:
                    print(
                        f"[relay_elevation] terrain horizon for {site_name}: "
                        f"NOT available (no horizon JSON, no DEM)"
                    )
        else:
            terrain_horizons = [None] * len(SITES)
            terrain_metas = [{"available": False}] * len(SITES)

        for s, (site_name, _, suffix) in enumerate(SITES):
            vis_frac = (elev[s] > 0.0).mean() * 100.0
            if terrain_horizons[s] is not None:
                hor_at = horizon_at_azimuth(terrain_horizons[s], az[s])
                vis_terr = (elev[s] > np.maximum(hor_at, 0.0)).mean() * 100.0
            else:
                vis_terr = vis_frac
            print(
                f"[relay_elevation]   {site_name}: "
                f"elev [{elev[s].min():.2f}, {elev[s].max():.2f}]°, "
                f"mean {elev[s].mean():.2f}°, "
                f"visible {float(vis_frac):.1f}%, "
                f"terrain-visible {float(vis_terr):.1f}%"
            )

        for s, (site_name, _, suffix) in enumerate(SITES):
            data = build_json(
                s, args.lon, args.lat, args.h,
                args.start, args.end, args.interval_h,
                t_mjd, elev, az, site_latlon,
                horizon=terrain_horizons[s], horizon_meta=terrain_metas[s],
            )
            json_path = os.path.join(args.output_dir, f"relay_elevation_{suffix}.json")
            with open(json_path, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
            print(f"[relay_elevation] wrote {json_path}")

        png = plot_timeseries(
            t_mjd, elev, args.output_dir, args.start, args.end
        )
        print(f"[relay_elevation] wrote {png}")
        print("[relay_elevation] DONE")
    finally:
        cleanup_spice()


if __name__ == "__main__":
    main()
