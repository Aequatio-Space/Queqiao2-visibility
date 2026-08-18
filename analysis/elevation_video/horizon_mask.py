#!/usr/bin/env python3
"""P0-1：沙克尔顿坑地形遮挡角（horizon mask）计算。

核心算法（与仿真器 ``illumination._horizon_mask_single_dem`` 独立，满足 P0 需求）：
- 预读 ROI 窗口到 numpy 数组（``rasterio.windows.from_bounds`` + ``out_shape`` 降采样）
- ``scipy.interpolate.RegularGridInterpolator`` 批量插值高程
- 射线外推至 ``max_radius_m=50000``（远山脊更高更遮挡）
- 30-50 km 范围内做月面曲率修正（``utils.spherical_horizon_angle_vec``）
- 方位角 0-359°（北=0°，顺时针），每方位对数间隔 300 个采样点

产出三个观测点（坑底中心 / 偏北 1 km / 偏南 1 km）的：
- horizon_mask_{site}.json — 地平线数据 + 统计
- horizon_mask_profile.png — 2D 剖面图（三站点地形地平线 + 地球高度角带，
  不含中继卫星叠加；中继分析见 relay_elevation.py / queqiao_relay.py）
- horizon_mask_polar.png — 极坐标地平线图（含地球高度角带，不含中继卫星）

可用 ``--sites`` 选择观测点子集；仅选一个站点时（如 ``--sites center``
只画坑底中心），剖面图/极坐标图输出为带站点角标的单站文件：
- horizon_mask_{site}_profile.png
- horizon_mask_{site}_polar.png

图表叠加的地球高度角范围来自兄弟分析产物 ``earth_elevation_timeseries.json``
（P0-2，按依赖顺序先运行即可自动加载；不存在时仅绘制地形地平线）。

图表均使用中文字体（Alibaba PuHuiTi 3.0，见 ``setup_chinese_font``）。

用法::

    python analysis/elevation_video/horizon_mask.py \\
        --dem /path/to/ldem_87s_5mpp.tif \\
        --lat -89.67 --lon 129.78 \\
        --max-radius-m 50000 --n-az 360 --n-samples 300 \\
        --output-dir figures-and-logs

只画坑底中心的极坐标图与剖面（输出文件带 center 角标）::

    python analysis/elevation_video/horizon_mask.py \\
        --dem /path/to/ldem_87s_5mpp.tif \\
        --lat -89.67 --lon 129.78 --sites center \\
        --output-dir figures-and-logs
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import warnings
from datetime import datetime, timezone

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.windows import from_bounds
from scipy.interpolate import RegularGridInterpolator

# 允许以脚本方式从仓库任意位置运行：将仓库根加入 sys.path
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from analysis.elevation_video.utils import (  # noqa: E402
    DEM_PATH,
    SHACKLETON_LAT_DEG,
    SHACKLETON_LON_DEG,
    lonlat_to_polar_stereo_xy,
    spherical_horizon_angle_vec,
)

# rasterio 窗口读取触发 NumPy 2.5 无害 shape-set 警告
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# 三个观测点：坑底中心 + 沿极射投影 +y（北）/ -y（南）平移 1 km
# (site_name, y_offset_m, 文件名后缀)
SITES: tuple[tuple[str, float, str], ...] = (
    ("shackleton_floor_center", 0.0, "center"),
    ("shackleton_floor_north1km", 1000.0, "north1km"),
    ("shackleton_floor_south1km", -1000.0, "south1km"),
)

# 站点中文名（图表图例 / 单站图标题）
SITE_LABELS: dict[str, str] = {
    "shackleton_floor_center": "坑底中心",
    "shackleton_floor_north1km": "坑底偏北 1 km",
    "shackleton_floor_south1km": "坑底偏南 1 km",
}

DEFAULT_OUTPUT_DIR = os.path.join(_REPO_ROOT, "logs", "elevation_video")

# ---------------------------------------------------------------------------
# 中文字体：Alibaba PuHuiTi 3.0（阿里巴巴普惠体）
# ---------------------------------------------------------------------------
# 优先注册用户字体目录中的 PuHuiTi 字重（55 Regular / 85 Bold），
# 找不到时回退到系统已安装的同一族名，再回退到 macOS 常见 CJK 字体。
_PUHUITI_FONT_PATHS = (
    os.path.expanduser("~/Library/Fonts/AlibabaPuHuiTi-3-55-Regular.ttf"),
    os.path.expanduser("~/Library/Fonts/AlibabaPuHuiTi-3-85-Bold.ttf"),
)


def setup_chinese_font() -> str | None:
    """配置 matplotlib 使用 Alibaba PuHuiTi 3.0，返回选中的字体名（找不到返回 None）。

    显式 ``fontManager.addfont`` 注册，使未安装该字体的环境也能直接加载；
    同时设置全局 rcParams，图例、坐标轴与标题统一使用中文字体。
    """
    for path in _PUHUITI_FONT_PATHS:
        if os.path.exists(path):
            try:
                fm.fontManager.addfont(path)
            except (OSError, RuntimeError, ValueError):
                pass

    avail = {f.name for f in fm.fontManager.ttflist}
    chosen = next(
        (c for c in ("Alibaba PuHuiTi 3.0",) if c in avail), None
    )
    if chosen is None:
        for fallback in ("Hiragino Sans GB", "Songti SC", "PingFang SC",
                         "STHeiti", "SimHei", "Arial Unicode MS",
                         "Noto Sans CJK SC"):
            if fallback in avail:
                chosen = fallback
                break
    if chosen is None:
        print("警告：未找到中文字体，中文可能显示为方块", file=sys.stderr)
        return None
    plt.rcParams["font.sans-serif"] = [chosen, "DejaVu Sans"]
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["axes.unicode_minus"] = False
    return chosen


_CHINESE_FONT = setup_chinese_font()


def load_roi_interpolator(
    ds: rasterio.io.DatasetReader,
    x0: float,
    y0: float,
    max_radius_m: float,
    roi_size: int,
) -> tuple[RegularGridInterpolator, np.ndarray, np.ndarray]:
    """预读以 (x0, y0) 为中心、半径 max_radius_m 的 ROI 窗口并构造插值器。

    :param ds: 已打开的 rasterio 数据集（DEM，极射投影）
    :param x0, y0: 观测点在极射投影坐标系中的位置（米）
    :param max_radius_m: ROI 半径（米），也是最大射线距离
    :param roi_size: 降采样后窗口长边的目标像素数
    :return: (interpolator, ys, xs) — ys/xs 为升序网格坐标（米）
    """
    bounds = ds.bounds
    left = max(x0 - max_radius_m, float(bounds.left))
    right = min(x0 + max_radius_m, float(bounds.right))
    bottom = max(y0 - max_radius_m, float(bounds.bottom))
    top = min(y0 + max_radius_m, float(bounds.top))
    if left >= right or bottom >= top:
        raise ValueError(
            f"ROI window empty at ({x0:.1f}, {y0:.1f}): observer outside DEM extent "
            f"{[bounds.left, bounds.bottom, bounds.right, bounds.top]}"
        )

    win = from_bounds(left, bottom, right, top, ds.transform).round_offsets()
    w = int(win.width)
    h = int(win.height)
    if w < 2 or h < 2:
        raise ValueError(f"ROI window too small ({w}x{h}) at ({x0:.1f}, {y0:.1f})")

    scale = min(roi_size / w, roi_size / h, 1.0)
    out_w = max(2, int(round(w * scale)))
    out_h = max(2, int(round(h * scale)))
    data = ds.read(
        1, window=win, out_shape=(out_h, out_w), resampling=Resampling.bilinear,
        masked=False,
    )
    data = np.asarray(data, dtype=np.float64)

    # 窗口子变换 + 降采样因子 → 每个降采样像素中心的 (x, y) 坐标
    sub_t = ds.window_transform(win)
    ratio_x = w / out_w
    ratio_y = h / out_h
    step_x = sub_t.a * ratio_x
    step_y = abs(sub_t.e) * ratio_y
    xs = sub_t.c + (np.arange(out_w) + 0.5) * step_x
    ys = (sub_t.f + 0.5 * sub_t.e) - np.arange(out_h) * step_y

    # RegularGridInterpolator 要求 points 升序；DEM 行自上而下 y 递减 → 翻转
    ys = ys[::-1]
    data = data[::-1, :]

    interp = RegularGridInterpolator(
        (ys, xs), data, method="linear", bounds_error=False, fill_value=np.nan
    )
    return interp, ys, xs


def compute_horizon_mask(
    dem_path: str,
    x0: float,
    y0: float,
    max_radius_m: float = 50000.0,
    n_az: int = 360,
    n_samples: int = 300,
    roi_size: int = 4096,
) -> dict:
    """计算单个观测点的地平线遮挡角（含曲率修正）。

    :param dem_path: DEM GeoTIFF 路径（南极极射投影）
    :param x0, y0: 观测点在极射投影坐标系中的位置（米）
    :return: dict 含 azimuths_deg, horizon_elev_deg, 统计量, 截断信息, 观测点高程
    """
    with rasterio.open(dem_path) as ds:
        interp, _, _ = load_roi_interpolator(ds, x0, y0, max_radius_m, roi_size)

        # 观测者高程（DEM 相对参考椭球）
        h0 = float(interp((y0, x0)))
        if not np.isfinite(h0):
            raise RuntimeError(
                f"Observer elevation is NaN at polar-stereo ({x0:.1f}, {y0:.1f})"
            )

        # ── 批量射线采样（numpy 向量化，一次插值全部 360×300 点）──
        az_deg = np.arange(n_az, dtype=float)
        az_rad = np.radians(az_deg)
        # 北=0° 顺时针：dx=sin(az)（东），dy=cos(az)（北）
        distances = np.geomspace(5.0, max_radius_m, n_samples)
        x_pts = x0 + np.outer(np.sin(az_rad), distances)  # (n_az, n_samples)
        y_pts = y0 + np.outer(np.cos(az_rad), distances)
        pts = np.stack([y_pts, x_pts], axis=-1).reshape(-1, 2)
        z = interp(pts).reshape(n_az, n_samples)

        # 曲率修正仰角（度）；无效采样点 → -inf（不影响 max）
        angles = spherical_horizon_angle_vec(distances[None, :], h0, z)
        valid = np.isfinite(angles)
        angles = np.where(valid, angles, -np.inf)
        horizon = angles.max(axis=1)

        # 截断检测：每个方位角最后一个有效采样距离
        d_last = np.where(valid, distances[None, :], 0.0).max(axis=1)
        truncated = d_last < 0.98 * max_radius_m
        truncated_azimuths = az_deg[truncated].tolist()

        # 观测点 XY（供图表/JSON 使用）
        dem_source = os.path.basename(dem_path)

    max_elev = float(horizon.max())
    return {
        "x_stereo_m": float(x0),
        "y_stereo_m": float(y0),
        "elev_m": float(h0),
        "azimuths_deg": az_deg.tolist(),
        "horizon_elev_deg": np.round(horizon, 4).tolist(),
        "max_elev_deg": round(max_elev, 4),
        "mean_elev_deg": round(float(horizon.mean()), 4),
        "median_elev_deg": round(float(np.median(horizon)), 4),
        "azimuth_of_max_deg": round(float(az_deg[int(np.argmax(horizon))]), 1),
        "truncated_azimuths_deg": truncated_azimuths,
        "n_truncated_az": len(truncated_azimuths),
        "dem_source": dem_source,
        "max_radius_m": float(max_radius_m),
        "n_az": int(n_az),
        "n_samples": int(n_samples),
        "roi_size": int(roi_size),
    }


def load_earth_elev_range(output_dir: str) -> list[float] | None:
    """从 P0-2 产物读取地球高度角当前范围 [min, max]（度）。

    若 ``earth_elevation_timeseries.json`` 尚不存在则返回 None
    （P0-1 可独立运行，交叉对比在 P0-2 完成后由 run_all 触发）。
    """
    path = os.path.join(output_dir, "earth_elevation_timeseries.json")
    if not os.path.exists(path):
        return None
    try:
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
        cp = data.get("current_period", {})
        lo = cp.get("min_elev_deg")
        hi = cp.get("max_elev_deg")
        if lo is None or hi is None:
            return None
        return [float(lo), float(hi)]
    except (json.JSONDecodeError, KeyError, TypeError):
        return None


def make_json_record(
    site_name: str,
    lon_deg: float,
    lat_deg: float,
    result: dict,
    earth_range: list[float] | None,
) -> dict:
    """组装 JSON 记录。"""
    max_elev = result["max_elev_deg"]
    record: dict = {
        "site_name": site_name,
        "lat_deg": lat_deg,
        "lon_deg": lon_deg,
        "x_stereo_m": result["x_stereo_m"],
        "y_stereo_m": result["y_stereo_m"],
        "elev_m": result["elev_m"],
        "dem_source": result["dem_source"],
        "max_radius_m": result["max_radius_m"],
        "n_az": result["n_az"],
        "azimuths_deg": result["azimuths_deg"],
        "horizon_elev_deg": result["horizon_elev_deg"],
        "max_elev_deg": max_elev,
        "mean_elev_deg": result["mean_elev_deg"],
        "median_elev_deg": result["median_elev_deg"],
        "azimuth_of_max_deg": result["azimuth_of_max_deg"],
        "truncated_azimuths_deg": result["truncated_azimuths_deg"],
        "n_truncated_az": result["n_truncated_az"],
        "earth_elev_range_deg": earth_range,
        "earth_always_below_horizon": (
            earth_range is not None and max_elev > earth_range[1]
        ),
        "evidence_level": "A",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }
    return record


def plot_profile(
    records: list[dict],
    output_dir: str,
    earth_range: list[float] | None,
    suffix: str | None = None,
) -> str:
    """2D 剖面图：横轴方位角 0-360°，纵轴仰角。

    多站点模式绘制全部所选站点（默认三站点叠加）+ 地球高度角带；
    单站点模式（如 ``--sites center``）仅绘制该站点，输出文件带站点角标
    （``horizon_mask_{suffix}_profile.png``，多站点仍为 ``horizon_mask_profile.png``）。
    （不含中继卫星叠加；中继分析见 relay_elevation.py / queqiao_relay.py。）
    """
    fig, ax = plt.subplots(figsize=(12, 5.5))
    colors = ["#1f77b4", "#d62728", "#2ca02c"]

    for rec, color in zip(records, colors):
        az = np.array(rec["azimuths_deg"])
        hor = np.array(rec["horizon_elev_deg"])
        label = SITE_LABELS.get(rec["site_name"], rec["site_name"])
        ax.plot(az, hor, color=color, lw=1.4, label=label)
        k = int(np.argmax(hor))
        ax.scatter([az[k]], [hor[k]], color=color, s=28, zorder=5)
        ax.annotate(
            f"最大 {hor[k]:.1f}° @ {az[k]:.0f}°",
            xy=(az[k], hor[k]),
            xytext=(az[k] + 6, hor[k] + 0.9),
            fontsize=15,
            color=color,
            arrowprops=dict(arrowstyle="->", lw=0.7, color=color),
        )
        k = int(np.argmin(hor))
        ax.scatter([az[k]], [hor[k]], color=color, s=28, zorder=5)
        ax.annotate(
            f"最小 {hor[k]:.1f}° @ {az[k]:.0f}°",
            xy=(az[k], hor[k]),
            xytext=(az[k] + 6, hor[k] + 0.9),
            fontsize=15,
            color=color,
            arrowprops=dict(arrowstyle="->", lw=0.7, color=color),
        )

    if earth_range is not None:
        lo, hi = earth_range
        ax.axhspan(lo, hi, color="tab:orange", alpha=0.22, zorder=0,
                   label=f"地球高度角范围 [{lo:.1f}°, {hi:.1f}°]")
        ax.axhline(hi, color="tab:orange", ls="--", lw=1.0, alpha=0.8)
        ax.axhline(lo, color="tab:orange", ls="--", lw=1.0, alpha=0.8)

    # ── y 轴范围自适应：完整显示遮挡角曲线，同时保留 0 线与地球高度角带 ──
    vals = np.concatenate(
        [np.asarray(rec["horizon_elev_deg"], dtype=float) for rec in records]
    )
    vals = vals[np.isfinite(vals)]
    y_bottom = -0.5
    if earth_range is not None:
        y_bottom = min(y_bottom, earth_range[0] - 1.0)
    y_top = float(np.max(vals)) + 2.5  # 为 max 标注留出空间
    if not np.isfinite(y_top):
        y_top = 15.0

    ax.set_xlabel("方位角（度，自北顺时针）", fontsize=18)
    ax.set_ylabel("地平线仰角（度，含月面曲率修正）", fontsize=18)
    if len(records) == 1:
        site_label = SITE_LABELS.get(
            records[0]["site_name"], records[0]["site_name"]
        )
        title = f"沙克尔顿坑{site_label}地形遮挡角剖面（极射投影，50 km 射线）"
    else:
        title = "沙克尔顿坑地形遮挡角剖面（极射投影，50 km 射线）"
    ax.set_title(title, fontsize=20)
    ax.set_xlim(0, 360)
    ax.set_ylim(y_bottom, y_top)
    ax.set_xticks(range(0, 361, 45))
    ax.grid(True, alpha=0.3)
    ax.legend(loc="center right", fontsize=18)
    ax.tick_params(axis='both', labelsize=18)
    fig.tight_layout()
    fname = (
        "horizon_mask_profile.png"
        if suffix is None
        else f"horizon_mask_{suffix}_profile.png"
    )
    path = os.path.join(output_dir, fname)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_polar(
    records: list[dict],
    output_dir: str,
    earth_range: list[float] | None,
    suffix: str | None = None,
) -> str:
    """极坐标地平线图（北=0° 顶部，顺时针），含地球高度角带（不含中继卫星）。

    单站点模式（如 ``--sites center``）输出带站点角标的
    ``horizon_mask_{suffix}_polar.png``；多站点仍为 ``horizon_mask_polar.png``。
    """
    fig, ax = plt.subplots(
        figsize=(7.5, 7.5), subplot_kw={"projection": "polar"}
    )
    colors = ["#1f77b4", "#d62728", "#2ca02c"]

    for rec, color in zip(records, colors):
        az = np.radians(np.array(rec["azimuths_deg"]))
        hor = np.array(rec["horizon_elev_deg"])
        label = SITE_LABELS.get(rec["site_name"], rec["site_name"])
        ax.plot(az, hor, color=color, lw=1.4, label=label)
    if earth_range is not None:
        lo, hi = earth_range
        az_full = np.linspace(0, 2 * np.pi, 361)
        ax.fill_between(
            az_full, lo, hi, color="tab:orange", alpha=0.22,
            label=f"地球高度角范围 [{lo:.1f}°, {hi:.1f}°]",
        )
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)  # 顺时针
    ax.set_rlabel_position(22)
    if len(records) == 1:
        site_label = SITE_LABELS.get(
            records[0]["site_name"], records[0]["site_name"]
        )
        title = f"沙克尔顿坑{site_label}地形遮挡角极坐标图（含地球高度角带）"
    else:
        title = "沙克尔顿坑地形遮挡角极坐标图（含地球高度角带）"
    ax.set_title(title, fontsize=20)
    ax.tick_params(axis='both', labelsize=18)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.16), fontsize=18, ncol=3)
    fig.tight_layout()
    fname = (
        "horizon_mask_polar.png"
        if suffix is None
        else f"horizon_mask_{suffix}_polar.png"
    )
    path = os.path.join(output_dir, fname)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--dem", default=DEM_PATH, help="DEM GeoTIFF 路径")
    p.add_argument("--lat", type=float, default=SHACKLETON_LAT_DEG,
                   help="观测点纬度（度），默认权威沙克尔顿坐标")
    p.add_argument("--lon", type=float, default=SHACKLETON_LON_DEG,
                   help="观测点经度（度）")
    p.add_argument("--max-radius-m", type=float, default=50000.0,
                   help="最大射线距离 / ROI 半径（米），默认 50000")
    p.add_argument("--n-az", type=int, default=360, help="方位角采样数，默认 360")
    p.add_argument("--n-samples", type=int, default=300,
                   help="每方位角射线采样点数，默认 300")
    p.add_argument("--roi-size", type=int, default=4096,
                   help="ROI 窗口降采样目标像素数，默认 4096")
    p.add_argument("--sites", default="center,north1km,south1km",
                   help="逗号分隔的观测点后缀子集（center/north1km/south1km），"
                        "默认全部三个；仅选一个时图表输出为带站点角标的单站文件")
    p.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR,
                   help="输出目录，默认 figures-and-logs")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    if not os.path.exists(args.dem):
        print(f"[horizon_mask] ERROR: DEM not found: {args.dem}")
        sys.exit(2)

    x_center, y_center = lonlat_to_polar_stereo_xy(args.lon, args.lat)
    print(
        f"[horizon_mask] center polar-stereo: ({x_center:.1f}, {y_center:.1f}) m; "
        f"DEM={args.dem}"
    )

    earth_range = load_earth_elev_range(args.output_dir)
    if earth_range is not None:
        print(f"[horizon_mask] loaded Earth elev range: {earth_range}")
    else:
        print("[horizon_mask] no Earth-elevation JSON yet; skipping cross-comparison")

    # --sites 选择观测点子集；保持 SITES 定义顺序（图例颜色与站点一致）
    site_by_suffix = {site[2]: (site[0], site[1]) for site in SITES}
    requested = [s.strip() for s in args.sites.split(",") if s.strip()]
    unknown = [s for s in requested if s not in site_by_suffix]
    if unknown:
        print(f"[horizon_mask] ERROR: unknown site suffix(es): {unknown}; "
              f"valid: {sorted(site_by_suffix)}")
        sys.exit(2)
    selected = [site[2] for site in SITES if site[2] in requested]

    records: list[dict] = []
    for suffix in selected:
        site_name, y_offset = site_by_suffix[suffix]
        # 站点沿极射投影 +y/-y 平移（极点附近不可用经纬度加减）
        x_site, y_site = x_center, y_center + y_offset
        # 平移后的经纬度（用于 JSON 记录/图表标注）
        from analysis.elevation_video.utils import polar_stereo_xy_to_lonlat
        lon_site, lat_site = polar_stereo_xy_to_lonlat(x_site, y_site)
        print(f"[horizon_mask] computing {site_name} at "
              f"stereo ({x_site:.1f}, {y_site:.1f}) ...")
        result = compute_horizon_mask(
            args.dem, x_site, y_site,
            max_radius_m=args.max_radius_m,
            n_az=args.n_az,
            n_samples=args.n_samples,
            roi_size=args.roi_size,
        )
        rec = make_json_record(site_name, lon_site, lat_site, result, earth_range)
        out_path = os.path.join(args.output_dir, f"horizon_mask_{suffix}.json")
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(rec, f, indent=2, ensure_ascii=False)
        records.append(rec)
        print(
            f"[horizon_mask]   {site_name}: max={rec['max_elev_deg']:.2f}°, "
            f"mean={rec['mean_elev_deg']:.2f}°, "
            f"az_of_max={rec['azimuth_of_max_deg']:.0f}°, "
            f"truncated_az={rec['n_truncated_az']}"
        )

    # 单站点模式：图表输出带站点角标（如 horizon_mask_center_profile.png）
    plot_suffix = selected[0] if len(selected) == 1 else None
    png1 = plot_profile(records, args.output_dir, earth_range, suffix=plot_suffix)
    png2 = plot_polar(records, args.output_dir, earth_range, suffix=plot_suffix)
    print(f"[horizon_mask] wrote {png1}")
    print(f"[horizon_mask] wrote {png2}")
    print("[horizon_mask] DONE")


if __name__ == "__main__":
    main()
