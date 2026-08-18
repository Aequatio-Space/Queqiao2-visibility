"""P0 科学分析共享工具函数。

提供四项核心能力：
- 南极极射投影坐标转换（与 ``illumination._lonlat_to_dem_xy`` 一致）
- 通用 ENU 高度角计算（复用 ``_spice_sun_vector`` 几何，目标天体可替换）
- 球面几何地平线仰角（含曲率修正）
- SPICE 内核加载 / 释放

所有函数为纯函数或模块级单例，无副作用；分析脚本通过 import 复用。
"""

from __future__ import annotations

import math
import os
import sys
from datetime import datetime, timedelta, timezone

import numpy as np

# 权威物理常数
MOON_R_EQ_M: float = 1737400.0  # 月球赤道半径（m，与 DEM 基准面 D_Moon_2000 一致）
MOON_R_KM: float = 1737.4  # 月球平均半径（km，鹊桥轨道计算用）
MU_MOON_KM3_S2: float = 4902.800066  # 月球引力常数（km³/s²）

# 权威沙克尔顿坐标（IAU/Gazetteer）
SHACKLETON_LAT_DEG: float = -89.67
SHACKLETON_LON_DEG: float = 129.78

# 仓库根（本文件位于 analysis/elevation_video/utils.py → 上溯三级）
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# SPICE 内核清单（相对仓库根目录）
_SPICE_KERNEL_REL: tuple[str, ...] = (
    "core/environment/data/latest_leapseconds.tls",
    "core/environment/data/moon_pa_de421_1900-2050.bpc",
    "core/environment/data/pck00011_n0066.tpc",
    "core/environment/data/de442s.bsp",
)
SPICE_KERNELS: tuple[str, ...] = tuple(
    os.environ.get("SPICE_KERNEL_DIR", _REPO_ROOT) + "/" + rel
    for rel in _SPICE_KERNEL_REL
)

# DEM 路径：必须通过环境变量 DEM_PATH 提供（本地路径因人而异，不硬编码）。
# 未设置时为空字符串，DEM 相关端点/分析会明确报错提示（见 viewer/server.py 的 503 分支）。
DEM_PATH: str = os.environ.get("DEM_PATH", "")

# MJD 参考历元（1858-11-17 00:00 UTC）
_MJD_EPOCH = datetime(1858, 11, 17, tzinfo=timezone.utc)


# ---------------------------------------------------------------------------
# 极射投影
# ---------------------------------------------------------------------------
def lonlat_to_polar_stereo_xy(
    lon_deg: float, lat_deg: float, R: float = MOON_R_EQ_M
) -> tuple[float, float]:
    """经纬度 → 南极极射投影 (x, y) 坐标（米）。

    与 ``illumination.py::_lonlat_to_dem_xy`` 公式一致：
      colat = radians(90 + lat), r = 2R·tan(colat/2)
      x = r·sin(lon), y = r·cos(lon)

    +x = 东（lon=90° 方向），+y = 北（赤道方向，lon=0 方向）。
    南极点 (lon, -90) → (0, 0)。
    """
    colat = math.radians(90.0 + lat_deg)
    r = 2.0 * R * math.tan(colat / 2.0)
    lon_rad = math.radians(lon_deg)
    return r * math.sin(lon_rad), r * math.cos(lon_rad)


def polar_stereo_xy_to_lonlat(
    x: float, y: float, R: float = MOON_R_EQ_M
) -> tuple[float, float]:
    """南极极射投影 (x, y) → 经纬度（度）。

    极射投影逆变换：
      r = hypot(x, y), colat = 2·atan(r / 2R)
      lat = colat_deg - 90, lon = degrees(atan2(x, y))
    """
    r = math.hypot(x, y)
    colat = 2.0 * math.atan2(r, 2.0 * R)
    lat_deg = math.degrees(colat) - 90.0
    lon_deg = math.degrees(math.atan2(x, y))
    return lon_deg, lat_deg


# ---------------------------------------------------------------------------
# 曲率修正
# ---------------------------------------------------------------------------
def spherical_horizon_angle(
    d_m: float, h0_m: float, h_m: float, R: float = MOON_R_EQ_M
) -> float:
    """球面几何地平线仰角（含曲率修正）。

    :param d_m:  月面距离（米，极射投影平面上的水平距离）
    :param h0_m: 观测者高程（米，参考椭球面以上）
    :param h_m:  目标点高程（米，参考椭球面以上）
    :param R:    月球半径（米）
    :return: 仰角（度）

    退化验证：d → 0 时退化为 flat-Earth 公式 atan2(h - h0, d)。
    数值验证：d=50 km, Δh=1000 m 时，曲率下降 ~720 m，修正后与
    flat-Earth 差异约 0.71°（不可忽略）。
    """
    theta = d_m / R  # 中心角（弧度）
    r_obs = R + h0_m  # 观测者到月心距离
    r_tgt = R + h_m  # 目标点到月心距离
    elevation = math.atan2(
        r_tgt * math.cos(theta) - r_obs,
        r_tgt * math.sin(theta),
    )
    return math.degrees(elevation)


def spherical_horizon_angle_vec(
    d_m: np.ndarray, h0_m: float, h_m: np.ndarray, R: float = MOON_R_EQ_M
) -> np.ndarray:
    """向量化球面几何地平线仰角（弧度，用于批量地平线射线采样）。"""
    theta = d_m / R
    r_obs = R + h0_m
    r_tgt = R + h_m
    return np.degrees(
        np.arctan2(r_tgt * np.cos(theta) - r_obs, r_tgt * np.sin(theta))
    )


# ---------------------------------------------------------------------------
# ENU 高度角
# ---------------------------------------------------------------------------
def _enu_basis(
    lon_deg: float, lat_deg: float, R_eq: float = MOON_R_EQ_M, R_pol: float | None = None
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """构造局部 ENU 正交基 (east, north, up)。

    与 ``illumination.py::_spice_sun_vector`` 完全一致：
    - Up = 月球椭球面外法线（非地心径向）
    - East = 纬圈切向，Gram-Schmidt 正交化
    - North = up × east，构成右手系
    在极点附近稳健（Gram-Schmidt 修正）。
    """
    if R_pol is None:
        R_pol = R_eq
    lon_rad = math.radians(lon_deg)
    lat_rad = math.radians(lat_deg)

    cos_lat, sin_lat = math.cos(lat_rad), math.sin(lat_rad)
    cos_lon, sin_lon = math.cos(lon_rad), math.sin(lon_rad)

    up_raw = np.array([cos_lat * cos_lon / R_eq, cos_lat * sin_lon / R_eq, sin_lat / R_pol])
    up = up_raw / np.linalg.norm(up_raw)

    east_raw = np.array([-sin_lon, cos_lon, 0.0])
    east = east_raw - np.dot(east_raw, up) * up
    east /= np.linalg.norm(east)

    north = np.cross(up, east)
    north /= np.linalg.norm(north)
    return east, north, up


def enu_elevation_angle(
    target_pos_iau_moon_m: np.ndarray,
    lon_deg: float,
    lat_deg: float,
    h_m: float,
    R_eq: float = MOON_R_EQ_M,
    R_pol: float | None = None,
) -> float:
    """通用 ENU 高度角（度），目标天体在 IAU_MOON 月固系中的位置 [x,y,z]（米）。

    参考面为月球参考椭球，不是坑底局部斜面。
    观测者位置按 ``r = R_eq + h_m`` 的球近似计算（LOLA DEM 约定）。

    :param target_pos_iau_moon_m: 目标在 IAU_MOON 系的位置（米）
    :param lon_deg: 观测点经度（度）
    :param lat_deg: 观测点纬度（度）
    :param h_m: 观测点高程（米，参考椭球面以上）
    :return: 高度角（度），>0 表示在地平线上方
    """
    if R_pol is None:
        R_pol = R_eq
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
    vec = np.asarray(target_pos_iau_moon_m, dtype=float) - obs

    east, north, up = _enu_basis(lon_deg, lat_deg, R_eq, R_pol)
    u_comp = float(np.dot(vec, up))
    e_comp = float(np.dot(vec, east))
    n_comp = float(np.dot(vec, north))
    horizontal = math.hypot(e_comp, n_comp)
    return math.degrees(math.atan2(u_comp, horizontal))


def enu_elevation_azimuth(
    target_pos_iau_moon_m: np.ndarray,
    lon_deg: float,
    lat_deg: float,
    h_m: float,
    R_eq: float = MOON_R_EQ_M,
    R_pol: float | None = None,
) -> tuple[float, float]:
    """通用 ENU 高度角 + 方位角（度），目标天体在 IAU_MOON 月固系中的位置 [x,y,z]（米）。

    与 :func:`enu_elevation_angle` 使用完全相同的观测者位置与局部 ENU 基，
    额外返回方位角（北=0°，顺时针增加，与地平线图约定一致）。

    :param target_pos_iau_moon_m: 目标在 IAU_MOON 系的位置（米）
    :param lon_deg: 观测点经度（度）
    :param lat_deg: 观测点纬度（度）
    :param h_m: 观测点高程（米，参考椭球面以上）
    :return: (elevation_deg, azimuth_deg) — 高度角 >0 表示在地平线上方；
        方位角 0°=北、90°=东、180°=南、270°=西
    """
    if R_pol is None:
        R_pol = R_eq
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
    vec = np.asarray(target_pos_iau_moon_m, dtype=float) - obs

    east, north, up = _enu_basis(lon_deg, lat_deg, R_eq, R_pol)
    u_comp = float(np.dot(vec, up))
    e_comp = float(np.dot(vec, east))
    n_comp = float(np.dot(vec, north))
    horizontal = math.hypot(e_comp, n_comp)
    elevation = math.degrees(math.atan2(u_comp, horizontal))
    # 方位角：atan2(east, north)，0°=北、90°=东、180°=南、270°=西
    azimuth = math.degrees(math.atan2(e_comp, n_comp)) % 360.0
    return elevation, azimuth


# ---------------------------------------------------------------------------
# SPICE 加载 / 释放
# ---------------------------------------------------------------------------
def load_spice_kernels() -> list[str]:
    """加载本地 SPICE 内核（4 个文件），返回已加载路径列表。

    幂等：重复调用 ``furnsh`` 同一内核会覆盖加载，无副作用。
    """
    import spiceypy

    for k in SPICE_KERNELS:
        spiceypy.furnsh(k)
    return list(SPICE_KERNELS)


def cleanup_spice() -> None:
    """释放全部 SPICE 内核。"""
    import spiceypy

    spiceypy.kclear()


def moon_radii_m() -> tuple[float, float]:
    """从 PCK 内核读取月球椭球半径（赤道, 极），单位米。"""
    import spiceypy

    _, radii_km = spiceypy.bodvrd("MOON", "RADII", 3)
    return radii_km[0] * 1000.0, radii_km[2] * 1000.0


# ---------------------------------------------------------------------------
# 时间转换
# ---------------------------------------------------------------------------
def mjd_to_utc_iso(t_mjd: float) -> str:
    """MJD → UTC ISO 字符串（与仿真器约定一致）。"""
    return (_MJD_EPOCH + timedelta(days=t_mjd)).isoformat()


def utc_iso_to_mjd(t_utc: str) -> float:
    """UTC ISO 字符串 → MJD。

    支持 ``"YYYY-MM-DD"``、``"YYYY-MM-DDTHH:MM:SS"``、
    ``"YYYY-MM-DDTHH:MM:SS+HH:MM"``（默认按 UTC 解释）。
    """
    s = t_utc.strip()
    if "T" in s or " " in s:
        # 处理空格分隔形式
        s = s.replace(" ", "T", 1)
    dt = datetime.fromisoformat(s)
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    dt_utc = dt.astimezone(timezone.utc)
    return (dt_utc - _MJD_EPOCH).total_seconds() / 86400.0


# ---------------------------------------------------------------------------
# 时间序列生成
# ---------------------------------------------------------------------------
def make_time_grid(
    start_mjd: float, end_mjd: float, interval_h: float
) -> np.ndarray:
    """生成 [start, end] 的均匀时间网格（MJD，浮点），步长 interval_h 小时。

    包含端点，点数 = floor((end-start)*24/interval) + 1。
    """
    n = int((end_mjd - start_mjd) * 24.0 / interval_h) + 1
    return start_mjd + np.arange(n, dtype=float) * (interval_h / 24.0)


# ---------------------------------------------------------------------------
# 中文字体（macOS 常见 CJK 字体，按优先级 fallback）
# ---------------------------------------------------------------------------
# 简体中文字体优先（与 queqiao_orbit_animation.py 的 _CJK_FONTS 一致）
_CJK_FONTS = [
    "Songti SC",
    "SimHei",
    "Heiti TC",
    "Arial Unicode MS",
    "Noto Sans CJK SC",
    "Hiragino Sans GB",
    "PingFang SC",
]


def setup_cjk_font() -> str | None:
    """配置 matplotlib 中文字体，返回选中的字体名（找不到返回 None）。

    供 2D 天空视角图/动画使用：标题、高度角、时间等中文标注需要 CJK 字体。
    幂等，可多次调用。
    """
    from matplotlib import font_manager as fm

    avail = {f.name for f in fm.fontManager.ttflist}
    chosen = next((c for c in _CJK_FONTS if c in avail), None)
    if chosen is None:
        print("警告：未找到中文字体，中文可能显示为方块", file=sys.stderr)
        return None
    plt = sys.modules.get("matplotlib.pyplot")
    if plt is not None:
        plt.rcParams["font.sans-serif"] = [chosen, "DejaVu Sans"]
        plt.rcParams["font.family"] = "sans-serif"
        plt.rcParams["axes.unicode_minus"] = False
    return chosen
