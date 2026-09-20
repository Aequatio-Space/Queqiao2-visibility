"""launch_window/geometry.py — shared SPICE geometry for the launch-window engine.

Pure-ephemeris quantities required by ``validate.py`` (B-7/B-8/B-9) and
``predict.py`` (P-2/P-3/P-5).  No calibration-derived constants live here;
they belong in ``calibration.json`` (see ``anti_cheat.py``).

Kernel note
-----------
The PRD lists ``de440.bsp`` / ``moon_pa_de440.bpc``; the repo ships
``de442s.bsp`` (1900-2050) and ``moon_pa_de421_1900-2050.bpc`` which cover
the same era and are sufficient for CE-7 (2026).  Swapping in de440 requires
no code change — only the kernel list below.

All functions are pure functions of the loaded SPICE kernels.  SPICE is not
thread-safe; these scripts are single-threaded CLIs, so no lock is needed.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable

import numpy as np

# 内核路径统一由仓库内的 analysis/elevation_video/utils.py 提供
# （默认 core/environment/data/，可用 SPICE_KERNEL_DIR 覆盖），
# 与 viewer/server.py 使用同一份内核清单，避免路径逻辑重复。
from analysis.elevation_video import utils as _u

# ---------------------------------------------------------------------------
# Kernel management
# ---------------------------------------------------------------------------
_REPO_ROOT = Path(__file__).resolve().parent.parent

KERNELS: tuple[str, ...] = tuple(_u.SPICE_KERNELS)

_loaded = False


def ensure_kernels(kernels: Iterable[str] = KERNELS) -> None:
    """Furnsh the SPICE kernels once (idempotent)."""
    global _loaded
    if _loaded:
        return
    import spiceypy

    missing = [k for k in kernels if not Path(k).exists()]
    if missing:
        raise FileNotFoundError(f"Missing SPICE kernels: {missing}")
    spiceypy.furnsh(list(kernels))
    _loaded = True


def kclear() -> None:
    """Unload all kernels (mainly for tests)."""
    global _loaded
    import spiceypy

    spiceypy.kclear()
    _loaded = False


# ---------------------------------------------------------------------------
# Time helpers
# ---------------------------------------------------------------------------
def utc_to_et(utc_str: str) -> float:
    """UTC ISO string -> SPICE ephemeris time (seconds past J2000)."""
    import spiceypy

    ensure_kernels()
    return spiceypy.str2et(utc_str)


def mjd_to_et(mjd: float) -> float:
    """MJD (UTC) -> SPICE ephemeris time."""
    import spiceypy

    ensure_kernels()
    return spiceypy.unitim(mjd + 2400000.5, "JDTDB", "ET")


def et_to_utc(et: float) -> str:
    """SPICE ephemeris time -> UTC ISO string (seconds precision)."""
    import spiceypy

    ensure_kernels()
    return spiceypy.et2utc(et, "ISOC", 0)


def et_to_mjd(et: float) -> float:
    """SPICE ephemeris time -> MJD (UTC)."""
    import spiceypy

    ensure_kernels()
    return spiceypy.unitim(et, "ET", "JDTDB") - 2400000.5


# ---------------------------------------------------------------------------
# Moon geometry
# ---------------------------------------------------------------------------
def moon_geocentric_declination_deg(et: float) -> float:
    """Geocentric declination of the Moon in the J2000 (Earth equatorial) frame.

    This is the quantity B-7 calibrates: the Moon's latitude above/below the
    Earth equatorial plane, which bounds the reachable parking-orbit
    inclination from Wenchang (19.6 N).
    """
    import spiceypy

    ensure_kernels()
    state, _lt = spiceypy.spkezr("MOON", et, "J2000", "NONE", "EARTH")
    r = float(np.linalg.norm(state[:3]))
    return math.degrees(math.asin(float(state[2]) / r))


def moon_geocentric_declination_rate_deg_per_day(et: float, dt_d: float = 1.0) -> float:
    """Central-difference rate of the Moon's geocentric declination (deg/day)."""
    dt_et = dt_d * 86400.0
    d0 = moon_geocentric_declination_deg(et - 0.5 * dt_et)
    d1 = moon_geocentric_declination_deg(et + 0.5 * dt_et)
    return (d1 - d0) / dt_d


def _moon_radii_m() -> tuple[float, float]:
    """(equatorial, polar) radius of the Moon in metres from the PCK kernel."""
    import spiceypy

    ensure_kernels()
    _, radii_km = spiceypy.bodvrd("MOON", "RADII", 3)
    return radii_km[0] * 1000.0, radii_km[2] * 1000.0


def sub_earth_lonlat_deg(et: float) -> tuple[float, float]:
    """(lon, lat) of the sub-Earth point on the Moon in the IAU_MOON frame.

    Longitude is east-positive, wrapped to [0, 360).  The sub-Earth meridian
    is the "12 o'clock" direction of the tidally-locked Moon; the lunar hour
    angle at a site is measured from it.
    """
    import spiceypy

    ensure_kernels()
    state, _lt = spiceypy.spkezr("EARTH", et, "IAU_MOON", "LT+S", "MOON")
    x, y, z = (float(v) for v in state[:3])
    lon = math.degrees(math.atan2(y, x)) % 360.0
    lat = math.degrees(math.asin(z / math.sqrt(x * x + y * y + z * z)))
    return lon, lat


def sub_solar_lonlat_deg(et: float) -> tuple[float, float]:
    """(lon, lat) of the sub-solar point on the Moon in the IAU_MOON frame.

    The sub-solar latitude drives the polar "season" (P-2): the south pole is
    preferentially lit when the sub-solar latitude is negative.
    """
    import spiceypy

    ensure_kernels()
    state, _lt = spiceypy.spkezr("SUN", et, "IAU_MOON", "LT+S", "MOON")
    x, y, z = (float(v) for v in state[:3])
    lon = math.degrees(math.atan2(y, x)) % 360.0
    lat = math.degrees(math.asin(z / math.sqrt(x * x + y * y + z * z)))
    return lon, lat


def lunar_hour_angle_deg(site_lon_deg: float, et: float) -> float:
    """Hour angle of the Moon at a landing site (deg, [0, 360)).

    Defined as the east-positive angular separation from the site meridian to
    the sub-Earth meridian: ``(site_lon - sub_earth_lon) mod 360``.  This is
    the longitude difference B-8 correlates with parking-orbit duration
    (CE-4's official 270 deg / 13.176 deg/day = 20.5 d).
    """
    sub_earth_lon, _ = sub_earth_lonlat_deg(et)
    return (site_lon_deg - sub_earth_lon) % 360.0


# ---------------------------------------------------------------------------
# Local ENU frame (matches core/environment/illumination.py conventions)
# ---------------------------------------------------------------------------
def _enu_basis(
    lon_deg: float, lat_deg: float, r_eq: float, r_pol: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Local (east, north, up) orthonormal basis at a lunar surface point.

    Up is the ellipsoid outward normal (not the geocentric radial direction);
    east is the latitude-circle tangent, Gram-Schmidt orthogonalised against
    up; north completes the right-handed system.  Robust near the poles.
    """
    lon_rad = math.radians(lon_deg)
    lat_rad = math.radians(lat_deg)
    cos_lat, sin_lat = math.cos(lat_rad), math.sin(lat_rad)
    cos_lon, sin_lon = math.cos(lon_rad), math.sin(lon_rad)

    up_raw = np.array(
        [cos_lat * cos_lon / r_eq, cos_lat * sin_lon / r_eq, sin_lat / r_pol]
    )
    up = up_raw / np.linalg.norm(up_raw)

    east_raw = np.array([-sin_lon, cos_lon, 0.0])
    east = east_raw - np.dot(east_raw, up) * up
    east /= np.linalg.norm(east)

    north = np.cross(up, east)
    north /= np.linalg.norm(north)
    return east, north, up


def sun_elev_az_deg(
    lon_deg: float, lat_deg: float, h_m: float, et: float
) -> tuple[float, float]:
    """Sun (elevation, azimuth) at a lunar surface point in the IAU_MOON frame.

    Elevation > 0 means the Sun is above the local horizon (flat-Moon
    geometry, no terrain).  Azimuth is 0 deg = North, 90 = East, 180 = South,
    270 = West (clockwise from North), matching the horizon-mask convention.
    """
    import spiceypy

    ensure_kernels()
    r_eq, r_pol = _moon_radii_m()
    lon_rad = math.radians(lon_deg)
    lat_rad = math.radians(lat_deg)
    radius = r_eq + h_m
    obs = np.array(
        [
            radius * math.cos(lat_rad) * math.cos(lon_rad),
            radius * math.cos(lat_rad) * math.sin(lon_rad),
            radius * math.sin(lat_rad),
        ]
    )

    state, _lt = spiceypy.spkezr("SUN", et, "IAU_MOON", "LT+S", "MOON")
    sun = np.asarray(state[:3], dtype=float) * 1000.0  # km -> m
    vec = sun - obs

    east, north, up = _enu_basis(lon_deg, lat_deg, r_eq, r_pol)
    u = float(np.dot(vec, up))
    e = float(np.dot(vec, east))
    n = float(np.dot(vec, north))
    elev = math.degrees(math.atan2(u, math.hypot(e, n)))
    az = math.degrees(math.atan2(e, n)) % 360.0
    return elev, az


def sun_elevation_deg(lon_deg: float, lat_deg: float, h_m: float, et: float) -> float:
    """Sun elevation only (deg) — convenience wrapper."""
    return sun_elev_az_deg(lon_deg, lat_deg, h_m, et)[0]


# ---------------------------------------------------------------------------
# Parking-orbit geometry (B-8)
# ---------------------------------------------------------------------------
def parking_longitude_diff_deg(
    site_lon_deg: float, et_oi: float, et_land: float
) -> dict[str, float]:
    """Longitude-difference quantities for the B-8 parking-duration model.

    The PRD hypothesis is ``parking_days ~= longitude_diff / 13.176`` where
    the longitude difference is between the lunar hour angle at orbit
    insertion and the landing-site longitude.  This function computes every
    candidate definition so the engine can test the hypothesis against the
    four historical missions (the PRD explicitly forbids trusting intuition).

    Definitions (all in degrees, east-positive, wrapped to [0, 360)):
    - ``site_lon``            : landing-site longitude (fixed in Moon frame)
    - ``sub_earth_lon_oi``    : sub-Earth longitude at orbit insertion
    - ``sub_solar_lon_oi``    : sub-solar longitude at orbit insertion
    - ``hour_angle_earth``    : site_lon - sub_earth_lon at OI
    - ``hour_angle_solar``    : site_lon - sub_solar_lon at OI
    - ``rotated_angle``       : site_lon - sub_earth_lon at landing
    - ``site_rotation_during_parking`` : parking_days * 13.176 (the angle the
      Moon rotates while parked; the site's inertial longitude advances by this)
    """
    se_oi, _ = sub_earth_lonlat_deg(et_oi)
    ss_oi, _ = sub_solar_lonlat_deg(et_oi)
    se_land, _ = sub_earth_lonlat_deg(et_land)
    park_d = (et_land - et_oi) / 86400.0
    return {
        "site_lon": site_lon_deg % 360.0,
        "sub_earth_lon_oi": se_oi,
        "sub_solar_lon_oi": ss_oi,
        "hour_angle_earth": (site_lon_deg - se_oi) % 360.0,
        "hour_angle_solar": (site_lon_deg - ss_oi) % 360.0,
        "hour_angle_earth_landing": (site_lon_deg - se_land) % 360.0,
        "site_rotation_during_parking": park_d * 13.176,
    }


# ---------------------------------------------------------------------------
# Transfer-plane geometry (PRD v2 — lambert_window.py)
# ---------------------------------------------------------------------------
def moon_position_j2000_km(et: float) -> np.ndarray:
    """Geocentric Moon position vector (km) in the J2000 frame.

    This is the arrival point of the TLI transfer leg.  Used to define the
    transfer plane spanned by the Earth centre, the TLI point, and the Moon.
    """
    import spiceypy

    ensure_kernels()
    state, _lt = spiceypy.spkezr("MOON", et, "J2000", "NONE", "EARTH")
    return np.asarray(state[:3], dtype=float)


def launch_site_vector_j2000_km(
    lon_deg: float, lat_deg: float, alt_m: float, et: float
) -> np.ndarray:
    """Geocentric position vector (km) of a launch site in the J2000 frame.

    Converts the site's geodetic (lon, lat, alt) coordinates from the
    Earth-fixed (ITRF) frame to the inertial J2000 frame at epoch ``et``.
    This is the vector whose presence in the transfer plane is the
    minute-level "crossing condition" (PRD v2 condition B).
    """
    import spiceypy

    ensure_kernels()
    # Geodetic -> geocentric rectangular in the Earth-fixed frame (ITRF).
    # SPICE's georec uses the IAU_EARTH body-fixed frame by default; we use
    # the WGS84 ellipsoid via the PCK radii.
    re, rp = _earth_radii_km()
    # Geodetic to geocentric (approximate, using the ellipsoid)
    lat_rad = math.radians(lat_deg)
    lon_rad = math.radians(lon_deg)
    e2 = 1.0 - (rp / re) ** 2
    N = re / math.sqrt(1.0 - e2 * math.sin(lat_rad) ** 2)
    x = (N + alt_m / 1000.0) * math.cos(lat_rad) * math.cos(lon_rad)
    y = (N + alt_m / 1000.0) * math.cos(lat_rad) * math.sin(lon_rad)
    z = (N * (1.0 - e2) + alt_m / 1000.0) * math.sin(lat_rad)
    geodetic = np.array([x, y, z], dtype=float)

    # Rotate from the Earth-fixed frame to J2000 using SPICE's pxform.
    # The body-fixed frame for the Earth is "IAU_EARTH" (available from the
    # PCK kernel; ITRF93 would require a separate Earth-orientation PCK).
    rot = spiceypy.pxform("IAU_EARTH", "J2000", et)
    return rot @ geodetic


def _earth_radii_km() -> tuple[float, float]:
    """(equatorial, polar) radius of the Earth in km from the PCK kernel."""
    import spiceypy

    ensure_kernels()
    _, radii_km = spiceypy.bodvrd("EARTH", "RADII", 3)
    return float(radii_km[0]), float(radii_km[2])


def transfer_plane_normal(
    tli_pos_km: np.ndarray, moon_pos_km: np.ndarray
) -> np.ndarray:
    """Unit normal of the transfer plane spanned by the TLI point and the Moon.

    The transfer plane contains the Earth centre (origin), the TLI point, and
    the Moon's arrival position.  Its normal is ``cross(tli, moon)``
    normalised.  The launch site's geocentric vector must lie in this plane
    for the transfer to be reachable (PRD v2 condition B).
    """
    n = np.cross(np.asarray(tli_pos_km, dtype=float),
                 np.asarray(moon_pos_km, dtype=float))
    norm = float(np.linalg.norm(n))
    if norm < 1e-9:
        raise ValueError("TLI point and Moon position are collinear")
    return n / norm


def plane_inclination_deg(normal: np.ndarray) -> float:
    """Inclination (deg) of a plane given its unit normal.

    The inclination is the angle between the plane's normal and the Earth's
    equatorial (J2000 Z) axis, i.e. ``acos(|n_z|)``.  A plane containing the
    launch site must have inclination >= the site's latitude for the site to
    be reachable (PRD v2 condition A).
    """
    n = np.asarray(normal, dtype=float)
    return math.degrees(math.acos(min(1.0, abs(float(n[2])))))


def site_in_plane_deg(
    site_vec_km: np.ndarray, normal: np.ndarray
) -> float:
    """Angular separation (deg) of the launch-site vector from the transfer plane.

    The signed distance of the site vector from the plane is
    ``dot(site, normal)``; the angular separation is
    ``asin(dot(site_hat, normal))``.  A value near 0 means the site lies in
    the plane (the crossing condition is satisfied).
    """
    site = np.asarray(site_vec_km, dtype=float)
    n = np.asarray(normal, dtype=float)
    site_hat = site / float(np.linalg.norm(site))
    return math.degrees(math.asin(float(np.clip(np.dot(site_hat, n), -1.0, 1.0))))