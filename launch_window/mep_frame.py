"""launch_window/mep_frame.py — MEP (Moon-Earth Perifocal) coordinate system (S-0).

Implements the MEP frame exactly as defined in the PRD v4 §S-0 (strictly per
the Yang Weilian paper):

| Item            | Definition                                                        |
|-----------------|-------------------------------------------------------------------|
| Origin          | Earth centre                                                      |
| Fundamental plane| Moon's orbit plane about Earth (intersects Earth & Moon spheres as their respective "equatorial planes", distinct from the true equator) |
| X axis          | Earth-Moon line, pointing from Moon to Earth (i.e. anti-Moon)     |
| Inertialisation | X-axis direction frozen at the instant of arrival at perilune     |
| Earth longitude | measured eastward from the X axis                                 |
| Lunar-surface longitude zero | the X-axis intersection with the lunar surface facing Earth |
| Latitude        | north-positive from the fundamental plane                         |

The frame is *inertialised*: the basis is computed once at ``t_arrive`` and
then held fixed, so that a state expressed in MEP coordinates is time
independent (the paper's decoupling relies on a circular lunar orbit; we do
NOT make that approximation — the Moon position is always the SPICE truth,
and the resulting loss of decoupling is handled explicitly by R-5).

Self-check (S-0 assertion): at ``t_arrive`` the Moon centre in MEP
coordinates must be at (lon, lat) = (180°, 0°) to within 1e-9.

All functions are pure functions of the loaded SPICE kernels (see
``geometry.ensure_kernels``).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from launch_window import geometry as g

# Moon mean radius (km) used for lunar-surface lon/lat conversions.  The PRD
# uses the SPICE PCK radii for the ellipsoid; we use the mean radius for the
# spherical lon/lat mapping (the paper's MEP lunar-surface frame is spherical).
MOON_MEAN_RADIUS_KM = 1737.4


@dataclass(frozen=True)
class MEPFrame:
    """An inertialised MEP basis frozen at a given arrival epoch.

    Attributes
    ----------
    t_arrive : float
        SPICE ephemeris time (ET) at which the frame is frozen.
    x, y, z : (3,) ndarray
        Orthonormal basis vectors (J2000 components) of the MEP frame.
    """

    t_arrive: float
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray

    @property
    def rot_j2000_to_mep(self) -> np.ndarray:
        """3x3 rotation matrix mapping J2000 vectors to MEP coordinates.

        Rows are the MEP basis vectors expressed in J2000, so
        ``v_mep = R @ v_j2000``.
        """
        return np.vstack([self.x, self.y, self.z])

    @property
    def rot_mep_to_j2000(self) -> np.ndarray:
        """3x3 rotation matrix mapping MEP vectors to J2000 (inverse)."""
        return self.rot_j2000_to_mep.T


def mep_at_epoch(t_arrive: float) -> MEPFrame:
    """Build the inertialised MEP frame frozen at arrival epoch ``t_arrive``.

    Parameters
    ----------
    t_arrive : float
        SPICE ephemeris time (ET) of arrival at perilune.

    Returns
    -------
    MEPFrame
        The frozen MEP basis.

    Notes
    -----
    - Fundamental-plane normal ``z`` = ``r_moon × v_moon`` (J2000, SPICE).
    - X axis = ``-r_hat_moon(t_arrive)`` (anti-Moon direction).
    - Y axis completes the right-handed system.
    """
    g.ensure_kernels()
    import spiceypy

    state, _lt = spiceypy.spkezr("MOON", t_arrive, "J2000", "NONE", "EARTH")
    r_moon = np.asarray(state[:3], dtype=float)
    v_moon = np.asarray(state[3:6], dtype=float)

    # Fundamental-plane normal (Moon's orbital angular momentum direction).
    z = np.cross(r_moon, v_moon)
    z = z / np.linalg.norm(z)

    # X axis: anti-Moon direction.
    x = -r_moon / np.linalg.norm(r_moon)

    # Y axis completes the right-handed system.
    y = np.cross(z, x)
    y = y / np.linalg.norm(y)

    return MEPFrame(t_arrive=t_arrive, x=x, y=y, z=z)


def to_mep(r_j2000: np.ndarray, frame: MEPFrame) -> np.ndarray:
    """Convert a J2000 position/velocity vector to MEP coordinates."""
    return frame.rot_j2000_to_mep @ np.asarray(r_j2000, dtype=float)


def from_mep(r_mep: np.ndarray, frame: MEPFrame) -> np.ndarray:
    """Convert a MEP position/velocity vector to J2000 coordinates."""
    return frame.rot_mep_to_j2000 @ np.asarray(r_mep, dtype=float)


# ---------------------------------------------------------------------------
# Earth-side lon/lat (origin at Earth centre)
# ---------------------------------------------------------------------------
def to_mep_lonlat(r_mep: np.ndarray) -> tuple[float, float]:
    """Earth-side MEP (lon, lat) in degrees from a MEP position vector.

    Longitude is measured eastward from the +X axis; latitude is
    north-positive from the fundamental (XY) plane.
    """
    x, y, z = (float(v) for v in r_mep)
    lon = math.degrees(math.atan2(y, x))
    lat = math.degrees(math.asin(z / math.sqrt(x * x + y * y + z * z)))
    return lon, lat


def from_mep_lonlat(lon_deg: float, lat_deg: float, r: float) -> np.ndarray:
    """MEP position vector from Earth-side (lon, lat, radius)."""
    lon = math.radians(lon_deg)
    lat = math.radians(lat_deg)
    cos_lat = math.cos(lat)
    return np.array(
        [
            r * cos_lat * math.cos(lon),
            r * cos_lat * math.sin(lon),
            r * math.sin(lat),
        ]
    )


# ---------------------------------------------------------------------------
# Lunar-surface side lon/lat (origin at Moon centre, same X axis, zero
# longitude at the sub-Earth point)
# ---------------------------------------------------------------------------
def to_mep_lunar_lonlat(r_mep: np.ndarray) -> tuple[float, float]:
    """Lunar-surface MEP (lon, lat) in degrees from a Moon-centred MEP vector.

    The origin is the Moon centre; the X axis is the same anti-Moon direction
    and the zero-longitude meridian is the sub-Earth point (the X-axis
    intersection facing Earth).  Because the Moon is tidally locked, the MEP
    lunar-surface longitude ≈ the Moon-fixed (IAU_MOON) longitude, ignoring
    libration (PRD note).
    """
    x, y, z = (float(v) for v in r_mep)
    lon = math.degrees(math.atan2(y, x))
    lat = math.degrees(math.asin(z / math.sqrt(x * x + y * y + z * z)))
    return lon, lat


def from_mep_lunar_lonlat(lon_deg: float, lat_deg: float, r: float) -> np.ndarray:
    """Moon-centred MEP position vector from lunar-surface (lon, lat, radius)."""
    lon = math.radians(lon_deg)
    lat = math.radians(lat_deg)
    cos_lat = math.cos(lat)
    return np.array(
        [
            r * cos_lat * math.cos(lon),
            r * cos_lat * math.sin(lon),
            r * math.sin(lat),
        ]
    )


def assert_moon_at_anti_x(frame: MEPFrame, tol: float = 1e-9) -> None:
    """S-0 self-check: the Moon centre in MEP is at (lon, lat) = (180°, 0°).

    Raises
    ------
    AssertionError
        If the Moon centre is not at (180°, 0°) within ``tol`` degrees.
    """
    g.ensure_kernels()
    import spiceypy

    state, _lt = spiceypy.spkezr("MOON", frame.t_arrive, "J2000", "NONE", "EARTH")
    r_moon_j2000 = np.asarray(state[:3], dtype=float)
    r_moon_mep = to_mep(r_moon_j2000, frame)
    lon, lat = to_mep_lonlat(r_moon_mep)
    # The Moon centre is at the anti-X direction: lon = 180°, lat = 0°.
    lon_err = abs(((lon - 180.0 + 180.0) % 360.0) - 180.0)
    assert lon_err < tol, f"MEP self-check failed: Moon lon={lon:.3e} (want 180)"
    assert abs(lat) < tol, f"MEP self-check failed: Moon lat={lat:.3e} (want 0)"
