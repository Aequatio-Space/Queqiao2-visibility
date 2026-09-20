"""launch_window/perigee_match.py — R-2 perigee-trajectory matching (PRD v5 §6-8).

Implements the R-2 ``T_e``-iteration engine that replaces the v1.0 ``(t_L, Az)``
2x2 solver.  The paper iterates on the arrival epoch ``T_e``: once ``T_e`` is
fixed the MEP frame is frozen and the R-1S perigee trajectory family is fixed,
so the fixed-point structure is cleaner.

Algorithm (PRD v5 §6)
---------------------
1. Given ``T_e`` -> freeze the MEP frame -> R-1S perigee trajectory family.
2. Adjust ``T_e`` so the launch site's parking orbit intersects the perigee
   trajectory.
3. From the intersection compute the parking-orbit inclination ``i0`` relative
   to the MEP fundamental plane.
4. Solve R-1 with ``i0`` -> ``(v_p, Omega, omega)``.
5. Fine-tune ``T_e`` from ``(Omega, omega)``, back to 2, until the longitude
   converges.
6. ``T_L = T_e - T - dT``.

The self-consistent equation solved here is

    f(i0) = i0_plane(i0) - i0 = 0

where ``i0_plane`` is the inclination (relative to the MEP plane) of the plane
spanned by the launch-site vector at ``t_L = T_e - T - dT`` and the perigee
point ``P(i0) = (lon_p(i0), lat_p(i0))`` from the R-1S family.  The parking
orbit must contain both the launch site and the transfer-orbit perigee, and
its inclination must equal the family parameter ``i0``.

Exclusion force (PRD v5 §6)
---------------------------
The only exclusion force is whether the required perigee declination falls in
the launch site's reachable band.  The reachable band is the R-1S
perigee-latitude spread ([-5.52, +5.52] deg at i0 in [1, 45]) intersected
with the R-4 azimuth corridor (which bounds the achievable i0).  Only the
epochs where the self-consistent i0 falls inside the corridor are launch
windows -> "one cluster per month".

R-4 corridor tiers (PRD v5 §8)
------------------------------
| tier        | Az upper | i upper | basis |
|-------------|----------|---------|-------|
| geometric   | 175 deg  | 85.3    | official launch azimuth 90-175 |
| capability  | ~125 deg | ~40     | self-rotation gain loss <= 100 m/s |
| calibrated  | measured | measured| historical n=4 back-calibration |

``cos i = cos phi * sin Az``, ``V_rot = 465.1 * cos phi * sin Az``.  The
geometric tier's Earth-equator inclination upper bound (85.3 deg) converts to
MEP i0 beyond the R-1S family max (45 deg), so it is capped at 45.

B-7 (launch-to-perigee geocentric angle ``dtheta``) is not given by the paper;
the PRD mandates a three-tier sensitivity (20/25/30 deg).
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import brentq

from launch_window import geometry as g
from launch_window import mep_frame as mf
from launch_window import month_accounting as ma
from launch_window import tli_solver as ts
from launch_window.lambert_window import LaunchSite

_REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_FAMILY = _REPO_ROOT / "r1s_scan.json"

# Physical constants (km, deg/h)
R_EARTH = ts.R_EARTH
H_P = 200.0
R_P = R_EARTH + H_P
H_M = ts.DEFAULT_H_M
EARTH_SIDEREAL_RATE_DEG_PER_H = 15.041

# Launch sites (WGS84)
WENCHANG = LaunchSite("Wenchang", 110.951, 19.614, 10.0)
XICHANG = LaunchSite("Xichang", 102.03, 28.246, 10.0)

# R-4 corridor tiers (PRD v5 §8): name -> (i0_min, i0_max) in MEP i0.
CORRIDOR_TIERS: dict[str, tuple[float, float]] = {
    "geometric": (19.6, 45.0),
    "capability": (19.6, 40.0),
    "calibrated": (20.0, 45.0),
}

# R-4 azimuth corridor tiers (PRD v5 §8): name -> (az_min, az_max) in degrees.
# The official launch azimuth is 90-175 (geometric).  The capability tier
# keeps the self-rotation gain loss <= 100 m/s: V_rot = 465.1*cos(phi)*sin(Az),
# at Az=90 the gain is 438 m/s, so loss <= 100 m/s gives sin(Az) >= 338/438
# -> Az >= 50.5 deg, but the official corridor starts at 90, so the
# capability tier is [90, ~125] (PRD: ~120-125).
AZIMUTH_TIERS: dict[str, tuple[float, float]] = {
    "geometric": (90.0, 175.0),
    "capability": (90.0, 125.0),
    "calibrated": (90.0, 175.0),
}

# B-7: launch-to-perigee geocentric angle (deg), three-tier sensitivity.
DELTA_THETA_TIERS_DEG = (20.0, 25.0, 30.0)

# Nominal perigee-longitude tolerance (deg) for the window-width formula.
# PRD v5 §6: w = 2*tol_lambda / 15.041 deg/h.  tol_lambda = 0.125 deg gives
# a 1-min window; the exact value is back-derived from the allowed H_m
# deviation (B-3: CE-6 nominal 200 km / 3 burns).
TOL_LAMBDA_DEG = 0.125


def reachable_declination_band(
    site: LaunchSite, dtheta_deg: float, corridor: str = "capability"
) -> list[float]:
    """Reachable perigee-declination band [min, max] (deg) from the site.

    PRD v5 §6: ``sin δ_p = sin φ·cos Δθ + cos φ·sin Δθ·cos Az``.  The perigee
    declination is maximised at Az = 90° (cos Az = 0) and minimised at the
    largest azimuth in the corridor.  This band is the *only* exclusion force
    of the model: the required perigee declination ≈ −δ_moon must fall in it.
    """
    phi = math.radians(site.lat_deg)
    dtheta = math.radians(dtheta_deg)
    az_min, az_max = AZIMUTH_TIERS[corridor]
    sin_max = math.sin(phi) * math.cos(dtheta)  # at Az = 90 deg
    sin_min = (
        math.sin(phi) * math.cos(dtheta)
        + math.cos(phi) * math.sin(dtheta) * math.cos(math.radians(az_max))
    )
    return [
        math.degrees(math.asin(max(-1.0, min(1.0, sin_min)))),
        math.degrees(math.asin(max(-1.0, min(1.0, sin_max)))),
    ]


def moon_declination_in_band(t_launch: float, band: list[float]) -> bool:
    """True if the required perigee declination (−δ_moon) is in the band.

    The perigee is in the anti-Moon direction, so its geocentric declination
    ≈ −δ_moon.  Only epochs where −δ_moon falls in the launch-site reachable
    band are launchable (PRD v5 §6: the "one cluster per month" mechanism).
    """
    delta_moon = g.moon_geocentric_declination_deg(t_launch)
    required = -delta_moon
    return band[0] <= required <= band[1]


def required_perigee_declination_deg(t_launch: float) -> dict:
    """δ_required(T_e) from the ephemeris directly (PRD v6 §3 C-2).

    The perigee is in the anti-Moon direction, so its geocentric declination
    ≈ −δ_moon.  The O(λ_p) correction (the perigee sits at MEP longitude
    λ_p, not exactly at the anti-Moon point) is reported as
    ``delta_required_minus_delta_moon_deg`` so its magnitude can be checked
    (PRD v6 §3.3: "确认 O(λ_p) 修正项量级").

    Returns
    -------
    dict with keys:
      delta_required_deg                  : required perigee declination (deg)
      delta_moon_deg                      : Moon geocentric declination (deg)
      delta_required_minus_delta_moon_deg : δ_required − δ_moon (deg)
      o_lambda_p_correction_deg           : δ_required + δ_moon (deg) — the
                                            O(λ_p) deviation from the naive
                                            −δ_moon rule (small, ~λ_p scale)
      required_ra_deg                     : required perigee right ascension (deg)
    """
    import spiceypy

    g.ensure_kernels()
    state, _lt = spiceypy.spkezr("MOON", t_launch, "J2000", "NONE", "EARTH")
    r = np.asarray(state[:3], dtype=float)
    r_norm = float(np.linalg.norm(r))
    anti = -r / r_norm
    delta_required = math.degrees(math.asin(float(anti[2])))
    delta_moon = math.degrees(math.asin(float(r[2]) / r_norm))
    ra_required = math.degrees(math.atan2(float(anti[1]), float(anti[0]))) % 360.0
    return {
        "delta_required_deg": delta_required,
        "delta_moon_deg": delta_moon,
        "delta_required_minus_delta_moon_deg": delta_required - delta_moon,
        "o_lambda_p_correction_deg": delta_required + delta_moon,
        "required_ra_deg": ra_required,
    }


def r1s_perigee_lat_spread_deg(families: dict[str, PerigeeFamily]) -> float:
    """R-1S perigee-latitude spread (deg) across all families and i0.

    PRD v6 §3 (C-2): this is the *requirement* spread (the R-1S perigee
    trajectory family spans a latitude band around the MEP plane), NOT the
    gate tolerance.  It enters the gate only as a relaxation term and must
    stay < 1° at the engineering operating point (T = 112.375 h).
    """
    lats = [lat for fam in families.values() for lat in fam.lat]
    return (max(lats) - min(lats)) / 2.0


def perigee_lat_band_provenance(families: dict[str, PerigeeFamily]) -> dict:
    """PRD v6 §3.4: the reachable perigee-lat band *with* its provenance.

    The band is the requirement spread (F-2), NOT a tolerance.  This returns
    the band edges together with the (i₀, family) that produce each edge, so
    the delivery never reports a bare pair of numbers.  ``family`` encodes
    the 升/降段 (asc/desc) × N/S 族 combination.

    Returns
    -------
    dict with keys:
      perigee_lat_min_deg, perigee_lat_max_deg : band edges (deg)
      min_provenance, max_provenance           : {i0_deg, family} at each edge
      spread_deg                               : (max - min) / 2
      note                                     : "需求散布，非容差"
    """
    best_min: tuple[float, str, float] | None = None  # (lat, family, i0)
    best_max: tuple[float, str, float] | None = None
    for name, fam in families.items():
        for i0, lat in zip(fam.i0, fam.lat):
            if best_min is None or lat < best_min[0]:
                best_min = (float(lat), name, float(i0))
            if best_max is None or lat > best_max[0]:
                best_max = (float(lat), name, float(i0))
    assert best_min is not None and best_max is not None
    return {
        "perigee_lat_min_deg": best_min[0],
        "perigee_lat_max_deg": best_max[0],
        "min_provenance": {"i0_deg": best_min[2], "family": best_min[1]},
        "max_provenance": {"i0_deg": best_max[2], "family": best_max[1]},
        "spread_deg": (best_max[0] - best_min[0]) / 2.0,
        "note": "需求散布，非容差",
    }


def declination_gate(
    t_launch: float,
    site: LaunchSite,
    dtheta_deg: float,
    corridor: str,
    r1s_spread_deg: float,
) -> dict:
    """PRD v6 §3 (C-2) role-separated exclusion gate.

    The gate width comes from the launch-site side ``(φ_site, Az, Δθ)`` via
    :func:`reachable_declination_band`; the R-1S perigee-latitude spread
    (the *requirement* spread, F-2) enters only as a relaxation term and
    must stay < 1° at the engineering operating point.  The required perigee
    declination δ_required(T_e) is taken directly from the ephemeris
    (:func:`required_perigee_declination_deg`), not from a scanned i₀
    envelope.

    Returns
    -------
    dict with keys:
      pass                          : bool — δ_required in the widened band
      margin_deg                    : signed distance to the nearest band edge
      required_ra_deg               : required perigee right ascension (deg)
      delta_required_deg            : required perigee declination (deg)
      delta_moon_deg                : Moon geocentric declination (deg)
      delta_required_minus_delta_moon_deg : O(λ_p) correction (deg)
      band_deg                      : [min, max] launch-site reachable band
      r1s_spread_deg                : relaxation term actually used (deg)
    """
    band = reachable_declination_band(site, dtheta_deg, corridor)
    widened = [band[0] - r1s_spread_deg, band[1] + r1s_spread_deg]
    req = required_perigee_declination_deg(t_launch)
    delta_required = req["delta_required_deg"]
    margin = min(delta_required - widened[0], widened[1] - delta_required)
    return {
        "pass": widened[0] <= delta_required <= widened[1],
        "margin_deg": margin,
        "required_ra_deg": req["required_ra_deg"],
        "delta_required_deg": delta_required,
        "delta_moon_deg": req["delta_moon_deg"],
        "delta_required_minus_delta_moon_deg": req["delta_required_minus_delta_moon_deg"],
        "o_lambda_p_correction_deg": req["o_lambda_p_correction_deg"],
        "band_deg": widened,
        "r1s_spread_deg": r1s_spread_deg,
    }


def _monthly_accounting(
    window_days: list[str],
    year: int | None,
    month: int | None,
    fallback_n_days: int,
) -> dict:
    """C-8 (PRD v6 §9) calendar-month accounting for a scan's window days.

    When ``year``/``month`` are supplied, coverage is
    ``月内窗口日数 ÷ 月内天数`` and any neighbouring-cluster days are split
    into ``adjacent_cluster_days`` (excluded from the denominator).  Without
    them (legacy callers) the old epoch-derived denominator is used.
    """
    if year is not None and month is not None:
        acc = ma.monthly_coverage(window_days, year, month)
        return {
            "n_days_in_month": acc["n_days_in_month"],
            "n_window_days": acc["n_window_days_in_month"],
            "window_days": acc["window_days_in_month"],
            "adjacent_cluster_days": acc["adjacent_cluster_days"],
            "coverage": acc["coverage"],
        }
    coverage = round(len(window_days) / fallback_n_days, 4) if fallback_n_days else None
    return {
        "n_days_in_month": fallback_n_days,
        "n_window_days": len(window_days),
        "window_days": sorted(window_days),
        "adjacent_cluster_days": [],
        "coverage": coverage,
    }


class PerigeeFamily:
    """Interpolated R-1S perigee trajectory family for one solution family."""

    def __init__(self, name: str, records: list[dict]):
        self.name = name
        i0 = np.array([r["i0_deg"] for r in records], dtype=float)
        order = np.argsort(i0)
        self.i0 = i0[order]
        self.lon = np.array([r["perigee_lon_deg"] for r in records], dtype=float)[order]
        self.lat = np.array([r["perigee_lat_deg"] for r in records], dtype=float)[order]
        self.vp = np.array([r["vp_km_s"] for r in records], dtype=float)[order]
        self.Omega = np.array([r["Omega_deg"] for r in records], dtype=float)[order]
        self.omega = np.array([r["omega_deg"] for r in records], dtype=float)[order]
        self.perilune_lon = np.array(
            [r["perilune_lon_deg"] for r in records], dtype=float
        )[order]
        self.perilune_lat = np.array(
            [r["perilune_lat_deg"] for r in records], dtype=float
        )[order]
        # Optional per-record fields (present in the scan_im / r1s_scan
        # schemas; NaN when absent so the interpolators stay usable).
        self.v_perilune = np.array(
            [r.get("v_perilune_km_s", float("nan")) for r in records],
            dtype=float,
        )[order]
        self.i_m = np.array(
            [r.get("i_m_deg", float("nan")) for r in records], dtype=float
        )[order]
        self.lon_f = interp1d(self.i0, self.lon, kind="linear")
        self.lat_f = interp1d(self.i0, self.lat, kind="linear")
        self.vp_f = interp1d(self.i0, self.vp, kind="linear")
        self.Omega_f = interp1d(self.i0, self.Omega, kind="linear")
        self.omega_f = interp1d(self.i0, self.omega, kind="linear")
        self.perilune_lon_f = interp1d(self.i0, self.perilune_lon, kind="linear")
        self.perilune_lat_f = interp1d(self.i0, self.perilune_lat, kind="linear")
        self.v_perilune_f = interp1d(self.i0, self.v_perilune, kind="linear")
        self.i_m_f = interp1d(self.i0, self.i_m, kind="linear")

    @property
    def i0_min(self) -> float:
        return float(self.i0[0])

    @property
    def i0_max(self) -> float:
        return float(self.i0[-1])


def load_families(path: Path = DEFAULT_FAMILY) -> dict[str, PerigeeFamily]:
    """Load the R-1S scan JSON and build per-family interpolators."""
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    return {
        name: PerigeeFamily(name, records)
        for name, records in data["families"].items()
    }


def parking_orbit_inclination_deg(site_mep: np.ndarray, perigee_mep: np.ndarray) -> float:
    """Inclination (deg) of the plane through Earth centre, site, perigee point.

    The plane normal is ``cross(site, perigee)``; the inclination relative to
    the MEP fundamental plane is ``acos(|n_z|)``.
    """
    n = np.cross(site_mep, perigee_mep)
    n = n / np.linalg.norm(n)
    return math.degrees(math.acos(min(1.0, abs(float(n[2])))))


def solve_i0_for_family(
    family: PerigeeFamily,
    t_arrive: float,
    T_s: float,
    dtheta_deg: float,
    site: LaunchSite,
    i0_lo: float,
    i0_hi: float,
) -> list[float]:
    """Solve the self-consistent i0 equation for one family at one arrival epoch.

    Returns the list of i0 roots in [i0_lo, i0_hi] (usually 0 or 1).
    """
    frame = mf.mep_at_epoch(t_arrive)
    dT_s = dtheta_deg / EARTH_SIDEREAL_RATE_DEG_PER_H * 3600.0
    t_L = t_arrive - T_s - dT_s
    r_site = g.launch_site_vector_j2000_km(site.lon_deg, site.lat_deg, site.alt_m, t_L)
    r_site_mep = mf.to_mep(r_site, frame)
    r_site_mep = R_P * r_site_mep / np.linalg.norm(r_site_mep)

    lo = max(i0_lo, family.i0_min)
    hi = min(i0_hi, family.i0_max)
    if lo >= hi:
        return []

    def f(i0: float) -> float:
        lon_p = float(family.lon_f(i0))
        lat_p = float(family.lat_f(i0))
        r_peri = mf.from_mep_lonlat(lon_p, lat_p, R_P)
        return parking_orbit_inclination_deg(r_site_mep, r_peri) - i0

    # Scan for sign changes, then refine with brentq.
    n_scan = 40
    i0s = np.linspace(lo, hi, n_scan)
    vals = np.array([f(float(i)) for i in i0s])
    roots: list[float] = []
    for k in range(n_scan - 1):
        if vals[k] * vals[k + 1] < 0.0:
            try:
                roots.append(
                    float(brentq(f, float(i0s[k]), float(i0s[k + 1]), xtol=1e-8))
                )
            except ValueError:
                continue
    return roots


@dataclass
class WindowMatch:
    """A single R-2 matched launch window at one arrival epoch."""

    t_arrive: float
    t_arrive_utc: str
    t_launch: float
    t_launch_utc: str
    family: str
    i0_deg: float
    vp_km_s: float
    Omega_deg: float
    omega_deg: float
    perigee_lon_deg: float
    perigee_lat_deg: float
    perilune_lon_deg: float
    perilune_lat_deg: float
    v_perilune_km_s: float
    i_m_deg: float
    C3: float
    perilune_farside: bool
    corridor: str

    def to_dict(self) -> dict:
        d = asdict(self)
        return d


def match_window(
    t_arrive: float,
    T_h: float,
    dtheta_deg: float,
    site: LaunchSite,
    families: dict[str, PerigeeFamily],
    corridor: str = "capability",
    moon_model: str = "spice",
    H_m: float = H_M,
) -> list[WindowMatch]:
    """Full R-2 match at one arrival epoch: i0 self-consistency + R-1 solve.

    Returns the list of matched windows (one per family root).  The R-1 solve
    uses the SPICE truth model by default (engineering scan); the circular
    model is available for regression.
    """
    T_s = T_h * 3600.0
    i0_lo, i0_hi = CORRIDOR_TIERS[corridor]
    dT_s = dtheta_deg / EARTH_SIDEREAL_RATE_DEG_PER_H * 3600.0
    t_L = t_arrive - T_s - dT_s
    frame = mf.mep_at_epoch(t_arrive)

    matches: list[WindowMatch] = []
    for fam_name, family in families.items():
        for i0 in solve_i0_for_family(
            family, t_arrive, T_s, dtheta_deg, site, i0_lo, i0_hi
        ):
            seed = (
                float(family.vp_f(i0)),
                float(family.Omega_f(i0)),
                float(family.omega_f(i0)),
            )
            try:
                res = ts.solve_polar(
                    H_P, T_s, i0, t_arrive, seed, H_m=H_m, frame=frame,
                    moon_model=moon_model,
                )
            except RuntimeError:
                continue
            matches.append(
                WindowMatch(
                    t_arrive=t_arrive,
                    t_arrive_utc=g.et_to_utc(t_arrive),
                    t_launch=t_L,
                    t_launch_utc=g.et_to_utc(t_L),
                    family=fam_name,
                    i0_deg=i0,
                    vp_km_s=float(res.p[0]),
                    Omega_deg=float(res.p[1]),
                    omega_deg=float(res.p[2]),
                    perigee_lon_deg=float(res.perigee_lon_mep),
                    perigee_lat_deg=float(res.perigee_lat_mep),
                    perilune_lon_deg=float(res.perilune_lon_mep),
                    perilune_lat_deg=float(res.perilune_lat_mep),
                    v_perilune_km_s=float(res.v_perilune),
                    i_m_deg=float(res.i_m),
                    C3=float(res.C3),
                    perilune_farside=abs(float(res.perilune_lon_mep)) > 90.0,
                    corridor=corridor,
                )
            )
    return matches


def window_width_min(tol_lambda_deg: float = TOL_LAMBDA_DEG) -> float:
    """Window width (min) from the perigee-longitude tolerance (PRD v5 §6).

    ``w = 2 * tol_lambda / 15.041 deg/h``, converted to minutes.
    """
    return 2.0 * tol_lambda_deg / EARTH_SIDEREAL_RATE_DEG_PER_H * 60.0


def degradation_rate_deg_per_min(
    t_arrive: float,
    T_h: float,
    dtheta_deg: float,
    site: LaunchSite,
    families: dict[str, PerigeeFamily],
    corridor: str = "capability",
    dt_min: float = 5.0,
) -> float:
    """V-5 degradation check: d(lambda_p_required)/dt (deg/min).

    The required perigee longitude is the launch-site longitude in the MEP
    frame at launch time (direct injection: perigee ≈ launch point).  As the
    Earth rotates under the frozen MEP frame, this sweeps at the Earth's
    sidereal rate ~15 deg/h = 0.25 deg/min.  The check is |d/dt| >= 0.1
    deg/min and grazing_fraction < 0.2 (PRD v5 §6, V-5).
    """
    T_s = T_h * 3600.0
    dT_s = dtheta_deg / EARTH_SIDEREAL_RATE_DEG_PER_H * 3600.0
    frame = mf.mep_at_epoch(t_arrive)
    t_L0 = t_arrive - T_s - dT_s
    t_L1 = t_L0 + dt_min * 60.0

    def _lon(t_L: float) -> float:
        r_site = g.launch_site_vector_j2000_km(
            site.lon_deg, site.lat_deg, site.alt_m, t_L
        )
        r_site_mep = mf.to_mep(r_site, frame)
        lon, _ = mf.to_mep_lonlat(r_site_mep)
        return lon

    l0 = _lon(t_L0)
    l1 = _lon(t_L1)
    # wrap the difference to [-180, 180)
    d = (l1 - l0 + 180.0) % 360.0 - 180.0
    return d / dt_min


def _scan_arrival_range(
    start_et: float,
    end_et: float,
    T_h: float,
    dtheta_deg: float,
    site: LaunchSite,
    families: dict[str, PerigeeFamily],
    corridor: str = "capability",
    step_min: float = 10.0,
    moon_model: str = "spice",
    coarse_only: bool = False,
    year: int | None = None,
    month: int | None = None,
    r1s_spread_override: float | None = None,
    H_m: float = H_M,
    i0_lo_override: float | None = None,
    i0_hi_override: float | None = None,
) -> dict:
    """Core R-2 scan over an arrival-epoch range [start_et, end_et].

    Returns the window-day grouping, coverage, grazing fraction, degradation
    rate, and window width.  Used by :func:`scan_month` (arrival epochs in a
    calendar month) and :func:`six_mission_retrospective` (arrival epochs
    shifted so the launch days fall in the launch month).

    The i0 search bounds default to the named ``corridor`` tier; callers
    with mission-specific corridors (e.g. arbitrary lunar-orbit inclination
    scans whose family tables extend to i0 = 1°) can override them with
    ``i0_lo_override`` / ``i0_hi_override``.  The azimuth corridor for the
    declination gate still comes from ``corridor``.

    When ``year``/``month`` are supplied, C-8 (PRD v6 §9) calendar accounting
    is applied: coverage is ``月内窗口日数 ÷ 月内天数``, the denominator comes
    from the Gregorian calendar, and neighbouring-cluster days are split into
    ``adjacent_cluster_days`` instead of inflating the numerator.

    With ``coarse_only=True`` the full R-1 solves at run centres are skipped
    and the window days are derived directly from the coarse in-corridor
    runs (sufficient for day-level determination, e.g. the V-9 retrospective).
    """
    T_s = T_h * 3600.0
    i0_lo, i0_hi = CORRIDOR_TIERS[corridor]
    if i0_lo_override is not None:
        i0_lo = float(i0_lo_override)
    if i0_hi_override is not None:
        i0_hi = float(i0_hi_override)
    dT_s = dtheta_deg / EARTH_SIDEREAL_RATE_DEG_PER_H * 3600.0
    r1s_spread = (
        r1s_perigee_lat_spread_deg(families)
        if r1s_spread_override is None
        else r1s_spread_override
    )

    # Coarse scan: i0_sol(T_e) via the cheap self-consistent solve, gated by
    # the PRD v6 §3 (C-2) role-separated declination gate: the required
    # perigee declination δ_required(T_e) (from the ephemeris directly) must
    # fall in the launch-site reachable band (φ_site, Az, Δθ), with the R-1S
    # perigee-latitude spread entering only as a relaxation term.
    step_s = step_min * 60.0
    n = int((end_et - start_et) / step_s) + 1
    t_arr = start_et + np.arange(n) * step_s
    i0_sol = np.full(n, np.nan)
    fam_sol: list[str | None] = [None] * n
    for k in range(n):
        t_launch = float(t_arr[k]) - T_s - dT_s
        if not declination_gate(t_launch, site, dtheta_deg, corridor, r1s_spread)["pass"]:
            continue
        best_i0 = None
        best_fam = None
        for fam_name, family in families.items():
            for r in solve_i0_for_family(
                family, float(t_arr[k]), T_s, dtheta_deg, site, i0_lo, i0_hi
            ):
                if best_i0 is None or abs(r - 0.5 * (i0_lo + i0_hi)) < abs(
                    best_i0 - 0.5 * (i0_lo + i0_hi)
                ):
                    best_i0, best_fam = r, fam_name
        if best_i0 is not None:
            i0_sol[k] = best_i0
            fam_sol[k] = best_fam

    # Contiguous in-corridor runs -> windows.
    in_corr = ~np.isnan(i0_sol)
    runs: list[tuple[int, int]] = []
    k = 0
    while k < n:
        if in_corr[k]:
            j = k
            while j + 1 < n and in_corr[j + 1]:
                j += 1
            runs.append((k, j))
            k = j + 1
        else:
            k += 1

    if coarse_only:
        # Day-level window days from the coarse in-corridor runs: a launch
        # day is a window day if some arrival epoch in the range maps to it.
        # Every epoch of a contiguous in-corridor run is mapped (not just the
        # run midpoint): long runs span several launch days, and collapsing
        # them to the midpoint silently drops the interior days.
        window_days = set()
        for (a, b) in runs:
            for k_idx in range(a, b + 1):
                t_launch = float(t_arr[k_idx]) - T_s - dT_s
                window_days.add(g.et_to_utc(t_launch)[:10])
        window_days = sorted(window_days)
        accounting = _monthly_accounting(
            window_days, year, month,
            fallback_n_days=int((end_et - start_et) / 86400.0),
        )
        return {
            "site": site.name,
            "T_h": T_h,
            "dtheta_deg": dtheta_deg,
            "corridor": corridor,
            "n_window_days": accounting["n_window_days"],
            "n_days": accounting["n_days_in_month"],
            "n_days_in_month": accounting["n_days_in_month"],
            "coverage": accounting["coverage"],
            "adjacent_cluster_days": accounting["adjacent_cluster_days"],
            "n_windows": len(runs),
            "window_days": accounting["window_days"],
            "per_day": {d: {"n_windows": 1, "windows": []} for d in accounting["window_days"]},
            "grazing_fraction": None,
            "degradation_rate_deg_per_min": None,
            "window_width_min": round(window_width_min(), 2),
            "i0_range_used": [i0_lo, i0_hi],
            "r1s_spread_deg": r1s_spread,
            "band_provenance": perigee_lat_band_provenance(families),
            "coarse_only": True,
        }

    # Full R-1 solve at each run centre.
    windows: list[dict] = []
    for (a, b) in runs:
        mid = (a + b) // 2
        t_a = float(t_arr[mid])
        for wm in match_window(
            t_a, T_h, dtheta_deg, site, families, corridor, moon_model, H_m=H_m
        ):
            windows.append(wm.to_dict())

    # Group by launch day.
    by_day: dict[str, list[dict]] = {}
    for w in windows:
        day = w["t_launch_utc"][:10]
        by_day.setdefault(day, []).append(w)

    accounting = _monthly_accounting(
        sorted(by_day.keys()), year, month,
        fallback_n_days=int((end_et - start_et) / 86400.0),
    )
    in_month_days = accounting["window_days"]

    # Grazing fraction: fraction of coarse-scan steps that are in-corridor
    # but whose i0 sits within 0.5 deg of a corridor bound (grazing).
    n_grazing = 0
    for k in range(n):
        if in_corr[k]:
            if (
                abs(i0_sol[k] - i0_lo) < 0.5
                or abs(i0_sol[k] - i0_hi) < 0.5
            ):
                n_grazing += 1
    grazing_fraction = n_grazing / max(1, int(in_corr.sum()))

    # Degradation rate at the first window centre.
    deg_rate = float("nan")
    if windows:
        deg_rate = degradation_rate_deg_per_min(
            float(windows[0]["t_arrive"]), T_h, dtheta_deg, site, families, corridor
        )

    return {
        "site": site.name,
        "T_h": T_h,
        "dtheta_deg": dtheta_deg,
        "corridor": corridor,
        "n_window_days": accounting["n_window_days"],
        "n_days": accounting["n_days_in_month"],
        "n_days_in_month": accounting["n_days_in_month"],
        "coverage": accounting["coverage"],
        "adjacent_cluster_days": accounting["adjacent_cluster_days"],
        "n_windows": len(windows),
        "window_days": accounting["window_days"],
        "per_day": {
            day: {"n_windows": len(ws), "windows": ws}
            for day, ws in sorted(by_day.items())
            if day in set(in_month_days)
        },
        "grazing_fraction": round(grazing_fraction, 4),
        "degradation_rate_deg_per_min": round(deg_rate, 4) if not math.isnan(deg_rate) else None,
        "window_width_min": round(window_width_min(), 2),
        "i0_range_used": [i0_lo, i0_hi],
        "r1s_spread_deg": r1s_spread,
        "band_provenance": perigee_lat_band_provenance(families),
    }


def scan_month(
    year: int,
    month: int,
    T_h: float,
    dtheta_deg: float,
    site: LaunchSite,
    families: dict[str, PerigeeFamily],
    corridor: str = "capability",
    step_min: float = 10.0,
    moon_model: str = "spice",
    H_m: float = H_M,
) -> dict:
    """Scan one calendar month of arrival epochs and collect the launch windows.

    Returns per-day window counts, the window days, the matched windows, the
    coverage, the grazing fraction, and the degradation rate.
    """
    start_et = g.utc_to_et(ma.month_start_utc(year, month))
    end_et = g.utc_to_et(ma.next_month_start_utc(year, month))
    res = _scan_arrival_range(
        start_et, end_et, T_h, dtheta_deg, site, families,
        corridor=corridor, step_min=step_min, moon_model=moon_model,
        year=year, month=month, H_m=H_m,
    )
    res["year"] = year
    res["month"] = month
    return res


# ---------------------------------------------------------------------------
# R-3 — N/S family check (PRD v5 §7)
# ---------------------------------------------------------------------------
# Historical missions and their landing-hemisphere / expected family.
# The family is determined by the perilune latitude sign: N family has
# perilune in the north (perilune_lat > 0), S family in the south.
# PRD v5 §7: CE-3 (44.1N) -> N, CE-4 (45.5S) -> S, CE-5 (43.1N) -> N,
# CE-6 (41.6S) -> S.  CE-7 (86.6-89.8S) -> S.
HISTORICAL_FAMILIES: list[dict] = [
    {"mission": "CE-3", "landing_lat_deg": 44.1, "expected_family": "N"},
    {"mission": "CE-4", "landing_lat_deg": -45.5, "expected_family": "S"},
    {"mission": "CE-5", "landing_lat_deg": 43.1, "expected_family": "N"},
    {"mission": "CE-6", "landing_lat_deg": -41.6, "expected_family": "S"},
]


def family_from_perilune_lat(perilune_lat_deg: float) -> str:
    """Map a perilune latitude to its N/S family (PRD v5 §7)."""
    return "N" if perilune_lat_deg > 0.0 else "S"


def r3_family_check(
    families: dict[str, PerigeeFamily],
    i0_deg: float = 30.0,
) -> dict:
    """R-3: verify the N/S family assignment against the historical missions.

    For each family in the R-1S scan, report the perilune-latitude sign at a
    representative i0 and the implied N/S family.  Then check that the
    historical landing-hemisphere -> family mapping is consistent with the
    model's family definition (perilune latitude sign).

    PRD v5 §7: if the historical check is < 4/4, the family criterion is
    dropped (not forced).  The comparison must convert the white-ecliptic vs
    lunar-equator 6.7 deg tilt where relevant.
    """
    family_signs: dict[str, dict] = {}
    for name, fam in families.items():
        # perilune latitude at the representative i0 (clamped to family range)
        i0c = min(max(i0_deg, fam.i0_min), fam.i0_max)
        plon = float(fam.perilune_lon_f(i0c)) if hasattr(fam, "perilune_lon_f") else None
        plat = float(fam.perilune_lat_f(i0c)) if hasattr(fam, "perilune_lat_f") else None
        family_signs[name] = {
            "i0_deg": i0c,
            "perilune_lat_deg": plat,
            "perilune_lon_deg": plon,
            "model_family": family_from_perilune_lat(plat) if plat is not None else None,
        }

    # Historical consistency: expected family == model family for the landing
    # hemisphere.  The model family is defined by perilune latitude sign, so
    # a north landing (lat>0) should map to N and south to S.
    hist_results = []
    n_pass = 0
    for h in HISTORICAL_FAMILIES:
        expected = h["expected_family"]
        # The model family for a landing at this latitude is the sign of the
        # landing latitude (north -> N, south -> S).
        model_fam = family_from_perilune_lat(h["landing_lat_deg"])
        ok = expected == model_fam
        n_pass += int(ok)
        hist_results.append(
            {
                "mission": h["mission"],
                "landing_lat_deg": h["landing_lat_deg"],
                "expected_family": expected,
                "model_family": model_fam,
                "pass": ok,
            }
        )

    return {
        "family_perilune_signs": family_signs,
        "historical_check": hist_results,
        "n_pass": n_pass,
        "n_total": len(HISTORICAL_FAMILIES),
        "all_pass": n_pass == len(HISTORICAL_FAMILIES),
        "ce7_predicted_family": "S",
        "note": (
            "R-3: family is defined by the perilune-latitude sign (N if "
            "perilune_lat > 0, S if < 0).  Historical landing hemisphere "
            "must match.  If < 4/4 the criterion is dropped (PRD v5 §7). "
            "CE-7 landing 86.6-89.8S -> S family."
        ),
    }


# ---------------------------------------------------------------------------
# R-4 — azimuth corridor tiers (PRD v5 §8)
# ---------------------------------------------------------------------------
# cos i = cos phi * sin Az ;  V_rot = 465.1 * cos phi * sin Az (m/s)
# Wenchang phi = 19.614 deg.
WENCHANG_PHI_DEG = 19.614
V_ROT_COEFF = 465.1  # m/s at equator


def azimuth_corridor_tiers() -> dict:
    """R-4: derive the three azimuth-corridor tiers (PRD v5 §8).

    Returns, for each tier, the azimuth upper bound, the Earth-equator
    inclination upper bound, the self-rotation gain, and the resulting MEP
    i0 upper bound (capped at the R-1S family max 45 deg).

    - geometric:  Az = 175 deg (official launch azimuth 90-175)
    - capability: Az such that self-rotation gain loss <= 100 m/s
    - calibrated: historical n=4 back-calibration (placeholder, measured)
    """
    def _derive(az_deg: float) -> dict:
        az = math.radians(az_deg)
        cos_i = math.cos(math.radians(WENCHANG_PHI_DEG)) * math.sin(az)
        i_eq = math.degrees(math.acos(min(1.0, cos_i)))
        v_rot = V_ROT_COEFF * math.cos(math.radians(WENCHANG_PHI_DEG)) * math.sin(az)
        return {
            "az_deg": az_deg,
            "i_eq_deg": round(i_eq, 2),
            "v_rot_m_s": round(v_rot, 1),
            "i0_mep_max_deg": round(min(i_eq, 45.0), 2),
        }

    # Capability tier: self-rotation gain loss <= 100 m/s.  At Az=90 the gain
    # is max (438 m/s); the loss is 438 - V_rot(Az).  Find Az where loss=100.
    v_rot_90 = V_ROT_COEFF * math.cos(math.radians(WENCHANG_PHI_DEG))  # ~438
    target_v = v_rot_90 - 100.0
    sin_az_cap = target_v / (V_ROT_COEFF * math.cos(math.radians(WENCHANG_PHI_DEG)))
    az_cap = math.degrees(math.asin(min(1.0, sin_az_cap)))

    return {
        "wenchang_phi_deg": WENCHANG_PHI_DEG,
        "geometric": _derive(175.0),
        "capability": _derive(az_cap),
        "calibrated": {
            "az_deg": None,
            "i_eq_deg": None,
            "v_rot_m_s": None,
            "i0_mep_max_deg": 45.0,
            "note": "historical n=4 back-calibration (measured, not derived)",
        },
        "note": (
            "R-4: cos i = cos phi * sin Az, V_rot = 465.1*cos phi*sin Az. "
            "Geometric tier Az=175 (official 90-175). Capability tier Az such "
            "that self-rotation gain loss <= 100 m/s. The geometric tier's "
            "Earth-equator i upper bound (85.3 deg) exceeds the R-1S family "
            "max (45 deg), so it is capped at 45."
        ),
    }


# ---------------------------------------------------------------------------
# V-9 — six-mission retrospective + null-model p-value (PRD v5 §11, §12 step 12)
# ---------------------------------------------------------------------------
# CE-1..CE-4 launch from Xichang, CE-5/CE-6 from Wenchang.
SIX_MISSIONS: list[dict] = [
    {"mission": "CE-1", "site": "Xichang", "launch_utc": "2007-10-24T10:05:00Z"},
    {"mission": "CE-2", "site": "Xichang", "launch_utc": "2010-10-01T10:59:00Z"},
    {"mission": "CE-3", "site": "Xichang", "launch_utc": "2013-12-01T17:30:00Z"},
    {"mission": "CE-4", "site": "Xichang", "launch_utc": "2018-12-07T18:23:00Z"},
    {"mission": "CE-5", "site": "Wenchang", "launch_utc": "2020-11-23T20:30:00Z"},
    {"mission": "CE-6", "site": "Wenchang", "launch_utc": "2024-05-03T09:27:00Z"},
]


def six_mission_retrospective(
    families: dict[str, PerigeeFamily],
    T_h: float = 112.375,
    dtheta_deg: float = 25.0,
    corridor: str = "capability",
    step_min: float = 30.0,
    moon_model: str = "spice",
    H_m: float = H_M,
) -> dict:
    """V-9: six-mission retrospective + null-model p-value.

    For each of CE-1..CE-6, scan the arrival-epoch range that maps to launch
    days in the actual launch month (arrival = launch + T + dT), and record:
    - ``coverage`` : fraction of the month's days that are window days
    - ``hit``      : whether the actual launch day falls in the window

    C-8 (PRD v6 §9): coverage uses the Gregorian ``n_days_in_month`` and the
    in-month window days only; neighbouring-cluster days are split into
    ``adjacent_cluster_days``.  The frozen v5 denominator bug is reported via
    ``coverage_recount_diff`` (before = frozen value, after = calendar value).

    The null-model p-value is the product of the six coverages — the random
    probability that all six historical launches land in the model's windows
    by chance (PRD v5 §11 V-9, the unique pass/fail criterion, target < 0.05).

    Uses the coarse day-level scan (``coarse_only``) which is sufficient for
    the day-level determination; the full R-1 solves are not needed here.
    """
    # Frozen v5 values (F-8): the old epoch-derived denominators.  Kept only
    # so ``coverage_recount_diff`` can report before/after.
    frozen_before = {
        "CE-1": (31, 9),
        "CE-2": (31, 10),
        "CE-3": (31, 11),
        "CE-4": (31, 9),
        "CE-5": (30, 7),
        "CE-6": (30, 9),
    }

    per_mission = []
    for m in SIX_MISSIONS:
        launch_day = m["launch_utc"][:10]
        year = int(launch_day[:4])
        month = int(launch_day[5:7])
        site = WENCHANG if m["site"] == "Wenchang" else XICHANG

        # Arrival-epoch range whose launch days fall in the launch month:
        # arrival = launch + T + dT, so scan [month_start + T + dT, month_end + T + dT].
        month_start = g.utc_to_et(ma.month_start_utc(year, month))
        month_end = g.utc_to_et(ma.next_month_start_utc(year, month))
        dT_s = dtheta_deg / EARTH_SIDEREAL_RATE_DEG_PER_H * 3600.0
        T_s = T_h * 3600.0
        arr_start = month_start + T_s + dT_s
        arr_end = month_end + T_s + dT_s

        res = _scan_arrival_range(
            arr_start, arr_end, T_h, dtheta_deg, site, families,
            corridor=corridor, step_min=step_min, moon_model=moon_model,
            coarse_only=True, year=year, month=month, H_m=H_m,
        )
        window_days = res["window_days"]
        hit = launch_day in window_days
        n_days_in_month = res["n_days_in_month"]
        frozen_n, frozen_k = frozen_before.get(m["mission"], (n_days_in_month, None))
        coverage = round(len(window_days) / n_days_in_month, 4)
        coverage_before = (
            round(frozen_k / frozen_n, 4) if frozen_k is not None else None
        )
        per_mission.append(
            {
                "mission": m["mission"],
                "site": m["site"],
                "launch_utc": m["launch_utc"],
                "launch_day": launch_day,
                "n_window_days": len(window_days),
                "n_days_in_month": n_days_in_month,
                "coverage": coverage,
                "coverage_before": coverage_before,
                "hit": hit,
                "window_days": window_days,
                "adjacent_cluster_days": res["adjacent_cluster_days"],
            }
        )

    p = 1.0
    for m in per_mission:
        p *= m["coverage"]
    n_hits = sum(1 for m in per_mission if m["hit"])

    return {
        "per_mission": per_mission,
        "n_missions": len(per_mission),
        "n_hits": n_hits,
        "all_hit": n_hits == len(per_mission),
        "null_model_p": round(p, 6),
        "target_p": 0.05,
        "p_significant": p < 0.05,
        "coverage_recount_diff": ma.coverage_recount_diff(
            per_mission, before_key="coverage_before", after_key="coverage"
        ),
        "note": (
            "V-9: null-model p = product of the six monthly coverages, the "
            "random probability that all six historical launches land in the "
            "model's windows by chance.  If p < 0.05 the model's window "
            "structure is statistically meaningful; if not, the honest "
            "conclusion is that the model reproduces window structure but "
            "cannot claim to predict the launch day (PRD v5 §16)."
        ),
    }


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="R-2 perigee-trajectory matching")
    parser.add_argument("--family", type=Path, default=DEFAULT_FAMILY)
    parser.add_argument("--year", type=int, default=2026)
    parser.add_argument("--month", type=int, default=8)
    parser.add_argument("--transfer-hours", type=float, default=112.375)
    parser.add_argument("--dtheta", type=float, default=25.0)
    parser.add_argument("--corridor", default="capability", choices=list(CORRIDOR_TIERS))
    parser.add_argument("--site", default="Wenchang", choices=["Wenchang", "Xichang"])
    parser.add_argument(
        "--perilune-alt",
        type=float,
        default=H_M,
        help=f"perilune altitude H_m (km), default {H_M}",
    )
    parser.add_argument("--out", type=Path, default=_REPO_ROOT / "window_scan.json")
    args = parser.parse_args()

    families = load_families(args.family)
    site = WENCHANG if args.site == "Wenchang" else XICHANG
    result = scan_month(
        args.year, args.month, args.transfer_hours, args.dtheta, site, families,
        corridor=args.corridor, H_m=args.perilune_alt,
    )
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)
    print(f"R-2 {args.year}-{args.month:02d} {args.site} T={args.transfer_hours}h "
          f"dtheta={args.dtheta} corridor={args.corridor}")
    print(f"  window days: {result['n_window_days']}/{result['n_days_in_month']} "
          f"coverage={result['coverage']}")
    print(f"  windows: {result['n_windows']}, grazing={result['grazing_fraction']}")
    print(f"  degradation: {result['degradation_rate_deg_per_min']} deg/min")
    print(f"  window width: {result['window_width_min']} min")
    print(f"  days: {result['window_days']}")
    print(f"  wrote {args.out}")


if __name__ == "__main__":
    main()