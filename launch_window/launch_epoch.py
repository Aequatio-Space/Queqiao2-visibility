"""launch_window/launch_epoch.py — L6 launch-epoch outer solver (PRD v7 §五).

Implements PRD v7 §五: given a target arrival epoch ``T_e``, find the launch
epoch ``t_L`` such that the parking-orbit plane (defined by the launch site's
geodetic coordinates and the orbit plane carrying the launch point and the
TLI/perigee point) is consistent with the Earth-fixed launch-site longitude
via GMST, with the parking-orbit coast placing TLI at the right argument of
latitude, and with ``t_TLI + T = T_e``.

The outer solver is a single-variable Newton step on

    G(T_e) = [t_TLI(T_e) + T_transfer] − T_e

using the analytical derivative hint from PRD v7 §五,
``dG/dT_e ≈ 1 − ω_moon/ω_earth ≈ 0.967``.  Because the relationship between
``T_e`` and ``t_TLI`` is mediated by the L5 perigee-family table (which gives
``i0`` for a chosen parking-orbit inclination) and a small inner fixed-point
between the parking-orbit inclination ``i_E`` and the family parameter
``i_0``, the inner link is built by reusing
``perigee_match.solve_i0_for_family`` (the same idea that ``match_window``
uses, R-2 §6).  Only the *lookup* is reused; no existing logic is rewritten.

Geometry chain (per PRD v7 §五)
------------------------------
1. ``T_e → M(Te)``  freeze the MEP frame at arrival.
2. Inner fixed-point ``i0 ↔ i_E`` (paper §6 / §5.2):
   * Given ``i0`` → ``(Ω_E, ω_E, v_p) = mep_orbit_to_equator(...)`` via the
     family table and the frozen rotation.
   * ``i_E`` from that plane; the launch-site constraint ``|φ_L| ≤ i_E ≤ 180°−|φ_L|``
     feeds back into the next ``i0`` guess.
3. Launch instant: ``GMST(t_L) = α_L − λ_L`` where ``α_L = Ω_E + Δλ`` with
   ``tan Δλ = tan u_L · cos i_E`` and ``sin u_L = sin φ_L / sin i_E``.
   Two branches (ascending / descending, ``u_L`` and ``180°−u_L``).
4. Coast: the TLI point is the parking-orbit perigee (true anomaly 0) and
   the insertion point sits at argument of latitude
   ``u_ins = u_L + ASCENT_ARC_DEG`` (the ascent geocentric arc; ``u_L`` is
   the launch-site argument of latitude).  On the
   elliptical CZ-5 parking orbit (200 × 380 km) the coast is the elapsed
   mean-anomaly time from insertion to the next perigee, computed with
   Kepler's equation (see :func:`coast_from_argument`).  Short and zero
   coasts are valid (``COAST_MIN_S = 0``, F0c: no whole-revolution
   extension / no long-coast branch); the upper bound ``COAST_MAX_S`` is
   the CZ-5 single-orbit cap.  The reachable-inclination envelope is
   evaluated at the insertion point (F0b), not the launch-site geodetic
   latitude.
5. ``t_TLI = t_L + T_ASCENT_S + coast``;
   residual ``G(T_e) = t_TLI + T − T_e``.

Outputs
-------
``Window`` dataclass with ``to_dict()`` returning the PRD v7 §一 schema
(excluding mission-derived keys: ``t_open``, ``t_close``, ``dv``, ``margin_kg``
which depend on the launcher capability ladder and are filled by
``launcher.py``).

Conventions
-----------
- All epochs are SPICE ET (seconds past J2000).  UTC strings are produced via
  ``spiceypy.et2utc``.
- GMST is computed via ``spiceypy.et2lst`` for the J2000-equatorial frame at
  the launch site's geodetic longitude; the ``SPICE`` angle is converted to
  degrees in [0, 360).  This matches SPICE to < 1 s equivalent (≈ 15″),
  well below the 1-s tolerance used in the PRD v7 §五 outer convergence.
- The "polar_ref" switch (PRD v7 §四.1) is honoured only by the inner
  ``solve_polar`` that re-solves the transfer orbit when computing the
  reported ``i_m_deg``; the L6 outer geometry uses the L5 table directly.
  When ``polar_ref='moon_equator'`` the L5 table is still the same
  (the table is built at scan time in the ``mep`` frame; the switching
  happens when validating the resulting orbit).  For PRD v7 §五 the
  ``polar_ref`` argument is propagated to a downstream solver hook for
  callers that want to validate / re-solve at the returned ``T_e``.

References
----------
- PRD v7 §五 (L6 solver algorithm)
- PRD v7 §四.1 (q3 polar-ref switch)
- Paper §5 (R-1S family) / §6 (R-2 launch-epoch iteration)
"""
from __future__ import annotations

import logging
import math
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Any, Iterable, Literal, Optional, Sequence

import numpy as np

from launch_window import geometry as g
from launch_window import mep_frame as mf
from launch_window import perigee_match as pm
from launch_window import tli_solver as ts
from launch_window import lambert as lambert_mod
from launch_window.lambert_window import LaunchSite

# Re-use the L5 perigee-family loader.  ``pm.DEFAULT_FAMILY`` defaults to
# ``r1s_scan.json`` (T=72 h, circular).  For PRD v7 §十 the engineering scan
# uses T=112.375 h, so callers should pass an explicit path to the right
# family table.

# ---------------------------------------------------------------------------
# Physical & methodological constants (whitelisted, not mission-derived)
# ---------------------------------------------------------------------------
# Ascent duration from launch to parking-orbit insertion (s).  PRD v7 §五
# names it ``T_ASCENT_S``.  The legacy module default (1800 s — matching
# ``lambert_window._PARKING_S``) is retained for generic callers; the
# CE-5/CE-6 engineering values are the measured launch→park times in
# ``T_TO_PARK_S`` (CZ-5 second-stage first cutoff).
T_ASCENT_S = 1800.0

# Measured launch→parking-orbit-insertion time (s) for the CZ-5 lunar
# missions: CE-5 second-stage first cutoff at T+740 s, CE-6 at T+720 s.
# Callers modelling the real CZ-5 ascent should pass this instead of the
# legacy generic 1800 s placeholder.
T_TO_PARK_S = 740.0

# Parking-orbit coast bounds (s).  F0c policy (2026-09): short and zero
# coasts are allowed — the coast is the free phase variable and a root may
# legitimately sit below the CE-5/CE-6 measured ~937/1020 s single-orbit
# value.  ``COAST_MIN_S = 0.0`` therefore disables the whole-revolution
# extension (no long-coast branch is generated); the upper bound remains the
# CZ-5 single-orbit cap.
COAST_MIN_S = 0.0
COAST_MAX_S = 1000.0

# Nominal parking-orbit coast (s) used to seed outer-solver arrival guesses.
COAST_NOMINAL_S = 1000.0

# Parking-orbit geometry.  The CZ-5/CE-5/CE-6 parking orbit is a 200 × 380 km
# ellipse rather than a circle: the second-stage first cutoff inserts the
# stack at argument of latitude ``u_L`` on that ellipse and the TLI burn is
# performed at perigee (``H_p = 200 km``).  ``H_PARK_KM`` is kept as the
# legacy perigee-altitude alias used by the C3 / Δv formulas.
H_PARK_PERI_KM = 200.0
H_PARK_APO_KM = 380.0
H_PARK_KM = H_PARK_PERI_KM

# Geocentric arc from the launch site to the parking-orbit insertion point
# (deg).  The coast to the TLI perigee starts at the *insertion* argument of
# latitude ``u_ins = u_L + ASCENT_ARC_DEG``, not at the launch-site argument
# ``u_L`` (the legacy convention implicitly treated the site as the insertion
# point).  Default from the physical ascent model
# (``perigee_reachability.py``: 740 s power-law ascent on the CZ-5 200x380 km
# parking orbit -> insertion ~255 km altitude, ~20.84 deg downrange).
# Set 0.0 to reproduce the legacy site-as-insertion convention.
ASCENT_ARC_DEG = 20.8377

# Earth sidereal rotation rate (deg/hour).  PRD v7 §五 uses 15.041 deg/h
# (matches the rate used in ``perigee_match.EARTH_SIDEREAL_RATE_DEG_PER_H``).
EARTH_SIDEREAL_RATE_DEG_PER_H = 15.041

# Newton derivative hint (PRD v7 §五: ``dG/dT_e ≈ 1 − ω_moon/ω_earth ≈ 0.967``).
# Used as the step factor when the outer Newton iteration falls back to its
# analytical linearisation.  Numerical derivatives are preferred when
# available (the inner fixed-point makes the function cheap to evaluate).
DG_DTE_HINT = 0.967

# Default convergence tolerance for the outer Newton iteration (seconds).
# PRD v7 §五: ``if abs(G) < 1.0: break   # 1 second``.
G_TOL_S = 1.0

# REQ-LW-2026-LAMBERT-IG-001 §八.3: master switch for the Lambert initial
# guess inside ``_G_te``'s self-consistency loop.  When ``True`` the
# Lambert-derived (i_park, raan, omega) replaces the family-table
# interpolation; failures fall back to the family-table path.  Setting
# ``False`` restores the baseline bit-identical numerical behaviour
# (REQ §九.1.4).
USE_LAMBERT_GUESS: bool = True

# Module logger for fallback warnings.
_log = logging.getLogger(__name__)

# Sidereal day length (seconds).  Used by ``scan_epochs`` to seed the
# Newton iteration on a grid of arrival epochs one sidereal day apart.
SIDEREAL_DAY_S = 86164.0905

# Maximum outer Newton iterations.
OUTER_MAX_IT = 12

# Maximum inner fixed-point iterations (i0 ↔ i_E).
INNER_FIX_MAX_IT = 12

# Dedup threshold for ``scan_epochs`` (seconds).  Two ``t_launch_et`` values
# within this threshold are considered the same launch opportunity and the
# second is dropped.
DEDUP_THRESHOLD_S = 600.0

# Max allowed delta (deg) between inbound and outbound parking-orbit
# planes; above this, the inner fixed-point is considered geometrically
# infeasible.  Used by ``_solve_inner_fixed_point_with_code`` to tag
# failures with ``"dtheta_exceeded"`` so callers can distinguish them
# from envelope / sign-change / GMST failures.
_DTHETA_MAX_DEG = 180.0

# Maximum |residual| (deg) accepted for a GMST inverse root.  The wrapped
# residual ``(gmst - target + 540) % 360 - 180`` jumps between +180 and
# -180 at ``gmst = target + 180``; Brent's method also brackets that
# discontinuity and returns a pseudo-root.  True roots have |residual| ~ 0,
# so the pseudo-roots are rejected by this tolerance.
GMST_ROOT_TOL_DEG = 1.0

# Vocabulary of fail codes returned by ``_solve_inner_fixed_point_with_code``
# and propagated through ``_G_te`` / ``_pocket_search`` / ``_make_reject``.
# Each code corresponds to a distinct geometric / numerical failure mode
# of the L6 outer solver chain.
_FAIL_CODES = frozenset({
    "inc_below_site_lat",     # inc_E < |phi_L|  (plane below launch site)
    "inc_above_supplement",   # inc_E > 180 - |phi_L|  (prograde-retrograde mix)
    "i0_out_of_family_range", # no sign-change in i0 scan / lo >= hi
    "dtheta_exceeded",        # dtheta_deg > _DTHETA_MAX_DEG (pre-check)
    "lambert_fail",           # _lambert_parking_elements returned None
    "gmst_fail",              # gmst_to_et raised ValueError (no root)
    "coast_exceeds_limit",    # parking-orbit coast above the solver cap
})


PolarRef = Literal["mep", "moon_equator"]


# ---------------------------------------------------------------------------
# Window dataclass (PRD v7 §一, L6 schema subset — dv / t_open filled by L7)
# ---------------------------------------------------------------------------
@dataclass
class Window:
    """A single L6 launch-window record (PRD v7 §一).

    ``to_dict()`` returns a flat dictionary with the schema fields used by
    downstream L7/L8 stages (launcher capability, period analysis).
    """

    t_launch_utc: str
    t_tli_utc: str
    t_arrive_utc: str
    T_h: float
    family: str
    branch: str  # 'ascending' | 'descending'
    i_park_deg: float
    raan_park_deg: float
    azimuth_deg: float
    coast_s: float
    vp_kms: float
    c3_km2s2: float
    i_m_deg: float
    perilune_alt_km: float
    polar_ref: str
    t_launch_et: float
    t_arrive_et: float

    def to_dict(self) -> dict:
        return asdict(self)


# ---------------------------------------------------------------------------
# WindowReject dataclass — per-rejection diagnostic record (side-channel).
# ---------------------------------------------------------------------------
@dataclass
class WindowReject:
    """Per-rejection diagnostic record from ``scan_epochs``.

    Mirrors :class:`landing_window.sites.RejectedSite` /
    :class:`launch_window_orekit.coarse_scan.CoarseArrivalCandidate`.
    Captures the in-scope geometry / family / branch context at every
    silent drop point of the L6 solver so callers can debug which
    sub-stage caused each rejection.

    ``to_dict()`` returns ``asdict(self)`` for downstream JSON serialisation
    (no specialised schema).
    """

    te0_et: float             # seed arrival epoch (SPICE ET)
    t_h_hours: float          # transfer time in hours
    family: str               # family ID
    branch: str               # 'ascending' | 'descending'
    reason: str               # short reason tag
    t_launch_et: Optional[float] = None
    inc_e_deg: Optional[float] = None
    raan_e_deg: Optional[float] = None
    omega_e_deg: Optional[float] = None
    vp_kms: Optional[float] = None
    azimuth_deg: Optional[float] = None
    u_l_deg: Optional[float] = None
    coast_s: Optional[float] = None
    c3_km2s2: Optional[float] = None
    g_residual_s: Optional[float] = None
    i0_deg: Optional[float] = None
    extra: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return asdict(self)


# ---------------------------------------------------------------------------
# GMST helpers — SPICE-consistent (forward + inverse)
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Forward ``gmst_deg`` uses SPICE's ``pxform('IAU_EARTH','J2000', et)``: the
# first column of the rotation matrix is the unit vector pointing from the
# Earth centre to the prime meridian in the J2000 inertial frame.  The GMST
# is the equatorial angle of that vector:
#
#     gmst_rad = atan2(M[1, 0], M[0, 0])
#
# This is the canonical GMST definition used by the SPICE PCK and is the
# ONLY definition that is internally consistent with every other frame
# transform in this engine (rotation of the launch-site vector, transfer-
# plane normal, etc.).  Earlier revisions used ``spiceypy.et2lst`` which,
# as of the v6 round, was discovered to return *local solar time* (LMST),
# not GMST — the offset was ~35°-110° at our test epochs and the per-UT-day
# rate was 0.046 deg/day instead of the expected 0.986 deg/day.  See
# ``gmst_et2lst_deg`` (kept for diagnostic use only) for the broken branch.
#
# The IAU 1982 closed-form formula ``gmst_analytic_deg`` is also kept for
# cross-check, but it is not consistent with the modern PCK (the PCK uses
# the IAU 2006 GMST; the IAU 1982 formula drifts by ~25 arcmin/yr from
# the truth at modern epochs).  The test suite compares pxform vs analytic
# only as a "no-wild-discrepancy" check.
#
# The inverse ``gmst_to_et`` uses a coarse sign-change scan on
# ``gmst_deg(et) - target`` plus a single Brent refinement.  The scan
# radius is ± 1.2 sidereal days (so at least one wrap is crossed); brentq
# refines to xtol = 1e-6 sec (well below 1 s of arc).

def gmst_deg(et: float) -> float:
    """GMST at epoch ``et`` (SPICE ET, seconds past J2000), degrees in [0, 360).

    Defined via ``pxform('IAU_EARTH', 'J2000', et)``: the equatorial angle
    of the IAU_EARTH X-axis (prime-meridian direction) in the J2000
    inertial frame.  This is the canonical GMST used by the rest of the
    engine — internally consistent with the launch-site J2000 vector,
    the transfer-plane normal, the Moon position, etc.
    """
    import spiceypy

    g.ensure_kernels()
    M = spiceypy.pxform("IAU_EARTH", "J2000", float(et))
    # M[:, 0] is the unit IAU_EARTH X-axis expressed in J2000.  Its
    # equatorial-plane projection angle relative to J2000 X-axis is GMST.
    angle_rad = math.atan2(float(M[1, 0]), float(M[0, 0]))
    return math.degrees(angle_rad) % 360.0


def gmst_pxform_deg(et: float) -> float:
    """Alias of :func:`gmst_deg` — pxform-based GMST, the production path."""
    return gmst_deg(et)


def gmst_et2lst_deg(et: float) -> float:
    """GMST via ``spiceypy.et2lst`` — kept for diagnostic use only.

    WARNING: as of the v6 round, ``spiceypy.et2lst(et, 399, 0, 'PLANETOCENTRIC')``
    returns local *solar* time at the prime meridian (LMST), NOT GMST.
    The two differ by the equation of time (≈ ±4 min ≈ ±1°) plus the
    sidereal/solar rate discrepancy (0.986 vs 1.000 deg/h).  At the 2026
    epoch the discrepancy vs the pxform-based GMST is ~35°-110°, drifting
    ~0.78 deg/UT-day.  The production path MUST use :func:`gmst_deg`
    (pxform-based), not this function.
    """
    import spiceypy

    g.ensure_kernels()
    result = spiceypy.et2lst(et, 399, 0.0, "PLANETOCENTRIC")
    hr, mn, sc = result[0], result[1], result[2]
    angle_deg = (hr + mn / 60.0 + sc / 3600.0) * 15.0  # 24h -> 360°
    return angle_deg % 360.0


def gmst_analytic_deg(et: float) -> float:
    """IAU 1982 closed-form GMST (deg, [0, 360)) — for cross-check only.

    IAU 1982 GMST at JD_UT1:

        GMST_deg = 280.46061837
                  + 360.98564736629 · (JD_UT1 − 2451545.0)
                  + 0.000387933 · T²
                  − T³ / 38710000

    where T = (JD_UT1 − 2451545.0) / 36525.  At modern (post-2010) epochs
    the IAU 1982 formula drifts from the IAU 2006/PCK truth by ~25 arcmin/yr,
    so we expose it for sanity tests only — never use it as the canonical GMST
    in production.  Use :func:`gmst_deg` (pxform-based) instead.
    """
    import spiceypy
    g.ensure_kernels()
    jd_tdb = float(spiceypy.unitim(et, "ET", "JDTDB"))
    d_days = jd_tdb - 2451545.0
    T = d_days / 36525.0
    gmst = (
        280.46061837
        + 360.98564736629 * d_days
        + 0.000387933 * T * T
        - T * T * T / 38710000.0
    )
    return gmst % 360.0


# IAU 1982 GMST polynomial coefficients (used by the test suite only).
_GMST_COEF_A = 280.46061837
_GMST_COEF_B = 360.98564736629
_GMST_COEF_C = 0.000387933
_GMST_COEF_D = 1.0 / 38710000.0


def gmst_to_et(gmst_target_deg: float, et_guess: float) -> float:
    """Inverse GMST: SPICE ET such that ``gmst_deg(et) == gmst_target_deg``.

    Searches in ``et_guess ± 1.2 sidereal days`` and refines with Brent's
    method on the wrapped difference ``gmst_deg(et) − target``.  Returns
    the solution closest to ``et_guess`` (within ± 1 sidereal day; the
    caller should pick a guess consistent with the launch geometry).

    Refinement xtol is 1e-6 sec of ET (≪ 1 arc-second of GMST, since
    GMST rate ≈ 360.9856 deg/UT-day = 1.5e-5 deg/sec; 1e-6 sec ×
    1.5e-5 deg/sec ≈ 5e-12 deg, far below the 1-sec/15″ outer
    tolerance).
    """
    target = gmst_target_deg % 360.0

    def residual(et: float) -> float:
        d = (gmst_deg(et) - target + 540.0) % 360.0 - 180.0  # signed wrap
        return d

    from scipy.optimize import brentq

    # ± 1.2 sidereal days: at least one GMST wrap is crossed.
    lo = et_guess - 1.2 * SIDEREAL_DAY_S
    hi = et_guess + 1.2 * SIDEREAL_DAY_S
    # Coarse sign-change scan (24 samples = 2.4 hr spacing).
    n = 24
    xs = np.linspace(lo, hi, n)
    vals = np.array([residual(float(x)) for x in xs])
    found: list[tuple[float, float, float]] = []
    for k in range(n - 1):
        if vals[k] * vals[k + 1] < 0.0:
            try:
                root = brentq(
                    residual,
                    float(xs[k]),
                    float(xs[k + 1]),
                    xtol=1e-6,
                )
                found.append((root, abs(residual(root)), abs(root - et_guess)))
            except ValueError:
                continue
    # Drop discontinuity pseudo-roots: the sawtooth residual also changes
    # sign across ``gmst = target + 180`` and Brent converges to that jump
    # (|residual| ~ 180 deg), which is NOT a solution of gmst == target.
    # Only roots with |residual| <= GMST_ROOT_TOL_DEG are true inverses.
    true_roots = [(root, resid, dist) for (root, resid, dist) in found
                  if resid <= GMST_ROOT_TOL_DEG]
    if not true_roots:
        raise ValueError(
            f"gmst_to_et: no true root for GMST {target:.6f} deg near "
            f"et={et_guess:.3f} (found {len(found)} pseudo-root(s))"
        )
    # Pick the root closest to ``et_guess`` (not the smallest residual: a
    # perfectly-solved root one parking-orbit period away is geometrically
    # wrong because the launch epoch is a single ET).
    true_roots.sort(key=lambda x: x[2])
    return true_roots[0][0]


# ---------------------------------------------------------------------------
# Parking-orbit geometry helpers
# ---------------------------------------------------------------------------
def parking_orbit_period_s(
    h_peri_km: float = H_PARK_PERI_KM,
    h_apo_km: Optional[float] = H_PARK_APO_KM,
) -> float:
    """Keplerian period of the (possibly elliptical) parking orbit (s).

    The engineering default is the CZ-5 200 × 380 km ellipse used by the
    CE-5/CE-6 missions.  Pass ``h_apo_km=None`` (or equal to ``h_peri_km``)
    to reproduce the legacy circular-parking-orbit period.
    """
    if h_apo_km is None:
        h_apo_km = h_peri_km
    rp = ts.R_EARTH + float(h_peri_km)
    ra = ts.R_EARTH + float(h_apo_km)
    a = 0.5 * (rp + ra)
    return 2.0 * math.pi * math.sqrt(a**3 / ts.MU_EARTH)


def parking_orbit_elements(
    h_peri_km: float = H_PARK_PERI_KM,
    h_apo_km: float = H_PARK_APO_KM,
) -> tuple[float, float, float]:
    """Return ``(a_km, e, n_rad_s)`` of the elliptical parking orbit.

    ``a`` is the semi-major axis (km), ``e`` the eccentricity, and ``n``
    the mean motion (rad/s) used by Kepler's equation in
    :func:`coast_from_argument`.
    """
    rp = ts.R_EARTH + float(h_peri_km)
    ra = ts.R_EARTH + float(h_apo_km)
    a = 0.5 * (rp + ra)
    e = (ra - rp) / (ra + rp)
    n = math.sqrt(ts.MU_EARTH / a**3)
    return a, e, n


def _parking_orbit_period_s(h_park_km: float = H_PARK_KM) -> float:
    """Legacy circular-parking-orbit period at altitude ``h_park`` (s)."""
    return parking_orbit_period_s(h_park_km, h_park_km)


def _mep_orbit_to_equator(
    family: pm.PerigeeFamily,
    i0: float,
    frame: mf.MEPFrame,
) -> tuple[float, float, float, float, np.ndarray]:
    """Convert a family entry (i0) to the parking orbit (i_E, Ω_E, ω_E, v_p, n̂).

    ``family`` returns (Ω_mep, ω_mep, v_p) at the given ``i0``.  These are
    defined in the frozen MEP frame at ``frame.t_arrive``; rotating back to
    J2000 (Earth-equatorial) via ``frame.rot_mep_to_j2000`` gives the
    classical elements of the parking orbit relative to the Earth equator.

    Parameters
    ----------
    family : PerigeeFamily
        L5 family interpolator (in MEP frame).
    i0 : float
        Family parameter (deg).
    frame : MEPFrame
        Frozen MEP frame at arrival.

    Returns
    -------
    i_E_deg : float
        Inclination of the parking orbit relative to the Earth equator,
        in [0°, 180°]; values above 90° are retrograde parking orbits.
    Omega_E_deg : float
        RAAN of the parking orbit (Earth equatorial J2000), measured to the
        ascending node (ẑ × n̂ direction).
    omega_E_deg : float
        Argument of perigee (Earth equatorial J2000).
    vp_kms : float
        Perigee velocity (km/s).
    n_hat : (3,) ndarray
        Unit normal of the parking orbit (J2000).  The sign is preserved
        (n_z < 0 for retrograde orbits); ``u_L`` branches select the
        ascending / descending pass.
    """
    Omega_mep = float(family.Omega_f(i0))
    omega_mep = float(family.omega_f(i0))
    vp = float(family.vp_f(i0))

    # Build the parking orbit's normal vector in the MEP frame from its
    # classical elements (MEP is the reference plane, so i_mep = i0):
    O = math.radians(Omega_mep)
    w = math.radians(omega_mep)
    i = math.radians(i0)
    # Normal in MEP = [sin i · sin Ω, −sin i · cos Ω, cos i]
    n_mep = np.array(
        [
            math.sin(i) * math.sin(O),
            -math.sin(i) * math.cos(O),
            math.cos(i),
        ],
        dtype=float,
    )
    n_mep /= np.linalg.norm(n_mep)
    # Rotate to J2000:
    n_hat = frame.rot_mep_to_j2000 @ n_mep
    n_hat /= np.linalg.norm(n_hat)

    # Inclination relative to the Earth equator (J2000 Z), [0°, 180°].
    # The sign of n_z MUST be preserved: folding with |n_z| would present
    # a retrograde parking orbit (i_E > 90°) as a prograde one and let it
    # slip through the prograde launch-azimuth corridor.
    inc_E = math.degrees(math.acos(max(-1.0, min(1.0, float(n_hat[2])))))

    # Ascending node direction = ẑ × n̂ (the orbit crosses the equator
    # going northward here).  n̂ × ẑ points to the *descending* node and
    # was the source of the historical 180° RAAN offset.
    z_hat = np.array([0.0, 0.0, 1.0])
    node = np.cross(z_hat, n_hat)
    node_norm = float(np.linalg.norm(node))
    if node_norm < 1e-12:
        # Equatorial orbit: RAAN is undefined; default to 0 and use the
        # J2000 X axis as the node reference.
        Omega_E = 0.0
        e_r = np.array([1.0, 0.0, 0.0])
    else:
        node /= node_norm
        e_r = node
        Omega_E = math.degrees(math.atan2(float(node[1]), float(node[0]))) % 360.0

    # Argument of perigee ω_E = angle, in the orbital plane, from the
    # ascending node to the perigee.  Compute the perigee direction in the
    # MEP frame from classical elements (Ω_mep, i_mep=i0, ω_mep), rotate to
    # J2000 via the frozen MEP frame, and project onto the J2000 perifocal
    # basis (e_r, e_w) to recover ω_E.
    O_mep = math.radians(Omega_mep)
    w_mep = math.radians(omega_mep)
    i_mep = math.radians(i0)
    r_peri_mep = np.array(
        [
            math.cos(w_mep) * math.cos(O_mep) - math.sin(w_mep) * math.sin(O_mep) * math.cos(i_mep),
            math.cos(w_mep) * math.sin(O_mep) + math.sin(w_mep) * math.cos(O_mep) * math.cos(i_mep),
            math.sin(w_mep) * math.sin(i_mep),
        ],
        dtype=float,
    )
    r_peri_mep /= np.linalg.norm(r_peri_mep)
    r_peri = frame.rot_mep_to_j2000 @ r_peri_mep

    # Right-handed perifocal basis (e_r × e_w = n̂) for the (Ω_E, inc_E)
    # node/inclination pair:
    e_w = np.cross(n_hat, e_r)
    e_w /= np.linalg.norm(e_w)

    # Project r_peri onto the orbital plane to compute ω_E.
    px = float(np.dot(r_peri, e_r))
    py = float(np.dot(r_peri, e_w))
    omega_E = math.degrees(math.atan2(py, px)) % 360.0

    return inc_E, Omega_E, omega_E, vp, n_hat


def _launch_azimuth_and_uL(
    i_E_deg: float,
    Omega_E_deg: float,
    site: LaunchSite,
    branch: str,
) -> tuple[float, float, float]:
    """Launch azimuth and launch argument-of-latitude for a parking orbit.

    Implements PRD v7 §五 constraint A:
        sin Az = cos i_E / cos φ_L   (standard launch azimuth, from North)
        sin u_L = sin φ_L / sin i_E

    Parameters
    ----------
    i_E_deg : float
        Parking-orbit inclination (deg, relative to Earth equator).
    Omega_E_deg : float
        Parking-orbit RAAN (deg).  Currently unused at this stage but kept
        in the signature for forward compatibility with the Δλ step.
    site : LaunchSite
        Launch site (WGS84).
    branch : {'ascending', 'descending'}
        ``'ascending'`` uses ``u_L`` (the site is on the ascending arc);
        ``'descending'`` uses ``180° − u_L``.

    Returns
    -------
    az_deg : float
        Launch azimuth (deg, 0–360; measured clockwise from North so that
        0° = North, 90° = East, 180° = South — the standard convention used
        by the R-2 / R-4 corridor tiers in ``perigee_match.AZIMUTH_TIERS``).
    u_L_deg : float
        Launch argument of latitude (deg, [0, 360]).
    alpha_L_deg : float
        Right ascension of the launch point (deg, [0, 360)).
    """
    phi = math.radians(site.lat_deg)
    i_E = math.radians(i_E_deg)

    sin_u = math.sin(phi) / math.sin(i_E)
    sin_u = max(-1.0, min(1.0, sin_u))
    u_pos = math.degrees(math.asin(sin_u))  # u_L in [0, +90]
    if branch == "ascending":
        u_L = u_pos
    elif branch == "descending":
        u_L = 180.0 - u_pos
    else:
        raise ValueError(f"branch must be 'ascending' or 'descending', got {branch!r}")

    # Standard launch azimuth (spherical triangle, from North clockwise):
    #   cos i_E = cos φ_L · sin Az
    # Ascending pass (site on the northward arc): Az = A; descending pass:
    # Az = 180° − A, with A = arcsin(cos i_E / cos φ_L) ∈ [−90°, 90°].
    # For retrograde orbits (i_E > 90°) A < 0, so the azimuth falls in
    # (180°, 360°) and is rejected by the prograde Wenchang corridor.
    sin_Az = math.cos(i_E) / math.cos(phi)
    sin_Az = max(-1.0, min(1.0, sin_Az))
    A_deg = math.degrees(math.asin(sin_Az))
    if branch == "ascending":
        az_deg = A_deg % 360.0
    else:
        az_deg = (180.0 - A_deg) % 360.0

    # Right ascension of the launch point: the direction of the point at
    # argument of latitude u_L on the orbit, expressed in the J2000
    # perifocal basis (e_r, e_w) of the (Ω_E, i_E) orbit.  Computing the
    # vector directly avoids the quadrant loss of the legacy
    # ``α_L = Ω_E + atan(tan u_L · cos i_E)`` form.
    cos_O = math.cos(math.radians(Omega_E_deg))
    sin_O = math.sin(math.radians(Omega_E_deg))
    e_r = np.array([cos_O, sin_O, 0.0])
    e_w = np.array(
        [-sin_O * math.cos(i_E), cos_O * math.cos(i_E), math.sin(i_E)]
    )
    r_L = (
        math.cos(math.radians(u_L)) * e_r
        + math.sin(math.radians(u_L)) * e_w
    )
    alpha_L = math.degrees(math.atan2(float(r_L[1]), float(r_L[0]))) % 360.0
    return az_deg, u_L, alpha_L


def _apply_rev_extension(
    coast_s: float, period_s: float, extend_revs: bool
) -> float:
    """Add whole parking-orbit revolutions to lift ``coast_s`` to bounds.

    When ``extend_revs`` is True, whole periods are added until the coast
    is at least ``COAST_MIN_S`` (capped at 5 revolutions).  Under the F0c
    policy ``COAST_MIN_S = 0.0``, so this is a **no-op**: short and zero
    coasts are kept as-is and no long-coast branch is generated.  The
    helper is retained for callers/tests that pass a positive
    ``COAST_MIN_S``.
    """
    if not extend_revs:
        return coast_s
    # For the 200 × 380 km parking orbit the period is ~5418 s, so the
    # loop only runs once or twice.
    n_rev = 0
    while coast_s < COAST_MIN_S and n_rev < 5:
        coast_s += period_s
        n_rev += 1
    return coast_s


def insertion_argument_deg(
    u_L_deg: float, ascent_arc_deg: float = ASCENT_ARC_DEG
) -> float:
    """Parking-orbit insertion argument of latitude (deg).

    ``u_ins = u_L + ascent_arc_deg`` (mod 360): the spacecraft is inserted at
    the launch-site argument of latitude ``u_L`` advanced by the ascent
    geocentric arc.  The argument of latitude increases along the direction
    of motion for both the ascending and descending branch, so the same
    ``+arc`` sign applies to either branch (the branch only selects
    ``u_L = u_pos`` or ``u_L = 180° − u_pos``).
    """
    return (float(u_L_deg) + float(ascent_arc_deg)) % 360.0


def _true_to_mean_anomaly(f_rad: float, e: float) -> float:
    """Mean anomaly (rad, [0, 2π)) for true anomaly ``f_rad`` on an ellipse.

    Solves Kepler's equation via the half-angle form
    ``tan(E/2) = sqrt((1−e)/(1+e)) · tan(f/2)`` followed by
    ``M = E − e·sin E``.
    """
    E = 2.0 * math.atan2(
        math.sqrt(1.0 - e) * math.sin(0.5 * f_rad),
        math.sqrt(1.0 + e) * math.cos(0.5 * f_rad),
    )
    E %= 2.0 * math.pi
    return (E - e * math.sin(E)) % (2.0 * math.pi)


def coast_from_argument(
    omega_E_deg: float,
    u_L_deg: float,
    h_peri_km: float = H_PARK_PERI_KM,
    h_apo_km: float = H_PARK_APO_KM,
) -> float:
    """Coast from the parking-orbit insertion point to TLI (s).

    The TLI burn is performed at the parking-orbit perigee (true anomaly
    ``f = 0`` — where the parking-orbit speed is highest, minimising the
    TLI Δv).  The insertion point sits at argument of latitude ``u_L``,
    i.e. true anomaly ``f_L = (u_L − ω_E) mod 360°``.  On an ellipse the
    elapsed time to the next perigee is the mean-anomaly difference
    divided by the mean motion (Kepler's equation), NOT a simple
    ``Δu/360° · P`` fraction:

        E_L = 2·atan2(sqrt(1−e)·sin(f_L/2), sqrt(1+e)·cos(f_L/2))
        M_L = E_L − e·sin(E_L)
        coast = (2π − M_L) / n

    For ``h_peri_km == h_apo_km`` (``e = 0``) this reduces exactly to the
    legacy circular formula ``((ω_E − u_L) mod 360°)/360° · P``.

    Parameters
    ----------
    omega_E_deg : float
        Parking-orbit argument of perigee (deg).
    u_L_deg : float
        Launch (insertion) argument of latitude (deg).
    h_peri_km, h_apo_km : float
        Parking-orbit perigee / apogee altitudes (km).  Defaults are the
        CZ-5 200 × 380 km ellipse.

    Returns
    -------
    coast_s : float
        Raw coast time, in ``[0, period)`` (no whole-revolution
        extension; use :func:`parking_orbit_coast_seconds` for that).
    """
    _a, e, n = parking_orbit_elements(h_peri_km, h_apo_km)
    f_launch = math.radians((u_L_deg - omega_E_deg) % 360.0)
    M_launch = _true_to_mean_anomaly(f_launch, e)
    return ((2.0 * math.pi - M_launch) % (2.0 * math.pi)) / n


def parking_orbit_coast_seconds(
    omega_E_deg: float,
    u_L_deg: float,
    extend_revs: bool = True,
    h_peri_km: float = H_PARK_PERI_KM,
    h_apo_km: float = H_PARK_APO_KM,
) -> float:
    """Coast (s) on the elliptical parking orbit, with rev extension.

    Computes the raw Kepler-equation coast via :func:`coast_from_argument`
    and, when ``extend_revs`` is True, adds whole parking-orbit
    revolutions until the coast is at least ``COAST_MIN_S`` (see
    :func:`_apply_rev_extension`).
    """
    coast_s = coast_from_argument(omega_E_deg, u_L_deg, h_peri_km, h_apo_km)
    period_s = parking_orbit_period_s(h_peri_km, h_apo_km)
    return _apply_rev_extension(coast_s, period_s, extend_revs)


def _coast_seconds(
    omega_E_deg: float,
    u_L_deg: float,
    period_s: float,
    extend_revs: bool,
) -> float:
    """Legacy circular-orbit coast from launch point to TLI (s).

    Kept for API compatibility: this is the ``e = 0`` special case of
    :func:`coast_from_argument`.  Production callers use
    :func:`parking_orbit_coast_seconds` (elliptical CZ-5 parking orbit);
    this helper is retained for callers/tests that pass an explicit
    circular period.

    Parameters
    ----------
    omega_E_deg, u_L_deg : float
        Argument of perigee and launch argument-of-latitude (deg).
    period_s : float
        Circular parking-orbit Keplerian period (s).
    extend_revs : bool
        If True, add whole revolutions so the coast falls within
        ``[COAST_MIN_S, COAST_MAX_S]`` whenever possible.  If False, return
        the raw (mod period) coast, which is in ``[0, period_s)``.

    Returns
    -------
    coast_s : float
        Coast time (s).  In ``[0, period_s)`` if ``extend_revs`` is False;
        otherwise biased to be in ``[COAST_MIN_S, COAST_MAX_S]`` when
        possible (if even the smallest positive coast violates the bounds
        the smallest positive value is returned).
    """
    d_theta = (omega_E_deg - u_L_deg) % 360.0
    coast = d_theta / 360.0 * period_s
    return _apply_rev_extension(coast, period_s, extend_revs)


# ---------------------------------------------------------------------------
# Inner fixed-point: i0 ↔ i_E
# ---------------------------------------------------------------------------
def _solve_inner_fixed_point(
    family: pm.PerigeeFamily,
    t_arrive: float,
    T_s: float,
    dtheta_deg: float,
    site: LaunchSite,
    i0_lo: float,
    i0_hi: float,
    t_L_override: Optional[float] = None,
    ascent_arc_deg: float = ASCENT_ARC_DEG,
    branch: str = "ascending",
) -> Optional[tuple[float, float, float, float, float, np.ndarray]]:
    """Legacy single-return wrapper around :func:`_solve_inner_fixed_point_with_code`.

    Preserves the public-facing single-return shape used by the
    ``refine_ce7_2026`` / ``refine_ce7_2026_T112_113`` modules — those
    callers ``if inner is not None: _i0, _inc_E, ... = inner`` against
    the bare 6-tuple, not the new ``(result, fail_code)`` tuple.
    """
    result, _code = _solve_inner_fixed_point_with_code(
        family=family,
        t_arrive=t_arrive,
        T_s=T_s,
        dtheta_deg=dtheta_deg,
        site=site,
        i0_lo=i0_lo,
        i0_hi=i0_hi,
        t_L_override=t_L_override,
        ascent_arc_deg=ascent_arc_deg,
        branch=branch,
    )
    return result


def _solve_inner_fixed_point_with_code(
    family: pm.PerigeeFamily,
    t_arrive: float,
    T_s: float,
    dtheta_deg: float,
    site: LaunchSite,
    i0_lo: float,
    i0_hi: float,
    t_L_override: Optional[float] = None,
    ascent_arc_deg: float = ASCENT_ARC_DEG,
    branch: str = "ascending",
) -> tuple[
    Optional[tuple[float, float, float, float, float, np.ndarray]],
    Optional[str],
]:
    """Run the i0 ↔ i_E inner fixed point at a frozen ``t_arrive``.

    Reuses ``pm.solve_i0_for_family`` to find the self-consistent ``i0``
    that matches the parking-orbit-plane inclination at the launch site,
    then converts (Ω, ω, v_p) to the Earth-equatorial (i_E, Ω_E, ω_E)
    representation.

    The launch site is evaluated at the *estimated launch epoch*
    ``t_L_est``:

    - ``t_L_override`` given (optional override): the site is placed at
      the caller-supplied epoch.  NOTE: feeding the GMST-recovered launch
      epoch back in this way was tried as a self-consistency iteration and
      rejected — the consistency residual ``t_L' − t_L`` has no root for
      the CE-5/CE-6 family tables (it stays ≈ +4 … +40 min), so the
      iteration wanders off.  The production path leaves this ``None``.
    - otherwise: ``t_L_est = t_arrive − T_s − dT_s`` with
      ``dT_s = dtheta_deg / 15.041 · 3600`` (the R-2 B-7 geocentric-arc
      convention; dtheta should be the measured launch→perigee angle,
      see ``verify_ce56.mission_dtheta_deg``).

    Returns ``(result, fail_code)`` where ``result`` is either the
    6-tuple ``(i0, inc_E, Omega_E, omega_E, vp, n_hat)`` on success or
    ``None`` on failure, and ``fail_code`` is a string from
    ``_FAIL_CODES`` on failure or ``None`` on success.  The fail-code
    vocabulary is:

    - ``"dtheta_exceeded"``: ``dtheta_deg > _DTHETA_MAX_DEG``.
    - ``"i0_out_of_family_range"``: ``lo >= hi`` (no overlap between
      caller-provided range and family table), or no sign-change
      detected in the i0 sign-change scan.
    - ``"inc_below_site_lat"``: ``inc_E < φ_ins`` (plane below the
      parking-orbit *insertion point* reachable band; see F0b note below).
    - ``"inc_above_supplement"``: ``inc_E > 180° − φ_ins`` (plane
      retrograde beyond the insertion point's supplement).

    F0b envelope convention (2026-09): the reachable inclination band is
    evaluated at the parking-orbit **insertion point**, not the launch-site
    geodetic latitude.  The insertion sits ``ascent_arc_deg`` downrange of
    the site, so its geocentric latitude is
    ``φ_ins = asin(sin i_E · sin u_ins)`` with
    ``u_ins = u_L + ascent_arc_deg``.  For a plane that contains the launch
    site (guaranteed by the inner fixed point) this band is structurally
    satisfied in the prograde direction; the check is retained to reject
    degenerate/retrograde cases and to keep the fail-code vocabulary stable.
    """
    # Pre-check: caller-supplied dtheta must be within the engine's
    # geometric envelope.  Beyond ~180° the inbound/outbound plane
    # mismatch is qualitatively different and the fixed-point no longer
    # converges to a physically meaningful family entry.
    if dtheta_deg > _DTHETA_MAX_DEG:
        return None, "dtheta_exceeded"

    frame = mf.mep_at_epoch(t_arrive)
    if t_L_override is not None:
        # Optional override: evaluate the site at the supplied epoch.
        t_L_est = float(t_L_override)
    else:
        # Initial guess: the R-2 B-7 convention
        # ``dT = dtheta / 15.041 * 3600`` (geocentric arc from launch to
        # perigee expressed as an Earth-rotation time).
        dT_s = dtheta_deg / EARTH_SIDEREAL_RATE_DEG_PER_H * 3600.0
        t_L_est = t_arrive - T_s - dT_s
    r_site = g.launch_site_vector_j2000_km(
        site.lon_deg, site.lat_deg, site.alt_m, t_L_est
    )
    r_site_mep = mf.to_mep(r_site, frame)
    r_site_mep = pm.R_P * r_site_mep / np.linalg.norm(r_site_mep)

    lo = max(i0_lo, family.i0_min)
    hi = min(i0_hi, family.i0_max)
    if lo >= hi:
        return None, "i0_out_of_family_range"

    def f(i0: float) -> float:
        lon_p = float(family.lon_f(i0))
        lat_p = float(family.lat_f(i0))
        r_peri = mf.from_mep_lonlat(lon_p, lat_p, pm.R_P)
        return pm.parking_orbit_inclination_deg(r_site_mep, r_peri) - i0

    from scipy.optimize import brentq

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
    if not roots:
        return None, "i0_out_of_family_range"

    # Convert the root to (i_E, Ω_E, ω_E, v_p, n̂) via _mep_orbit_to_equator.
    i0 = roots[0]
    inc_E, Omega_E, omega_E, vp, n_hat = _mep_orbit_to_equator(family, i0, frame)

    # Envelope check at the parking-orbit insertion point (F0b; split into
    # the two physically distinct failure modes so the reject record can
    # carry the right ``fail_code``).  ``u_ins`` is the insertion argument
    # of latitude and the insertion geocentric latitude is
    # ``φ_ins = asin(sin i_E · sin u_ins)``; the launch-site geodetic
    # latitude is deliberately NOT used any more (it over-rejects planes
    # whose downrange insertion lies at a lower latitude, cf. Wenchang
    # i_E ≈ 19.5° < φ_gd ≈ 19.614°).
    _az_env, u_L_env, _alpha_env = _launch_azimuth_and_uL(
        inc_E, Omega_E, site, branch=branch
    )
    u_ins_env = insertion_argument_deg(u_L_env, ascent_arc_deg)
    sin_phi_ins = abs(
        math.sin(math.radians(inc_E)) * math.sin(math.radians(u_ins_env))
    )
    phi_ins = math.degrees(
        math.asin(max(-1.0, min(1.0, sin_phi_ins)))
    )
    if inc_E < phi_ins:
        return None, "inc_below_site_lat"
    if inc_E > 180.0 - phi_ins:
        return None, "inc_above_supplement"

    return (i0, inc_E, Omega_E, omega_E, vp, n_hat), None


# ---------------------------------------------------------------------------
# Lambert-guess helpers (REQ-LW-2026-LAMBERT-IG-001 §六 / §八)
# ---------------------------------------------------------------------------
def _r1_at_tli_km(
    inc_E_deg: float,
    Omega_E_deg: float,
    omega_E_deg: float,
    r_park_km: float,
) -> np.ndarray:
    """TLI-side position vector (J2000) on the parking orbit at perigee.

    The TLI point is, by construction (PRD v7 §五), at the parking-orbit
    perigee.  Given ``(i_E, Ω_E, ω_E)`` the J2000 perigee direction is the
    first row of the classical-elements rotation matrix; multiplying by
    ``r_park`` (km) gives r1 in km.  REQ §八.3 / R4 / R6: this is the
    correct r1 to feed Lambert — NOT the launch-site vector.
    """
    inc = math.radians(inc_E_deg)
    Om = math.radians(Omega_E_deg)
    w = math.radians(omega_E_deg)
    cos_O = math.cos(Om)
    sin_O = math.sin(Om)
    cos_i = math.cos(inc)
    sin_i = math.sin(inc)
    cos_w = math.cos(w)
    sin_w = math.sin(w)
    return np.array(
        [
            r_park_km * (cos_O * cos_w - sin_O * sin_w * cos_i),
            r_park_km * (sin_O * cos_w + cos_O * sin_w * cos_i),
            r_park_km * (sin_w * sin_i),
        ],
        dtype=float,
    )


def _lambert_parking_elements(
    *,
    t_arrive: float,
    T_s: float,
    family: pm.PerigeeFamily,
    site: LaunchSite,
    branch: str,
    dtheta_deg: float,
    i0_lo: float,
    i0_hi: float,
    t_L_seed: float,
    t_ascent_s: float,
    ascent_arc_deg: float = ASCENT_ARC_DEG,
) -> tuple[Optional[tuple[float, float, float, float]], Optional[str]]:
    """Lambert-driven replacement for the family-table inner link.

    Returns ``(result, fail_code)`` where ``result`` is
    ``(inc_E_deg, Omega_E_deg, omega_E_deg, vp_kms)`` on success or
    ``None`` on failure.  ``fail_code`` is the string tag identifying
    which sub-step failed (always ``"lambert_fail"`` at present — kept
    separate from ``_solve_inner_fixed_point_with_code``'s vocabulary so
    the L6 outer-solver reject records can distinguish Lambert from
    fixed-point failures).

    REQ §六 / §八.3: invoked inside the ``_G_te`` self-consistency loop
    using the current best ``t_TLI`` estimate; on ``RuntimeError`` from
    ``lambert_initial_guess`` (REQ R1 / R2) we log a warning and return
    ``(None, "lambert_fail")``.
    """
    # Seed t_TLI estimate from family-table path so r1 is well-defined
    # even on the very first iteration (REQ §十三.2 — "先用 family-table
    # 算一次 t_TLI^0 作为初值").  The site epoch follows the same R-2
    # dtheta convention as ``_G_te`` (see the note there on why the
    # recovered-t_L feedback was rejected).
    inner, inner_code = _solve_inner_fixed_point_with_code(
        family=family,
        t_arrive=t_arrive,
        T_s=T_s,
        dtheta_deg=dtheta_deg,
        site=site,
        i0_lo=i0_lo,
        i0_hi=i0_hi,
        ascent_arc_deg=ascent_arc_deg,
        branch=branch,
    )
    if inner is None:
        # Propagate the inner-fixed-point fail code through so the
        # outer reject record captures the *root* cause (envelope /
        # sign-change / dtheta), not a generic "lambert_fail" that would
        # mis-attribute the failure to the Lambert step.
        return None, inner_code or "lambert_fail"
    _i0, inc_E0, Omega_E0, omega_E0, vp0, _n_hat0 = inner

    # Build r1 from the family-table parking elements + t_TLI estimate.
    # coast / t_TLI use the same convention as _G_te's tail (coast starts at
    # the insertion argument u_ins = u_L + ascent arc).
    az_seed, u_L_seed, _alpha = _launch_azimuth_and_uL(
        inc_E0, Omega_E0, site, branch
    )
    coast_seed = parking_orbit_coast_seconds(
        omega_E0, insertion_argument_deg(u_L_seed, ascent_arc_deg),
        extend_revs=False,
    )
    t_TLI_seed = t_L_seed + t_ascent_s + coast_seed

    r_park_km = ts.R_EARTH + H_PARK_KM
    r1 = _r1_at_tli_km(inc_E0, Omega_E0, omega_E0, r_park_km)
    r2 = np.asarray(g.moon_position_j2000_km(t_arrive), dtype=float)
    dt_s = t_arrive - t_TLI_seed

    try:
        v1 = lambert_mod.lambert_initial_guess(
            r1_eci_km=r1,
            r2_eci_km=r2,
            dt_s=float(dt_s),
            mu_earth_km3s2=ts.MU_EARTH,
            prograde=True,
            low_path=True,
        )
    except (RuntimeError, ValueError) as exc:
        _log.warning(
            "lambert_initial_guess failed at Te=%s (family=%s branch=%s): %s; "
            "falling back to family-table path",
            f"{t_arrive:.3f}",
            family.name,
            branch,
            exc,
        )
        return None, "lambert_fail"

    try:
        i_park, raan, omega, _a_km = lambert_mod.classical_from_state(
            r1, v1, ts.MU_EARTH
        )
    except ValueError as exc:
        _log.warning(
            "classical_from_state failed at Te=%s (family=%s branch=%s): %s; "
            "falling back to family-table path",
            f"{t_arrive:.3f}",
            family.name,
            branch,
            exc,
        )
        return None, "lambert_fail"

    # Envelope check (REQ §三: 不修改几何约束 A/B/C).  Sub-codes mirror
    # the inner-fixed-point vocabulary so callers can deduplicate by
    # envelope vs Lambert-only failure modes.
    phi = abs(site.lat_deg)
    if i_park < phi:
        return None, "inc_below_site_lat"
    if i_park > 180.0 - phi:
        return None, "inc_above_supplement"

    # v_p equivalent: the Lambert v1 is at r1 (TLI point), which IS the
    # parking-orbit perigee by construction.  Its magnitude is the
    # perigee velocity on the parking orbit.  This makes vp_kms
    # numerically interchangeable with the family-table vp_kms for
    # downstream C3 / window bookkeeping (REQ R7).
    vp_kms = float(np.linalg.norm(v1))
    return (float(i_park), float(raan) % 360.0, float(omega) % 360.0, vp_kms), None


# ---------------------------------------------------------------------------
# Outer residual G(T_e)
# ---------------------------------------------------------------------------
def _G_te(
    t_arrive: float,
    T_s: float,
    family: pm.PerigeeFamily,
    site: LaunchSite,
    branch: str,
    dtheta_deg: float,
    i0_lo: float,
    i0_hi: float,
    extend_revs: bool,
    *,
    return_full: bool = False,
    t_ascent_s: float = T_ASCENT_S,
    coast_max_s: Optional[float] = None,
    ascent_arc_deg: float = ASCENT_ARC_DEG,
) -> tuple[
    float, Optional[float], Optional[float], Optional[float],
    Optional[float], Optional[float], Optional[float], Optional[str],
]:
    """Compute the outer residual ``G(t_arrive) = t_TLI + T − t_arrive``.

    The inner fixed point runs at a frozen ``t_arrive`` to give the
    parking-orbit plane (``inc_E, Omega_E, omega_E, vp``).  The launch
    azimuth and argument-of-latitude ``u_L`` are derived from
    ``inc_E`` and the launch site.  ``GMST(t_L) = alpha_L − lon_site``
    gives the launch epoch ``t_L``.  The coast starts at the insertion
    argument ``u_ins = u_L + ascent_arc_deg`` (``ascent_arc_deg=0``
    reproduces the legacy site-as-insertion convention).  Finally
    ``t_TLI = t_L + T_ASCENT_S + coast`` and
    ``G = (t_TLI + T) − t_arrive``.

    t_L self-consistency: the inner fixed point initialises with a
    ``t_L_est = t_arrive − T − dT`` (using the geocentric arc ``dT``
    in lieu of the geodetic longitude).  Once ``t_L`` is recovered
    from GMST, the launch-site vector changes by ~16 m/s-of-arc per
    17 s of ``t_L`` drift; we re-evaluate the inner fixed point with
    the new ``t_L`` until ``|Δt_L| ≤ 1 s`` (max 4 iterations).  The
    residual ~17 s observed in earlier rounds came from this
    inconsistency, not from the outer solve.

    Returns the 8-tuple
    ``(G, t_L, az_deg, u_L_deg, inc_E, omega_E, vp, fail_code)``.
    On success ``fail_code`` is ``None`` and the other fields carry the
    converged geometry.  On any failure (inner-fixed-point, GMST
    inverse, Lambert guess) ``G`` is the sentinel ``1e9``, the
    auxiliary fields are ``None`` and ``fail_code`` is a string from
    ``_FAIL_CODES`` (e.g. ``"inc_below_site_lat"``, ``"gmst_fail"``,
    ``"lambert_fail"``) identifying which sub-step aborted.
    """
    SENTINEL = 1e9

    def _one_pass(t_arrive_x: float, t_L_seed: Optional[float]):
        inner, inner_code = _solve_inner_fixed_point_with_code(
            family=family,
            t_arrive=t_arrive_x,
            T_s=T_s,
            dtheta_deg=dtheta_deg,
            site=site,
            i0_lo=i0_lo,
            i0_hi=i0_hi,
            ascent_arc_deg=ascent_arc_deg,
            branch=branch,
        )
        if inner is None:
            return None, inner_code
        _i0, inc_E, Omega_E, omega_E, vp, _n_hat = inner
        az_deg, u_L_deg, alpha_L_deg = _launch_azimuth_and_uL(
            inc_E, Omega_E, site, branch
        )
        target_gmst = (alpha_L_deg - site.lon_deg) % 360.0
        # If we have a t_L seed, use it to anchor gmst_to_et closer to
        # the true root.
        guess = t_L_seed if t_L_seed is not None else t_arrive_x - T_s - t_ascent_s
        try:
            t_L = gmst_to_et(target_gmst, guess)
        except ValueError:
            return None, "gmst_fail"
        return (inc_E, Omega_E, omega_E, vp, az_deg, u_L_deg, t_L), None

    # First pass with the R-2 dtheta estimate for the launch site.
    # NOTE (self-consistency attempt, 2026-09): feeding the recovered t_L
    # back into the inner fixed point (``t_L_override``) was tried and
    # rejected: the consistency residual h(t_L) = t_L' − t_L is strictly
    # positive over the feasible site-epoch range (h ≈ +4 … +40 min) for
    # the CE-5/CE-6 family tables — the fixed point does not exist, so the
    # naive iteration wanders eastward until ``inc_below_site_lat``.  The
    # R-2 B-7 ``dtheta`` is therefore a *physical bias term* (launch →
    # perigee time, dtheta = Δt_L·ω_E), not a removable initial guess;
    # callers should set it from the measured ascent+coast timeline.
    out, out_code = _one_pass(t_arrive, None)
    if out is None:
        return (SENTINEL, None, None, None, None, None, None, out_code)
    inc_E, Omega_E, omega_E, vp, az_deg, u_L_deg, t_L = out

    # REQ-LW-2026-LAMBERT-IG-001 §八.3: if the master switch is on, replace
    # the family-table-derived (i_E, Ω_E, ω_E, vp) with a Lambert solve
    # using the converged t_L.  On any failure (geometry degeneracy,
    # izzo non-convergence) the helper returns ``(None, fail_code)`` and
    # we keep the family-table values bit-identically (REQ §九.1.4).
    # The ``fail_code`` is captured here for diagnostic propagation only;
    # the outer residual ``G`` is still well-defined with the family
    # values, so we do NOT return the sentinel in that case.
    if USE_LAMBERT_GUESS:
        lam, _lam_code = _lambert_parking_elements(
            t_arrive=t_arrive,
            T_s=T_s,
            family=family,
            site=site,
            branch=branch,
            dtheta_deg=dtheta_deg,
            i0_lo=i0_lo,
            i0_hi=i0_hi,
            t_L_seed=t_L,
            t_ascent_s=t_ascent_s,
            ascent_arc_deg=ascent_arc_deg,
        )
        if lam is not None:
            inc_E, Omega_E, omega_E, vp = lam

    coast = parking_orbit_coast_seconds(
        omega_E, insertion_argument_deg(u_L_deg, ascent_arc_deg),
        extend_revs=extend_revs,
    )
    # Solver-level coast cap: the parking-orbit coast is part of the launch
    # timeline (it decides the TLI epoch), so a timeline constraint must be
    # applied *inside* the G evaluation — invalidating the grid point /
    # pocket rather than post-filtering finished windows (a post-filter can
    # hide a valid short-coast root that lives elsewhere in the pocket).
    if coast_max_s is not None and coast > float(coast_max_s):
        return (SENTINEL, t_L, az_deg, u_L_deg, inc_E, omega_E, vp,
                "coast_exceeds_limit")
    t_TLI = t_L + t_ascent_s + coast
    G = (t_TLI + T_s) - t_arrive
    if return_full:
        return (G, t_L, az_deg, u_L_deg, inc_E, omega_E, vp, None)
    return (G, t_L, az_deg, u_L_deg, inc_E, omega_E, vp, None)


def solve_epoch(
    Te_guess: float,
    site: LaunchSite,
    T_h: float,
    family: pm.PerigeeFamily,
    polar_ref: str = "mep",
    dtheta_deg: float = 25.0,
    i0_lo: float = 19.6,
    i0_hi: float = 114.0,
    branch: str = "ascending",
    extend_revs: bool = True,
    i0_init: Optional[float] = None,
    H_m: float = 200.0,
    t_ascent_s: float = T_ASCENT_S,
    rejects: Optional[list[WindowReject]] = None,
    te0_et: Optional[float] = None,
    coast_max_s: Optional[float] = None,
    ascent_arc_deg: float = ASCENT_ARC_DEG,
) -> Optional[Window]:
    """Solve the L6 outer problem ``G(T_e) = t_TLI + T − T_e = 0`` for one combo.

    Algorithm: pocket-search root finding.

    The feasible domain of ``T_e`` consists of 1–7 h wide discrete *pockets*
    separated by ``+86164 s`` (sidereal-day) gaps where the inner fixed
    point returns ``None``.  Inside each pocket ``G`` is continuous with
    ``dG/dT_e ≈ −0.96 s/s`` (the 1 − ω_moon/ω_earth factor); each pocket
    therefore contains at most one root.  The GMST anchor at the
    sidereal-day boundary creates a ``+86164 s`` apparent sign jump that
    is a *discontinuity artifact*, not a true root: such pseudo-roots are
    discarded by the ``|G| ≤ G_TOL_S`` test after Brent refinement.

    Implementation
    --------------
    1. Build a coarse grid ``[Te_guess − 12 h, Te_guess + 12 h]`` at
       ``_POCKET_STEP_S`` (0.5 h) and evaluate ``G`` at every point,
       dropping any point where the inner fixed point fails.
    2. Adjacent grid points that are *both* defined AND have opposite
       signs → refine with ``brentq`` (xtol=1e-6 s) to ``|G| ≤ G_TOL_S``.
       Any bracket that does not refine below ``G_TOL_S`` (a
       discontinuity pseudo-root) is dropped.
    3. Among the qualifying roots pick the one closest to ``Te_guess``
       and build the ``Window``.
    4. If no root is found, expand the search to
       ``[Te_guess − 24 h, Te_guess + 24 h]`` once (the spec allows
       one expansion).  Still nothing → return ``None``.

    ``DG_DTE_HINT`` (≈ 0.967) is documented as a side-information
    constant only — the grid + Brent approach subsumes the analytical
    step.

    Parameters
    ----------
    Te_guess : float
        Initial guess for the arrival epoch (SPICE ET).
    site : LaunchSite
        Launch site (WGS84).
    T_h : float
        Transfer time (hours).
    family : PerigeeFamily
        L5 perigee-family interpolator for one solution family (e.g.
        ``"S_asc"``, ``"N_desc"``, ...).
    polar_ref : {'mep', 'moon_equator'}
        Polar reference axis (PRD v7 §四.1).  Recorded on the Window.
    dtheta_deg : float
        Geocentric arc launch→perigee (deg); used to initialise the
        inner fixed-point estimate.
    i0_lo, i0_hi : float
        Search range for the family parameter (deg).
    branch : {'ascending', 'descending'}
        Branch to solve.
    extend_revs : bool
        If True, add whole parking-orbit revolutions to the coast when
        it would otherwise be below ``COAST_MIN_S``.
    i0_init : float, optional
        Initial family-parameter guess (deg).  Currently unused by the
        pocket search; kept for API compatibility.
    H_m : float
        Target perilune altitude (km).  Recorded on the Window.

    Returns
    -------
    Window or None
        The converged window record, or ``None`` if no root is found.
    """
    if polar_ref not in ("mep", "moon_equator"):
        raise ValueError(
            f"polar_ref must be 'mep' or 'moon_equator', got {polar_ref!r}"
        )
    if branch not in ("ascending", "descending"):
        raise ValueError(
            f"branch must be 'ascending' or 'descending', got {branch!r}"
        )

    T_s = T_h * 3600.0
    Te0 = float(Te_guess)
    if te0_et is None:
        te0_et = Te0

    # Two-stage pocket search: first ±12 h, then ±24 h if empty.
    for radius_s in (_POCKET_RADIUS_S, _POCKET_EXPAND_RADIUS_S):
        win = _pocket_search(
            Te0=Te0,
            radius_s=radius_s,
            T_s=T_s,
            family=family,
            site=site,
            branch=branch,
            dtheta_deg=dtheta_deg,
            i0_lo=i0_lo,
            i0_hi=i0_hi,
            extend_revs=extend_revs,
            polar_ref=polar_ref,
            H_m=H_m,
            t_ascent_s=t_ascent_s,
            rejects=rejects,
            te0_et=te0_et,
            coast_max_s=coast_max_s,
            ascent_arc_deg=ascent_arc_deg,
        )
        if win is not None:
            return win
    # Both pocket stages exhausted with no qualifying root.
    if rejects is not None:
        lo = max(i0_lo, family.i0_min)
        hi = min(i0_hi, family.i0_max)
        i0_deg_val = lo if lo < hi else None
        # Enrich the ``no_root_found`` record with the moon ephemeris at
        # the seed arrival epoch and the geometric envelope bounds, so
        # the diagnostic tells the caller exactly which Earth-Moon
        # geometry blocked the launch opportunity.
        try:
            moon_pos = np.asarray(g.moon_position_j2000_km(float(te0_et)), dtype=float)
            dist_km = float(np.linalg.norm(moon_pos))
            dec_deg = float(np.degrees(np.arcsin(moon_pos[2] / dist_km)))
            ra_deg = float(np.degrees(np.arctan2(moon_pos[1], moon_pos[0])) % 360.0)
            moon_dec, moon_ra, moon_dist = dec_deg, ra_deg, dist_km
        except Exception:
            moon_dec = moon_ra = moon_dist = float("nan")
        envelope = abs(site.lat_deg)
        rejects.append(_make_reject(
            te0_et=float(te0_et),
            T_h=T_h,
            family=family.name,
            branch=branch,
            reason="no_root_found",
            i0_deg=i0_deg_val,
            extra={
                "pocket_radius_s": _POCKET_EXPAND_RADIUS_S,
                "i0_lo": float(i0_lo),
                "i0_hi": float(i0_hi),
                "i0_min": float(family.i0_min),
                "i0_max": float(family.i0_max),
                "moon_dec_deg": moon_dec,
                "moon_ra_deg": moon_ra,
                "earth_moon_dist_km": moon_dist,
                "inc_e_min_allowed": envelope,
                "inc_e_max_allowed": 180.0 - envelope,
                "dtheta_actual_deg": float(dtheta_deg),
            },
        ))
    return None


# Pocket-search grid parameters (PRD v7 §五 / 工作1 spec).
_POCKET_STEP_S = 1800.0           # 0.5 h grid step
_POCKET_RADIUS_S = 12.0 * 3600.0  # ±12 h search radius (primary)
_POCKET_EXPAND_RADIUS_S = 24.0 * 3600.0  # ±24 h fallback (single expansion)
_BRENT_XTOL_S = 1e-6              # brentq xtol (≪ 1 s of GMST arc)


def _pocket_search(
    *,
    Te0: float,
    radius_s: float,
    T_s: float,
    family: pm.PerigeeFamily,
    site: LaunchSite,
    branch: str,
    dtheta_deg: float,
    i0_lo: float,
    i0_hi: float,
    extend_revs: bool,
    polar_ref: str,
    H_m: float,
    t_ascent_s: float = T_ASCENT_S,
    rejects: Optional[list[WindowReject]] = None,
    te0_et: Optional[float] = None,
    coast_max_s: Optional[float] = None,
    ascent_arc_deg: float = ASCENT_ARC_DEG,
) -> Optional[Window]:
    """Pocket-search root finder on ``[Te0 − radius, Te0 + radius]``.

    Returns the ``Window`` of the root closest to ``Te0`` (with
    ``|G| ≤ G_TOL_S``), or ``None`` if no qualifying root exists in
    the search interval.

    When ``rejects`` is provided, per-rejection ``WindowReject`` records
    are appended at each drop point (no sign change, brentq failure,
    discontinuity artifact, ``_G_te`` sentinel).  Each reject record
    carries an enriched ``extra`` field with the moon ephemeris at the
    arrival epoch and the geometric envelope bounds
    (``inc_e_min/max_allowed``, ``dtheta_actual_deg``) so callers can
    debug which Earth-Moon geometry combination caused the rejection.

    ``te0_et`` is the seed ``Te0`` passed by the caller for context;
    defaults to ``Te0``.
    """
    from scipy.optimize import brentq as _brentq

    SENTINEL_THRESHOLD = 1e8  # values >= this count as _G_te failures

    if te0_et is None:
        te0_et = float(Te0)

    # Build the coarse grid.
    n_half = int(radius_s / _POCKET_STEP_S)
    xs = np.arange(Te0 - n_half * _POCKET_STEP_S, Te0 + n_half * _POCKET_STEP_S + 0.5 * _POCKET_STEP_S, _POCKET_STEP_S)

    # Evaluate G at every grid point (caching inner-failure → None).
    n = len(xs)
    G_vals: list[Optional[float]] = [None] * n
    t_L_vals: list[Optional[float]] = [None] * n
    aux_vals: list = [None] * n  # (az, uL, inc_E, omega_E, vp)
    # Parallel fail-code array populated by ``_G_te``'s 8th return
    # element.  ``None`` means "no failure recorded here" (which is
    # either success or a non-tagged drop such as no-sign-change).
    fail_codes: list[Optional[str]] = [None] * n

    def _moon_geom(et: float) -> dict[str, float]:
        """Compute moon RA/Dec/distance (J2000) at ``et`` (SPICE ET).

        Used to enrich every ``WindowReject.extra`` field with the
        geometric context that explains why a launch opportunity was
        rejected.  Catches SPICE errors so a missing kernel / bad epoch
        never breaks the reject-recording path.
        """
        try:
            moon_pos = np.asarray(g.moon_position_j2000_km(et), dtype=float)
            dist_km = float(np.linalg.norm(moon_pos))
            if dist_km <= 0.0:
                raise ValueError("moon position has zero norm")
            dec_deg = float(np.degrees(np.arcsin(moon_pos[2] / dist_km)))
            ra_deg = float(np.degrees(np.arctan2(moon_pos[1], moon_pos[0])) % 360.0)
            return {
                "moon_dec_deg": dec_deg,
                "moon_ra_deg": ra_deg,
                "earth_moon_dist_km": dist_km,
            }
        except Exception as exc:  # pragma: no cover (defensive)
            _log.debug("moon ephemeris lookup failed at et=%s: %s", et, exc)
            return {
                "moon_dec_deg": float("nan"),
                "moon_ra_deg": float("nan"),
                "earth_moon_dist_km": float("nan"),
            }

    def _geom_extra(et: float, **extra: Any) -> dict[str, Any]:
        """Build an enriched ``extra`` dict combining moon ephemeris
        with the caller-supplied geometric context (inc_E, dtheta, ...).

        All keys are present in the dict even if NaN — the schema is
        stable across reject sites so downstream JSON consumers can
        rely on every key being available.
        """
        moon = _moon_geom(et)
        envelope = abs(site.lat_deg)
        out: dict[str, Any] = {
            "moon_dec_deg": moon["moon_dec_deg"],
            "moon_ra_deg": moon["moon_ra_deg"],
            "earth_moon_dist_km": moon["earth_moon_dist_km"],
            "inc_e_min_allowed": envelope,
            "inc_e_max_allowed": 180.0 - envelope,
            "dtheta_actual_deg": float(dtheta_deg),
            "i0_lo": float(i0_lo),
            "i0_hi": float(i0_hi),
        }
        out.update(extra)
        return out

    def _record_reject(reason: str, t_L: Optional[float], G: Optional[float],
                       aux: Optional[tuple], grid_et: Optional[float] = None,
                       **kwargs: Any) -> None:
        if rejects is None:
            return
        az_, uL_, inc_, om_, vp_ = (None, None, None, None, None) if aux is None else aux
        g_resid = 1e9 if (G is None or G >= SENTINEL_THRESHOLD) else G
        # If the rejection reason is tagged by a fail code from the inner
        # solver chain (e.g. "inc_below_site_lat", "gmst_fail", ...),
        # prefer that as the structured ``reason`` field and stash the
        # human-readable tag in ``extra`` for backward compatibility.
        fail_code = kwargs.pop("fail_code", None)
        effective_reason = fail_code if fail_code in _FAIL_CODES else reason
        # Geometry context for this grid point (only when we know the ET).
        et_for_extra = kwargs.pop("et_for_extra", grid_et)
        extra_kwargs = dict(kwargs)
        if et_for_extra is not None:
            extra_kwargs.setdefault(
                "inc_e_required_deg",
                float(inc_) if inc_ is not None else None,
            )
            extra_kwargs.setdefault("reason_tag", reason)
            enriched = _geom_extra(et_for_extra, **extra_kwargs)
        else:
            extra_kwargs.setdefault("reason_tag", reason)
            enriched = dict(extra_kwargs)
        rejects.append(_make_reject(
            te0_et=te0_et,
            T_h=T_s / 3600.0,
            family=family.name,
            branch=branch,
            reason=effective_reason,
            t_launch_et=t_L,
            inc_e_deg=inc_,
            raan_e_deg=None,
            omega_e_deg=om_,
            vp_kms=vp_,
            azimuth_deg=az_,
            u_l_deg=uL_,
            g_residual_s=g_resid,
            extra=enriched,
        ))

    for k, x in enumerate(xs):
        G_, t_L_, az_, uL_, inc_, om_, vp_, code_ = _G_te(
            x, T_s, family, site, branch, dtheta_deg, i0_lo, i0_hi, extend_revs,
            t_ascent_s=t_ascent_s, coast_max_s=coast_max_s,
            ascent_arc_deg=ascent_arc_deg,
        )
        G_vals[k] = G_
        t_L_vals[k] = t_L_
        aux_vals[k] = (az_, uL_, inc_, om_, vp_)
        fail_codes[k] = code_
        # Sentinel from _G_te: classify as inner/GMST/lambert failure.
        if t_L_ is None or (G_ is not None and G_ >= SENTINEL_THRESHOLD):
            _record_reject(
                "inner_solve_failed",
                t_L_,
                G_,
                aux_vals[k],
                grid_index=k,
                grid_et=float(x),
                fail_code=code_,
                inc_e_required_deg=(
                    float(inc_) if inc_ is not None else None
                ),
            )

    # Find adjacent valid pairs with opposite signs → refine with brentq.
    qualifying_roots: list[tuple[float, float]] = []  # (|Te_root − Te0|, Te_root)
    for k in range(n - 1):
        G1, t1 = G_vals[k], t_L_vals[k]
        G2, t2 = G_vals[k + 1], t_L_vals[k + 1]
        if t1 is None or t2 is None:
            continue
        if G1 * G2 >= 0.0:
            _record_reject(
                "no_sign_change",
                t1,
                G1,
                aux_vals[k],
                bracket_index=k,
                g1=G1,
                g2=G2,
                grid_et=float(xs[k]),
                inc_e_required_deg=(
                    float(aux_vals[k][2])
                    if aux_vals[k] and aux_vals[k][2] is not None
                    else None
                ),
            )
            continue
        a, b = float(xs[k]), float(xs[k + 1])
        try:
            r = float(_brentq(
                lambda x: _G_te(
                    x, T_s, family, site, branch, dtheta_deg, i0_lo, i0_hi,
                    extend_revs, t_ascent_s=t_ascent_s,
                    coast_max_s=coast_max_s,
                    ascent_arc_deg=ascent_arc_deg,
                )[0],
                a, b,
                xtol=_BRENT_XTOL_S,
            ))
        except (ValueError, RuntimeError):
            _record_reject(
                "brentq_fail",
                t1,
                G1,
                aux_vals[k],
                bracket_index=k,
                a=a,
                b=b,
                grid_et=float(xs[k]),
                inc_e_required_deg=(
                    float(aux_vals[k][2])
                    if aux_vals[k] and aux_vals[k][2] is not None
                    else None
                ),
            )
            continue
        G_r, t_L_r, az_r, uL_r, inc_r, om_r, vp_r, code_r = _G_te(
            r, T_s, family, site, branch, dtheta_deg, i0_lo, i0_hi, extend_revs,
            t_ascent_s=t_ascent_s, coast_max_s=coast_max_s,
            ascent_arc_deg=ascent_arc_deg,
        )
        # Discard discontinuity artifacts (GMST wrap): the apparent root
        # at the pocket boundary is huge in |G| and never survives this
        # tolerance test.
        if t_L_r is None or abs(G_r) > G_TOL_S:
            _record_reject(
                "discontinuity_artifact",
                t_L_r,
                G_r,
                (az_r, uL_r, inc_r, om_r, vp_r),
                r=r,
                G_r=G_r,
                grid_et=float(r),
                fail_code=code_r,
                inc_e_required_deg=(
                    float(inc_r) if inc_r is not None else None
                ),
            )
            continue
        win = _build_window(
            t_L=t_L_r,
            Te=r,
            T_h=T_s / 3600.0,
            family=family,
            branch=branch,
            inc_E_deg=inc_r,
            Omega_E_deg=om_r,
            az_deg=az_r,
            u_L_deg=uL_r,
            omega_E_deg=om_r,
            vp_kms=vp_r,
            polar_ref=polar_ref,
            H_m=H_m,
            t_ascent_s=t_ascent_s,
            rejects=rejects,
            te0_et=te0_et,
            site=site,
            dtheta_deg=dtheta_deg,
            i0_lo=i0_lo,
            i0_hi=i0_hi,
            ascent_arc_deg=ascent_arc_deg,
        )
        if win is not None:
            qualifying_roots.append((abs(r - Te0), win))

    if not qualifying_roots:
        return None

    # Pick the root closest to Te0.
    qualifying_roots.sort(key=lambda x: x[0])
    return qualifying_roots[0][1]


def _build_window(
    t_L: float,
    Te: float,
    T_h: float,
    family: pm.PerigeeFamily,
    branch: str,
    inc_E_deg: float,
    Omega_E_deg: float,
    az_deg: float,
    u_L_deg: float,
    omega_E_deg: float,
    vp_kms: float,
    polar_ref: str,
    H_m: float,
    t_ascent_s: float = T_ASCENT_S,
    rejects: Optional[list[WindowReject]] = None,
    te0_et: Optional[float] = None,
    *,
    site: Optional[LaunchSite] = None,
    dtheta_deg: float = 25.0,
    i0_lo: float = 19.6,
    i0_hi: float = 114.0,
    ascent_arc_deg: float = ASCENT_ARC_DEG,
) -> Optional[Window]:
    """Build the Window dataclass from the converged L6 geometry.

    ``omega_E_deg`` is here both the RAAN-axis label and the actual value;
    the function signature keeps the redundant name to make the dataflow
    obvious at the call site.

    When ``rejects`` is provided, drop points (coast-out-of-bounds,
    residual-too-large) append a ``WindowReject`` capturing the in-scope
    geometry before returning ``None``.

    ``site``, ``dtheta_deg``, ``i0_lo``, ``i0_hi`` are diagnostic-only
    fields used to enrich the reject record's ``extra`` dict with the
    moon ephemeris at ``Te`` and the geometric envelope bounds; they
    are keyword-only and optional so existing callers that don't pass
    them get the legacy (un-enriched) reject records.
    """
    coast_s = parking_orbit_coast_seconds(
        omega_E_deg, insertion_argument_deg(u_L_deg, ascent_arc_deg),
        extend_revs=True,
    )
    if coast_s < COAST_MIN_S or coast_s > COAST_MAX_S:
        # No number of revolutions brings the coast into bounds: the launch
        # opportunity is geometrically infeasible for this site / orbit.
        if rejects is not None:
            reject_extra: dict[str, Any] = {
                "coast_min_s": COAST_MIN_S,
                "coast_max_s": COAST_MAX_S,
            }
            # Enrich with the same geometry fields used in ``_pocket_search``
            # so downstream JSON consumers see a stable schema.
            if site is not None:
                _merge_moon_extra(reject_extra, float(Te), site, dtheta_deg,
                                  i0_lo, i0_hi)
                reject_extra["inc_e_required_deg"] = float(inc_E_deg)
            rejects.append(_make_reject(
                te0_et=float(te0_et if te0_et is not None else Te),
                T_h=T_h,
                family=family.name,
                branch=branch,
                reason="coast_out_of_bounds",
                t_launch_et=float(t_L),
                inc_e_deg=float(inc_E_deg),
                raan_e_deg=float(Omega_E_deg),
                omega_e_deg=float(omega_E_deg),
                vp_kms=float(vp_kms),
                azimuth_deg=float(az_deg) % 360.0,
                u_l_deg=float(u_L_deg),
                coast_s=float(coast_s),
                extra=reject_extra,
            ))
        return None

    t_TLI = t_L + t_ascent_s + coast_s
    # Sanity: the residual G should be small by construction; if it is not,
    # drop the window.
    residual = (t_TLI + T_h * 3600.0) - Te
    if abs(residual) > 60.0:
        if rejects is not None:
            r_park = ts.R_EARTH + H_PARK_KM
            c3 = vp_kms**2 - 2.0 * ts.MU_EARTH / r_park
            residual_extra: dict[str, Any] = {"t_TLI": float(t_TLI)}
            if site is not None:
                _merge_moon_extra(residual_extra, float(Te), site, dtheta_deg,
                                  i0_lo, i0_hi)
                residual_extra["inc_e_required_deg"] = float(inc_E_deg)
            rejects.append(_make_reject(
                te0_et=float(te0_et if te0_et is not None else Te),
                T_h=T_h,
                family=family.name,
                branch=branch,
                reason="residual_too_large",
                t_launch_et=float(t_L),
                inc_e_deg=float(inc_E_deg),
                raan_e_deg=float(Omega_E_deg),
                omega_e_deg=float(omega_E_deg),
                vp_kms=float(vp_kms),
                azimuth_deg=float(az_deg) % 360.0,
                u_l_deg=float(u_L_deg),
                coast_s=float(coast_s),
                c3_km2s2=float(c3),
                g_residual_s=float(residual),
                extra=residual_extra,
            ))
        return None

    # C3 (PRD v7 §六): C3 = v_p^2 − 2μ_E/(R_E + h_park).
    r_park = ts.R_EARTH + H_PARK_KM
    c3 = vp_kms**2 - 2.0 * ts.MU_EARTH / r_park

    return Window(
        t_launch_utc=g.et_to_utc(t_L),
        t_tli_utc=g.et_to_utc(t_TLI),
        t_arrive_utc=g.et_to_utc(Te),
        T_h=T_h,
        family=family.name,
        branch=branch,
        i_park_deg=float(inc_E_deg),
        raan_park_deg=float(Omega_E_deg),
        azimuth_deg=float(az_deg) % 360.0,
        coast_s=float(coast_s),
        vp_kms=float(vp_kms),
        c3_km2s2=float(c3),
        i_m_deg=float("nan"),   # populated by downstream re-solve
        perilune_alt_km=float(H_m),
        polar_ref=polar_ref,
        t_launch_et=float(t_L),
        t_arrive_et=float(Te),
    )


# ---------------------------------------------------------------------------
# Scan entry point
# ---------------------------------------------------------------------------
def scan_epochs(
    t_start_et: float,
    t_end_et: float,
    T_h_list: Sequence[float],
    families: dict[str, pm.PerigeeFamily],
    site: LaunchSite,
    polar_ref: str = "mep",
    branches: Sequence[str] = ("ascending", "descending"),
    dtheta_deg: float = 25.0,
    i0_lo: float = 19.6,
    i0_hi: float = 114.0,
    seed_step_s: float = SIDEREAL_DAY_S,
    return_rejects: bool = False,
    reject_limit: Optional[int] = None,
    t_ascent_s: float = T_ASCENT_S,
    coast_max_s: Optional[float] = None,
    ascent_arc_deg: float = ASCENT_ARC_DEG,
) -> list[Window] | tuple[list[Window], list[WindowReject]]:
    """Scan ``[t_start_et, t_end_et]`` for L6 launch windows.

    Seed ``T_e`` on a ``seed_step_s`` grid (default: one sidereal day),
    run :func:`solve_epoch` for every (T_h, family, branch) combination
    (explicit enumeration of the 4 L5 families × 2 升/降 branches per
    family), and dedup the resulting launches by ``t_launch_et`` within
    ``DEDUP_THRESHOLD_S`` (600 s).

    ``t_ascent_s`` is the launch→parking-orbit insertion time threaded to
    every :func:`solve_epoch` call (CE-5/CE-6: 740/738 s); the default is
    the legacy module constant ``T_ASCENT_S``.

    ``coast_max_s``, when given, is enforced *inside* the G evaluation
    (solver-level timeline constraint): grid points whose parking-orbit
    coast exceeds the cap are invalidated with the ``coast_exceeds_limit``
    fail code, so the pocket search only ever sees cap-compliant
    timelines.

    ``ascent_arc_deg`` is threaded to every :func:`solve_epoch` call: the
    parking-orbit coast starts at the insertion argument
    ``u_ins = u_L + ascent_arc_deg`` instead of the launch-site argument
    ``u_L`` (default :data:`ASCENT_ARC_DEG`; pass 0.0 for the legacy
    site-as-insertion convention).

    The seed grid step is decoupled from the pocket search radius: each
    ``solve_epoch`` call performs its own ±12 h / ±24 h pocket search
    around the seed, so seeds one sidereal day apart cannot produce
    duplicate windows — but the explicit dedup is retained as a safety
    net for callers that tighten ``seed_step_s`` below one sidereal day.
    """
    windows: list[Window] = []
    # ``rejects`` is mutated in-place by ``solve_epoch`` only when it is a
    # list; when ``return_rejects=False`` we keep it ``None`` so the
    # downstream ``if rejects is not None`` guards in ``solve_epoch`` /
    # ``_pocket_search`` / ``_make_window`` skip their construction work
    # entirely (zero overhead on the production path).
    rejects: list[WindowReject] | None = [] if return_rejects else None
    if seed_step_s <= 0.0:
        raise ValueError(f"seed_step_s must be positive, got {seed_step_s!r}")
    n = int((t_end_et - t_start_et) / seed_step_s) + 1
    Te_guesses = [t_start_et + i * seed_step_s for i in range(n)]

    for Te0 in Te_guesses:
        for T_h in T_h_list:
            for fam_name, family in families.items():
                for branch in branches:
                    win = solve_epoch(
                        Te_guess=Te0,
                        site=site,
                        T_h=T_h,
                        family=family,
                        polar_ref=polar_ref,
                        dtheta_deg=dtheta_deg,
                        i0_lo=i0_lo,
                        i0_hi=i0_hi,
                        branch=branch,
                        extend_revs=True,
                        rejects=rejects,
                        te0_et=Te0,
                        t_ascent_s=t_ascent_s,
                        coast_max_s=coast_max_s,
                        ascent_arc_deg=ascent_arc_deg,
                    )
                    if win is not None:
                        windows.append(win)

    # Dedup by t_launch_et (kept as the canonical ordering: ascending by
    # launch epoch, breaking ties by family/branch so the dedup window
    # keeps the first-encountered entry per launch).
    windows.sort(key=lambda w: (w.t_launch_et, w.family, w.branch))
    deduped: list[Window] = []
    for w in windows:
        if deduped and abs(w.t_launch_et - deduped[-1].t_launch_et) < DEDUP_THRESHOLD_S:
            # Record a "dedup" reject: the launch epoch was already
            # produced by an earlier (T_h, family, branch) combination and
            # therefore does not represent an independent opportunity.
            # ``te0_et`` is not strictly the seed arrival epoch here (the
            # dedup loop runs after the seed loop has ended), so we
            # substitute the actual recovered ``t_arrive_et`` from the
            # dropped ``Window`` and document the substitution in
            # ``extra``.  ``T_h``, ``family`` and ``branch`` come from the
            # dropped window itself (PRD v7 §一 fields on the Window
            # dataclass).
            if rejects is not None:
                dedup_extra: dict[str, Any] = {
                    "gap_s": w.t_launch_et - deduped[-1].t_launch_et,
                    "prev_t_launch_et": deduped[-1].t_launch_et,
                    "note": "te0_et substituted with arrival epoch",
                }
                # Enrich with moon ephemeris at the actual arrival epoch.
                _merge_moon_extra(
                    dedup_extra,
                    float(w.t_arrive_et),
                    site,
                    float(dtheta_deg),
                    float(i0_lo),
                    float(i0_hi),
                )
                dedup_extra["inc_e_required_deg"] = float(w.i_park_deg)
                rejects.append(_make_reject(
                    te0_et=w.t_arrive_et,
                    T_h=w.T_h,
                    family=w.family,
                    branch=w.branch,
                    reason="dedup",
                    t_launch_et=w.t_launch_et,
                    inc_e_deg=w.i_park_deg,
                    raan_e_deg=w.raan_park_deg,
                    vp_kms=w.vp_kms,
                    azimuth_deg=w.azimuth_deg,
                    coast_s=w.coast_s,
                    c3_km2s2=w.c3_km2s2,
                    extra=dedup_extra,
                ))
            continue
        deduped.append(w)

    if return_rejects:
        if reject_limit is not None and len(rejects) > reject_limit:
            rejects = rejects[:reject_limit]
        return deduped, rejects
    return deduped


# ---------------------------------------------------------------------------
# Rejection-construction helper (centralises WindowReject boilerplate).
# ---------------------------------------------------------------------------
def _make_reject(
    *,
    te0_et: float,
    T_h: float,
    family: str,
    branch: str,
    reason: str,
    t_launch_et: Optional[float] = None,
    inc_e_deg: Optional[float] = None,
    raan_e_deg: Optional[float] = None,
    omega_e_deg: Optional[float] = None,
    vp_kms: Optional[float] = None,
    azimuth_deg: Optional[float] = None,
    u_l_deg: Optional[float] = None,
    coast_s: Optional[float] = None,
    c3_km2s2: Optional[float] = None,
    g_residual_s: Optional[float] = None,
    i0_deg: Optional[float] = None,
    extra: Optional[dict[str, Any]] = None,
) -> WindowReject:
    """Construct a :class:`WindowReject` from in-scope locals.

    Centralises construction so each drop point is a single line.  Pass
    only the fields that are meaningful at the call site; everything else
    stays ``None``.
    """
    return WindowReject(
        te0_et=float(te0_et),
        t_h_hours=float(T_h),
        family=family,
        branch=branch,
        reason=reason,
        t_launch_et=(float(t_launch_et) if t_launch_et is not None else None),
        inc_e_deg=(float(inc_e_deg) if inc_e_deg is not None else None),
        raan_e_deg=(float(raan_e_deg) if raan_e_deg is not None else None),
        omega_e_deg=(float(omega_e_deg) if omega_e_deg is not None else None),
        vp_kms=(float(vp_kms) if vp_kms is not None else None),
        azimuth_deg=(float(azimuth_deg) if azimuth_deg is not None else None),
        u_l_deg=(float(u_l_deg) if u_l_deg is not None else None),
        coast_s=(float(coast_s) if coast_s is not None else None),
        c3_km2s2=(float(c3_km2s2) if c3_km2s2 is not None else None),
        g_residual_s=(float(g_residual_s) if g_residual_s is not None else None),
        i0_deg=(float(i0_deg) if i0_deg is not None else None),
        extra=(dict(extra) if extra is not None else {}),
    )


def _merge_moon_extra(
    target: dict[str, Any],
    et: float,
    site: LaunchSite,
    dtheta_deg: float,
    i0_lo: float,
    i0_hi: float,
) -> None:
    """In-place: merge moon ephemeris + envelope bounds into ``target``.

    Adds the keys ``moon_dec_deg``, ``moon_ra_deg``, ``earth_moon_dist_km``,
    ``inc_e_min_allowed``, ``inc_e_max_allowed``, ``dtheta_actual_deg``,
    ``i0_lo``, ``i0_hi`` to ``target`` (any pre-existing keys with the
    same name are overwritten — callers should not pre-set them).
    Catches SPICE errors so a missing kernel / bad epoch never breaks
    the reject-recording path; in that case the keys are present with
    NaN values to keep the schema stable.
    """
    try:
        moon_pos = np.asarray(g.moon_position_j2000_km(et), dtype=float)
        dist_km = float(np.linalg.norm(moon_pos))
        if dist_km <= 0.0:
            raise ValueError("moon position has zero norm")
        dec_deg = float(np.degrees(np.arcsin(moon_pos[2] / dist_km)))
        ra_deg = float(np.degrees(np.arctan2(moon_pos[1], moon_pos[0])) % 360.0)
    except Exception as exc:  # pragma: no cover (defensive)
        _log.debug("moon ephemeris lookup failed at et=%s: %s", et, exc)
        dec_deg = float("nan")
        ra_deg = float("nan")
        dist_km = float("nan")
    envelope = abs(site.lat_deg)
    target["moon_dec_deg"] = dec_deg
    target["moon_ra_deg"] = ra_deg
    target["earth_moon_dist_km"] = dist_km
    target["inc_e_min_allowed"] = envelope
    target["inc_e_max_allowed"] = 180.0 - envelope
    target["dtheta_actual_deg"] = float(dtheta_deg)
    target["i0_lo"] = float(i0_lo)
    target["i0_hi"] = float(i0_hi)
