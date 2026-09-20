"""launch_window/launcher.py — L7 capability ladder & Δv-budget (PRD v7 §六, §七).

Provides:

- :class:`Launcher`         — dataclass describing a launch vehicle config.
- :func:`cz5_capability_ladder` — 13-tap CZ-5 capability ladder.
- :func:`c3_curve_linear`  — analytical payload mass curve, rocket-equation linearised.
- :func:`feasible`         — feasibility check at a single window.
- :func:`apply_capability_ladder` — batch feasibility over a capability ladder.
- :func:`width_from_margin`— window opening width from Δv-yaw margin.
- :func:`annotate_width`   — batch width annotation.
- :class:`DVBudget`        — Δv budget breakdown (PRD v7 §七).
- :func:`dv_budget`        — compute Δv breakdown for a single window dict.

Physical constants (km, km/s, kg, s) are imported from
``launch_window.tli_solver`` (single source of truth).  Only one new constant
is defined locally — the Earth sidereal angular velocity ω_E (which the PRD
v7 §六 requires for the Δv-yaw inversion but which is not present in
``tli_solver.py``).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Iterable, Optional, Sequence

from launch_window import tli_solver

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
# Earth sidereal angular velocity (rad/s).  PRD v7 §六 uses this in the
# Δv-yaw inversion.  Not defined in tli_solver.py, hence defined here once.
OMEGA_EARTH = 7.2921159e-5  # rad/s  (Earth's sidereal rotation rate)

# Re-export for ergonomics so callers can ``from launch_window.launcher import MU_EARTH``.
MU_EARTH = tli_solver.MU_EARTH  # km^3 / s^2
MU_MOON = tli_solver.MU_MOON
R_EARTH = tli_solver.R_EARTH
R_MOON = tli_solver.R_MOON

G0 = 9.80665e-3  # km/s^2  (standard gravity in km/s^2 units)

# Default YF-75D vacuum specific impulse for the CZ-5 upper stage (s).
# Documented in PRD v7 §六 as the linearisation reference.
ISP_YF75D = 442.0  # s

# Reference C3 (km^2/s^2) and TLI reference payload (kg) for the 13-tap ladder.
# PRD v7 §六: ladder spans 8200 … 9400 kg in 100 kg steps at C3_ref = −2.0 km^2/s^2.
C3_REF = -2.0  # km^2/s^2
LADDER_MIN_KG = 8200.0
LADDER_MAX_KG = 9400.0
LADDER_STEP_KG = 100.0

DEFAULT_H_PARK_KM = 200.0
DEFAULT_H_PARK_PERI_KM = 200.0
DEFAULT_H_PARK_APO_KM = 380.0

DEFAULT_COAST_MIN_S = 600.0
DEFAULT_COAST_MAX_S = 6000.0

# Nominal coast (s) used for seeding / window-width bookkeeping.
DEFAULT_COAST_NOMINAL_S = 1000.0

# Reference perilune altitude (km) used by dv_budget for v_circ_moon.
DEFAULT_H_PERILUNE_KM = 300.0

# Moon equator vs. MEP plane offset (deg).  PRD v7 §七 dv_trim estimate.
MOON_EQUATOR_MEP_OFFSET_DEG = 6.7


# ---------------------------------------------------------------------------
# Launcher dataclass
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Launcher:
    """Capability description of a single launch-vehicle configuration.

    Parameters
    ----------
    name : str
        Identifier (e.g. ``"CZ-5/YZ-2 @ 8200 kg"``).
    m_payload_kg : float
        Probe wet mass to be delivered (kg).
    c3_curve : Callable[[float], float]
        Maps a C3 (km^2/s^2) to the *available* payload mass (kg).  See
        :func:`c3_curve_linear` for a rocket-equation linearisation that
        passes the L7 monotonicity test.
    dv_margin_mps : float, default 30.0
        Δv (m/s) the vehicle can spend on yaw steering / dogleg before TLI.
        Used by :func:`width_from_margin` to bound the launch window.
    coast_min_s, coast_max_s : float
        Allowed parking-orbit coast duration bounds (seconds).  Corrected
        default 600 s ≤ coast ≤ 6000 s (attitude-stabilisation minimum /
        battery-thermal limit; CE-5/CE-6 single-orbit coast ≈ 937/1020 s).
    coast_nominal_s : float, default 1000.0
        Nominal (measured) parking-orbit coast used for seeding the outer
        launch-epoch search and window-width bookkeeping.
    h_park_km : float, default 200.0
        Parking-orbit perigee altitude (km); legacy alias of
        ``h_park_peri_km``.
    h_park_peri_km, h_park_apo_km : float, default 200.0 / 380.0
        Parking-orbit perigee / apogee altitudes (km) of the elliptical
        CZ-5 parking orbit.  TLI is performed at perigee.
    """

    name: str
    m_payload_kg: float
    c3_curve: Callable[[float], float]
    dv_margin_mps: float = 30.0
    coast_min_s: float = DEFAULT_COAST_MIN_S
    coast_max_s: float = DEFAULT_COAST_MAX_S
    coast_nominal_s: float = DEFAULT_COAST_NOMINAL_S
    h_park_km: float = DEFAULT_H_PARK_KM
    h_park_peri_km: float = DEFAULT_H_PARK_PERI_KM
    h_park_apo_km: float = DEFAULT_H_PARK_APO_KM


# ---------------------------------------------------------------------------
# C3 → payload capability (rocket-equation linearisation)
# ---------------------------------------------------------------------------
def _dv_for_c3(c3_kms2: float, mu_over_r: float) -> float:
    """Δv (km/s) to escape from a circular parking orbit to hyperbolic C3.

    ``Δv(C3) = sqrt(C3 + 2 μ_E / r_park) - sqrt(μ_E / r_park)``.

    Equivalently ``Δv = sqrt(C3 + 2 μ/r) - sqrt(μ/r)`` with the second
    term being the circular-park velocity.  We accept ``mu_over_r`` so
    callers can amortise the constant if they call this in a tight loop.
    """
    v_p = math.sqrt(c3_kms2 + 2.0 * mu_over_r)
    v_circ = math.sqrt(mu_over_r)
    return v_p - v_circ


def c3_curve_linear(
    m_ref_kg: float,
    c3_ref_kms2: float = C3_REF,
    isp_s: float = ISP_YF75D,
    h_park_km: float = DEFAULT_H_PARK_KM,
    mu_earth: float = MU_EARTH,
    g0_kms2: float = G0,
) -> Callable[[float], float]:
    """Return a monotonically-decreasing ``C3 → m_avail`` callable.

    The payload curve is obtained by linearising the Tsiolkovsky rocket
    equation around a reference C3 reference point::

        Δv(C3) = sqrt(C3 + 2 μ_E / r_park) − sqrt(μ_E / r_park)
        m_avail(C3) = m_ref · exp(−(Δv(C3) − Δv(C3_ref)) / (Isp · g0))

    Interpretation
    --------------
    ``Δv(C3) − Δv(C3_ref)`` is the *incremental* hyperbolic injection Δv
    relative to the reference.  Multiplying by the standard Tsiolkovsky
    exponent ``1 / (Isp · g0)`` converts the Δv into a logarithmic mass
    penalty.  This is the standard rocket-equation linearisation used by
    launch-vehicle performance handbooks for small C3 swings (the PRD v7
    notes that ``v_p`` differs by only a few m/s across the launch-window
    family, so a linear-in-Δv model is well within 0.1 % of truth).

    Parameters
    ----------
    m_ref_kg : float
        Available payload (kg) at the reference C3.
    c3_ref_kms2 : float, default -2.0
        Reference C3 (km^2/s^2) at which ``m_avail(c3_ref) = m_ref_kg``.
    isp_s : float, default 442.0
        YF-75D vacuum specific impulse (s).
    h_park_km : float, default 200.0
        Parking-orbit altitude (km).
    mu_earth : float, default 398600.4418
        Earth GM (km^3/s^2) from :mod:`tli_solver`.
    g0_kms2 : float, default 9.80665e-3
        Standard gravity (km/s^2).

    Returns
    -------
    Callable[[float], float]
        ``m_avail(C3)`` in kg.
    """
    r_park = R_EARTH + h_park_km
    mu_over_r = mu_earth / r_park
    dv_ref = _dv_for_c3(c3_ref_kms2, mu_over_r)
    inv_isp_g0 = 1.0 / (isp_s * g0_kms2)

    def m_avail(c3_kms2: float) -> float:
        dv = _dv_for_c3(c3_kms2, mu_over_r)
        return m_ref_kg * math.exp(-(dv - dv_ref) * inv_isp_g0)

    return m_avail


def cz5_capability_ladder(
    ladder_min_kg: float = LADDER_MIN_KG,
    ladder_max_kg: float = LADDER_MAX_KG,
    ladder_step_kg: float = LADDER_STEP_KG,
    c3_ref_kms2: float = C3_REF,
    isp_s: float = ISP_YF75D,
) -> list[Launcher]:
    """Return the 13-tap CZ-5 capability ladder (PRD v7 §六).

    The ladder indexes capability by the available-payload mass at the
    reference C3 ``C3_REF = -2.0 km^2/s^2``::

        m_ref ∈ {8200, 8300, …, 9400} kg

    Each tap is a :class:`Launcher` whose ``c3_curve`` is the
    rocket-equation-linearised curve anchored at ``(C3_ref, m_ref)``.
    """
    n_taps = int(round((ladder_max_kg - ladder_min_kg) / ladder_step_kg)) + 1
    ladder: list[Launcher] = []
    for k in range(n_taps):
        m_ref = ladder_min_kg + k * ladder_step_kg
        c3_curve = c3_curve_linear(
            m_ref_kg=m_ref, c3_ref_kms2=c3_ref_kms2, isp_s=isp_s
        )
        ladder.append(
            Launcher(
                name=f"CZ-5/YZ-2 @{int(round(m_ref))}kg",
                m_payload_kg=m_ref,
                c3_curve=c3_curve,
            )
        )
    return ladder


# ---------------------------------------------------------------------------
# Feasibility per window
# ---------------------------------------------------------------------------
def feasible(launcher: Launcher, window_dict: dict) -> tuple[bool, float]:
    """Feasibility check at a single window.

    Parameters
    ----------
    launcher : Launcher
        Launch-vehicle configuration.
    window_dict : dict
        Must contain ``"c3"`` (km^2/s^2).  Optionally ``"margin_kg"``
        (will be overwritten).

    Returns
    -------
    (ok, margin_kg) : (bool, float)
        ``ok = m_avail(C3) >= m_payload_kg``;
        ``margin_kg = m_avail(C3) - m_payload_kg``.
    """
    c3 = float(window_dict["c3"])
    m_avail = float(launcher.c3_curve(c3))
    margin = m_avail - launcher.m_payload_kg
    window_dict["margin_kg"] = margin
    window_dict["m_avail_kg"] = m_avail
    return margin >= 0.0, margin


def apply_capability_ladder(
    windows: Sequence[dict], ladder: Sequence[Launcher]
) -> dict:
    """Apply the full capability ladder to a list of windows.

    For each window ``w`` and each tier ``L``::

        w[f"margin_kg_{int(L.m_payload_kg)}"]   = m_avail(w.c3) - L.m_payload_kg
        w[f"feasible_{int(L.m_payload_kg)}"]    = (margin >= 0)

    A consolidated field ``w["margin_kg"]`` (margin relative to the
    *lightest* tier — i.e. the smallest ``m_payload_kg`` in the ladder,
    conventionally 8200 kg) is also added, matching the PRD dataclass.

    Returns
    -------
    summary : dict
        ``{"per_window": [...], "feasible_count": {tier_kg: n, ...},
            "total_windows": int}``.
    """
    if not ladder:
        raise ValueError("ladder is empty")
    base_kg = int(min(L.m_payload_kg for L in ladder))

    feasible_count: dict[int, int] = {int(L.m_payload_kg): 0 for L in ladder}
    for w in windows:
        c3 = float(w["c3"])
        # m_avail_kg at the lightest (base) tier is the canonical reference.
        m_avail_base = float(ladder[0].c3_curve(c3))
        w["m_avail_kg"] = m_avail_base
        for L in ladder:
            m_avail = float(L.c3_curve(c3))
            margin = m_avail - L.m_payload_kg
            tier_kg = int(L.m_payload_kg)
            w[f"margin_kg_{tier_kg}"] = margin
            w[f"feasible_{tier_kg}"] = margin >= 0.0
            if margin >= 0.0:
                feasible_count[tier_kg] += 1
        # PRD-aligned base margin (smallest tier):
        w["margin_kg"] = w.get(f"margin_kg_{base_kg}", float("nan"))

    return {
        "per_window": list(windows),
        "feasible_count": feasible_count,
        "total_windows": len(windows),
    }


# ---------------------------------------------------------------------------
# Window width from Δv-yaw margin (PRD v7 §六)
# ---------------------------------------------------------------------------
def width_from_margin(
    i_park_deg: float,
    site_lat_deg: float,
    dv_margin_mps: float,
    v_park_kms: float,
    omega_earth: float = OMEGA_EARTH,
) -> dict:
    """Compute launch-window opening width from a Δv-yaw budget.

    Model (PRD v7 §六)::

        Δv_yaw ≈ 2 v_park · sin(δΨ / 2)
        δΨ    ≈ ω_E · cos(φ_L) · sin(i_E) · δt

    Solving ``Δv_yaw = dv_margin`` for ``δt``::

        δt_max = 2 · asin(dv_margin / (2 v_park)) / (ω_E · cos(φ_L) · sin(i_E))

    The full window opening is ``width_min = 2 · δt_max`` (symmetric about
    the geometric optimum).

    Parameters
    ----------
    i_park_deg : float
        Parking-orbit inclination relative to the Earth equator (deg).
    site_lat_deg : float
        Launch-site geodetic latitude (deg, |φ_L| ≤ 90).
    dv_margin_mps : float
        Δv budget for yaw steering (m/s).
    v_park_kms : float
        Circular parking-orbit speed (km/s).

    Returns
    -------
    dict with keys ``delta_t_max_s``, ``width_min``, ``delta_psi_max_deg``,
    ``feasible`` (bool).  When ``|φ_L| > i_E`` the launch site cannot reach
    the orbit plane at all, so the width is 0.
    """
    i_rad = math.radians(i_park_deg)
    phi_rad = math.radians(site_lat_deg)
    # Convert dv margin (m/s) → (km/s) for unit consistency with v_park.
    dv_margin_kms = dv_margin_mps * 1.0e-3
    if v_park_kms <= 0.0:
        raise ValueError("v_park_kms must be positive")
    if dv_margin_kms <= 0.0:
        return {
            "delta_t_max_s": 0.0,
            "width_min": 0.0,
            "delta_psi_max_deg": 0.0,
            "feasible": True,
        }
    # Geometric accessibility (PRD v7 §五 constraint A):
    if abs(math.sin(phi_rad)) > abs(math.sin(i_rad)) + 1e-12:
        return {
            "delta_t_max_s": 0.0,
            "width_min": 0.0,
            "delta_psi_max_deg": 0.0,
            "feasible": False,
        }
    sin_i = math.sin(i_rad)
    cos_phi = math.cos(phi_rad)
    if sin_i * cos_phi <= 0.0:
        return {
            "delta_t_max_s": 0.0,
            "width_min": 0.0,
            "delta_psi_max_deg": 0.0,
            "feasible": False,
        }
    # δΨ from Δv_yaw = 2 v_park sin(δΨ/2):
    arg = dv_margin_kms / (2.0 * v_park_kms)
    arg = max(-1.0, min(1.0, arg))
    delta_psi = 2.0 * math.asin(arg)
    delta_t_max = delta_psi / (omega_earth * cos_phi * sin_i)
    return {
        "delta_t_max_s": delta_t_max,
        "width_min": 2.0 * delta_t_max,
        "delta_psi_max_deg": math.degrees(delta_psi),
        "feasible": True,
    }


def annotate_width(
    windows: Sequence[dict],
    site_lat_deg: float,
    dv_margin_mps: float,
    omega_earth: float = OMEGA_EARTH,
) -> Sequence[dict]:
    """Annotate each window dict in-place with ``width_min`` etc.

    Reads from each window dict:
        - ``"i_park_deg"``      parking-orbit inclination (deg)
        - ``"v_park_kms"``      parking-orbit speed (km/s); defaults to
                                7.79 km/s (200-km LEO) if missing.

    Writes:
        - ``"width_min"``       full opening width (seconds)
        - ``"delta_t_max_s"``   half-width (seconds)
        - ``"delta_psi_max_deg"`` (deg)
        - ``"width_feasible"``  (bool) — synonym of ``"feasible"``.
    """
    for w in windows:
        i_park_deg = float(w.get("i_park_deg", 0.0))
        v_park_kms = float(w.get("v_park_kms", 7.79))
        result = width_from_margin(
            i_park_deg=i_park_deg,
            site_lat_deg=site_lat_deg,
            dv_margin_mps=dv_margin_mps,
            v_park_kms=v_park_kms,
            omega_earth=omega_earth,
        )
        w["width_min"] = result["width_min"]
        w["delta_t_max_s"] = result["delta_t_max_s"]
        w["delta_psi_max_deg"] = result["delta_psi_max_deg"]
        w["width_feasible"] = result["feasible"]
    return windows


# ---------------------------------------------------------------------------
# Δv budget (PRD v7 §七)
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class DVBudget:
    """Δv budget breakdown (PRD v7 §七).  All values in km/s.

    Attributes
    ----------
    dv_tli    : TLI burn (parking circular → hyperbolic).
    dv_mcc    : Mid-course correction budget (default 0; caller may override).
    dv_loi    : Lunar-orbit insertion burn (hyperbolic → circular at 200 km).
    dv_trim   : Plane-change trim (0 when ``polar_ref='moon_equator'``,
                ≈ 2·v·sin(3.35°) when ``polar_ref='mep'``).
    dv_total  : Sum of the four components.
    v_park    : Circular parking speed used (km/s).
    v_circ_moon : Circular lunar-orbit speed at the perilune altitude (km/s).
    notes     : Free-form annotations (``None`` when not set).
    """

    dv_tli: Optional[float]
    dv_mcc: Optional[float]
    dv_loi: Optional[float]
    dv_trim: Optional[float]
    dv_total: Optional[float]
    v_park: Optional[float] = None
    v_circ_moon: Optional[float] = None
    notes: Optional[str] = None


def v_park_kms(h_park_km: float = DEFAULT_H_PARK_KM, mu: float = MU_EARTH) -> float:
    """Circular parking-orbit speed (km/s) at altitude ``h_park_km``."""
    r = R_EARTH + h_park_km
    return math.sqrt(mu / r)


def v_circ_moon_kms(
    h_perilune_km: float = DEFAULT_H_PERILUNE_KM, mu: float = MU_MOON
) -> float:
    """Circular lunar-orbit speed (km/s) at ``h_perilune_km`` altitude."""
    r = R_MOON + h_perilune_km
    return math.sqrt(mu / r)


def dv_budget(
    window_dict: dict,
    polar_ref: str = "moon_equator",
    dv_mcc: float = 0.0,
    h_park_km: float = DEFAULT_H_PARK_KM,
    h_perilune_km: float = DEFAULT_H_PERILUNE_KM,
    v_perilune_kms: Optional[float] = None,
    vp_kms: Optional[float] = None,
    c3_kms2: Optional[float] = None,
) -> DVBudget:
    """Build the Δv budget for a single window (PRD v7 §七).

    Parameters
    ----------
    window_dict : dict
        Must contain enough info to compute dv_tli, dv_loi and dv_trim.
        Two equivalent input conventions are supported:

        - Direct:  ``window_dict["vp_kms"]`` (or ``"v_park_kms"``) for
          dv_tli and ``window_dict["v_perilune_kms"]`` for dv_loi/dv_trim.
        - Indirect: ``window_dict["c3"]`` (km²/s²) is used to back out
          ``vp_kms`` (the TLI perigee speed) and ``v_perilune_kms`` is
          taken from the explicit argument or, if absent, left ``None``.

    polar_ref : {'mep', 'moon_equator'}
        ``'mep'`` includes the ≈6.7° Moon-equator plane change
        (``dv_trim = 2 v_perilune sin(3.35°)``); ``'moon_equator'`` sets
        ``dv_trim = 0`` because the q3 constraint already targets the
        true lunar polar axis.
    dv_mcc : float, default 0.0
        Mid-course correction budget (km/s).  Kept at 0 here because
        the PRD leaves it as an explicit parameter (caller can override).
    h_park_km, h_perilune_km : float
        Reference orbit altitudes (km).
    v_perilune_kms, vp_kms, c3_kms2 : float, optional
        Direct overrides for the three underlying scalar inputs.

    Notes
    -----
    When ``v_perilune_kms`` is not available from the window dict or
    explicit argument, ``dv_loi`` and ``dv_trim`` are returned as
    ``None`` and ``dv_total`` is also ``None``.  Callers (typically the
    R-1 solver) are expected to populate these values downstream.
    """
    # ----- dv_tli ----------------------------------------------------------
    if vp_kms is None:
        vp_kms = window_dict.get("vp_kms")
    if vp_kms is None:
        vp_kms = window_dict.get("v_park_kms")
    if vp_kms is None and c3_kms2 is None:
        c3_kms2 = window_dict.get("c3")
    if vp_kms is None and c3_kms2 is not None:
        # C3 = vp^2 - 2 μ_E / r_park → vp = sqrt(C3 + 2 μ_E / r_park)
        r_park = R_EARTH + h_park_km
        vp_kms = math.sqrt(float(c3_kms2) + 2.0 * MU_EARTH / r_park)
    if vp_kms is None:
        dv_tli: Optional[float] = None
        v_park_used: Optional[float] = None
    else:
        v_park_used = v_park_kms(h_park_km)
        dv_tli = float(vp_kms) - v_park_used

    # ----- dv_loi & dv_trim (need v_perilune) ------------------------------
    if v_perilune_kms is None:
        v_perilune_kms = window_dict.get("v_perilune_kms")
    v_circ_moon = v_circ_moon_kms(h_perilune_km)
    if v_perilune_kms is None:
        dv_loi: Optional[float] = None
        dv_trim: Optional[float] = None
    else:
        dv_loi = float(v_perilune_kms) - v_circ_moon
        if polar_ref == "mep":
            dv_trim = 2.0 * float(v_perilune_kms) * math.sin(
                math.radians(MOON_EQUATOR_MEP_OFFSET_DEG / 2.0)
            )
        else:
            dv_trim = 0.0

    # ----- total -----------------------------------------------------------
    if dv_tli is None or dv_loi is None or dv_trim is None:
        dv_total: Optional[float] = None
    else:
        dv_total = float(dv_tli) + float(dv_mcc) + float(dv_loi) + float(dv_trim)

    notes = None
    if dv_loi is None:
        notes = (
            "dv_loi / dv_trim are None because v_perilune_kms was not "
            "provided.  Callers (R-1 solver) must populate this key."
        )

    return DVBudget(
        dv_tli=dv_tli,
        dv_mcc=dv_mcc,
        dv_loi=dv_loi,
        dv_trim=dv_trim,
        dv_total=dv_total,
        v_park=v_park_used,
        v_circ_moon=v_circ_moon,
        notes=notes,
    )
