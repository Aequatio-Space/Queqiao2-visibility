"""launch_window/tli_solver.py — R-1 lunar-side boundary value solver.

Implements the lunar-side boundary value problem of the PRD v4 §R-1.

Variables (paper eq. 4.1-4.2)
-----------------------------
``p = (v_p, Ω, ω)``  perigee velocity, ascending-node longitude, argument of
perigee (all in the MEP frame).

``q = (q1, q2, q3)``
    q1 = H_m                       perilune altitude (km)
    q2 = cos f = (r1·v1)/(r1 v1)   Moon-centred position-velocity angle cosine
    q3 = cos θ = k̂·(r1×v1)/|r1×v1| Moon-centred orbital-plane normal vs MEP Z

``r1, v1`` are the Moon-centred position and velocity at arrival.

Given parameters: perigee altitude ``h_p``, transfer time ``T``, initial orbit
inclination ``i`` (relative to the MEP fundamental plane), arrival epoch
``t_arrive``.

Targets
-------
- Coplanar (validation):  q_target = (H_m, 0, ±1)   (+1 prograde, −1 retrograde)
- Polar lunar (CE-7):     q_target = (H_m, 0, 0)    (q3 = 0 ⇒ i_m = 90°)

Dynamics
--------
Restricted four-body problem (Earth + Moon + Sun point masses) in the
Earth-centred J2000 inertial frame.  The equation of motion includes:

- Earth direct gravity
- Moon direct gravity + Moon indirect term (third-body correction for
  using an Earth-centred inertial frame)
- Solar third-body perturbation (opt-in via ``use_sun=True``; validated
  against validate_mode_c.py:277-295)

Integrator DOP853, rtol=1e-12, atol=1e-9.

The Newton iteration uses central finite differences for the Jacobian
(Δv_p = 1e-4 km/s, ΔΩ = Δω = 1e-3 deg), q1 normalised by H_m/1000 km, and
``scipy.optimize.root(method='hybr')``.  No state-transition-matrix
integration is implemented (the PRD explicitly defers it as an efficiency
optimisation).

Outputs (all returned, downstream steps depend on them)
-------------------------------------------------------
perigee_lon_mep, perigee_lat_mep       ← R-2 input
perilune_lon_mep, perilune_lat_mep     ← R-3/R-6 input
v_perilune, i_m, C3
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal, Optional

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

from launch_window import mep_frame as mf

# Polar reference axis selector (PRD v7 §四.1 / §一 Mission).
# 'mep' (default, backward-compatible) uses the MEP Z axis (white-pole
# approximation; the paper's convention that ignores the 6.7° Moon-equator
# tilt).  'moon_equator' uses the IAU_MOON Z axis (true lunar polar axis)
# expressed in J2000 at the arrival epoch.  This is the engineering switch
# that lets the solver target the actual polar-lunar orbit without a
# post-solve plane correction.
PolarRef = Literal["mep", "moon_equator"]
_VALID_POLAR_REFS = ("mep", "moon_equator")


def _polar_ref_axis(polar_ref: str, t_arrive: float, frame: mf.MEPFrame) -> np.ndarray:
    """Return the unit vector against which q3 is computed.

    Parameters
    ----------
    polar_ref : {'mep', 'moon_equator'}
        ``'mep'`` returns ``frame.z`` (the frozen white-pole direction).
        ``'moon_equator'`` returns the IAU_MOON Z axis at ``t_arrive``
        expressed in J2000 (the true lunar polar direction), via SPICE's
        ``pxform('IAU_MOON', 'J2000', et)``.
    t_arrive : float
        SPICE ephemeris time at arrival (only used for ``'moon_equator'``).
    frame : MEPFrame
        The frozen MEP frame (used to return ``frame.z`` for ``'mep'``).

    Returns
    -------
    (3,) ndarray
        Unit polar-reference axis in J2000.
    """
    if polar_ref == "mep":
        return frame.z
    if polar_ref == "moon_equator":
        import spiceypy
        from launch_window import geometry as g
        g.ensure_kernels()
        rot = spiceypy.pxform("IAU_MOON", "J2000", float(t_arrive))
        # rot is 3x3 whose columns are IAU_MOON axes in J2000; the 3rd column
        # is the IAU_MOON Z axis (true polar axis) in J2000.
        return np.asarray(rot[:, 2], dtype=float)
    raise ValueError(
        f"polar_ref must be one of {_VALID_POLAR_REFS}, got {polar_ref!r}"
    )

# ---------------------------------------------------------------------------
# Physical constants (km, km/s, km^3/s^2)
# ---------------------------------------------------------------------------
MU_EARTH = 398600.4418
MU_MOON = 4902.800066
MU_SUN = 1.32712440018e11  # km³/s², DE440
R_EARTH = 6378.137
R_MOON = 1737.4

# Integration settings (PRD R-1)
RTOL = 1e-12
ATOL = 1e-9

# Finite-difference steps (PRD R-1)
DV_P = 1e-4  # km/s
D_ANGLE = 1e-3  # deg

# Default perilune altitude (km).
# B-3 nominal is CE-6's 200 km / 3 burns, but the R-1V golden regression
# against the paper's Table 1 calibrates to H_m = 300 km (the perilune
# longitude matches to < 0.1° with this value).
DEFAULT_H_M = 300.0

# q1 normalisation scale (PRD R-1: q1 normalised by H_m/1000 km)
Q1_SCALE_KM = 1000.0

# Circular-orbit model parameters (paper Table 1 assumption)
# The paper's Table 1 assumes a circular lunar orbit.  The PRD instructs that
# if the SPICE-truth comparison exceeds tolerance, re-run with the mean
# Earth-Moon distance 384,400 km before doubting the solver.
#
# NOTE (R-1V finding): the paper's ~12°/day figure is the Moon's *synodic*
# motion (relative to the Sun), NOT the sidereal mean motion on its orbit.
# Using 12°/day for the circular-orbit model degrades the perilune-longitude
# match from <0.5° to 2-6°, so the sidereal mean motion (13.18°/day,
# 27.32-day period) is used here.
MOON_MEAN_DISTANCE_KM = 384400.0
MOON_SIDEREAL_PERIOD_S = 27.321661 * 86400.0  # sidereal month
MOON_MEAN_MOTION = 2.0 * math.pi / MOON_SIDEREAL_PERIOD_S  # rad/s


# ---------------------------------------------------------------------------
# Moon state helper (SPICE truth, non-circular orbit)
# ---------------------------------------------------------------------------
def _moon_state_j2000(t: float) -> tuple[np.ndarray, np.ndarray]:
    """Return (r_moon, v_moon) in J2000 (km, km/s) at epoch ``t`` (ET)."""
    import spiceypy

    state, _lt = spiceypy.spkezr("MOON", t, "J2000", "NONE", "EARTH")
    return np.asarray(state[:3], dtype=float), np.asarray(state[3:6], dtype=float)


def _moon_state_circular(t: float, t_arrive: float, frame: mf.MEPFrame) -> tuple[np.ndarray, np.ndarray]:
    """Return (r_moon, v_moon) in J2000 for a circular lunar orbit.

    The Moon moves on a circle of radius ``MOON_MEAN_DISTANCE_KM`` in the MEP
    fundamental plane.  In the frozen MEP frame the Moon is at longitude
    ``π + ω_m·(t − t_arrive)`` (it is at 180° at ``t_arrive``, the anti-X
    direction), moving counterclockwise (prograde) with mean motion ``ω_m``.
    The state is then rotated to J2000 via the frozen MEP frame.
    """
    lon = math.pi + MOON_MEAN_MOTION * (t - t_arrive)
    r_mep = MOON_MEAN_DISTANCE_KM * np.array([math.cos(lon), math.sin(lon), 0.0])
    v_mep = MOON_MEAN_DISTANCE_KM * MOON_MEAN_MOTION * np.array(
        [-math.sin(lon), math.cos(lon), 0.0]
    )
    return mf.from_mep(r_mep, frame), mf.from_mep(v_mep, frame)


def _rhs(
    t: float,
    y: np.ndarray,
    t_arrive: float,
    frame: mf.MEPFrame,
    moon_model: str,
    use_sun: bool = False,
) -> np.ndarray:
    """Restricted four-body acceleration in J2000 (Earth+Moon+Sun point masses).

    The Earth-centred equation of motion includes:

    - Earth direct gravity:  -μ_E · r / |r|³
    - Moon direct gravity:   -μ_M · (r − r_moon) / |r − r_moon|³
    - Moon indirect term:    -μ_M · r_moon / |r_moon|³
      (third-body correction for using an Earth-centred inertial frame)
    - Sun tidal (opt-in):    -μ_S · [(r − r_sun)/|r − r_sun|³ + r_sun/|r_sun|³]
      (validated against validate_mode_c.py:277-295)

    Parameters
    ----------
    t : float
        SPICE ephemeris time (ET).
    use_sun : bool
        If True, include solar third-body perturbation (default False).
    """
    r = y[:3]
    v = y[3:6]
    if moon_model == "circular":
        r_moon, _ = _moon_state_circular(t, t_arrive, frame)
    else:
        r_moon, _ = _moon_state_j2000(t)
    r_rel = r - r_moon
    r_norm = np.linalg.norm(r)
    r_rel_norm = np.linalg.norm(r_rel)
    a = (
        -MU_EARTH * r / r_norm**3
        - MU_MOON * r_rel / r_rel_norm**3
    )

    # Moon indirect term (third-body perturbation correction for
    # Earth-centred inertial frame).  This is a bug-fix: the direct
    # acceleration -μ_M · (r − r_moon)/|r − r_moon|³ alone is correct
    # only in a Moon-centred frame; in an Earth-centred frame the
    # indirect term -μ_M · r_moon/|r_moon|³ must be added.
    r_moon_norm = np.linalg.norm(r_moon)
    a += -MU_MOON * r_moon / (r_moon_norm ** 3)

    if use_sun:
        import spiceypy

        state_s, _ = spiceypy.spkezr("SUN", t, "J2000", "NONE", "EARTH")
        r_sun = np.asarray(state_s[:3], dtype=float)
        r_s2sc = r - r_sun
        r_sun_norm = np.linalg.norm(r_sun)
        r_s2sc_norm = np.linalg.norm(r_s2sc)
        a += -MU_SUN * (r_s2sc / (r_s2sc_norm ** 3) + r_sun / (r_sun_norm ** 3))

    return np.concatenate([v, a])


# ---------------------------------------------------------------------------
# Initial-state construction (perigee, true anomaly 0)
# ---------------------------------------------------------------------------
def initial_state_mep(
    p: tuple[float, float, float],
    h_p: float,
    i_deg: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Construct the perigee state (r, v) in MEP coordinates.

    Parameters
    ----------
    p : (v_p, Ω, ω)
        Perigee velocity (km/s, signed), ascending-node longitude (deg),
        argument of perigee (deg).
    h_p : float
        Perigee altitude (km).
    i_deg : float
        Initial orbit inclination (deg) relative to the MEP fundamental plane.

    Returns
    -------
    (r_mep, v_mep) : (3,) arrays
        Perigee position and velocity in MEP coordinates (km, km/s).
    """
    v_p, Omega_deg, omega_deg = p
    r_p = R_EARTH + h_p
    O = math.radians(Omega_deg)
    w = math.radians(omega_deg)
    i = math.radians(i_deg)
    cO, sO = math.cos(O), math.sin(O)
    cw, sw = math.cos(w), math.sin(w)
    ci, si = math.cos(i), math.sin(i)

    r_mep = r_p * np.array(
        [
            cO * cw - sO * sw * ci,
            sO * cw + cO * sw * ci,
            sw * si,
        ]
    )
    v_mep = v_p * np.array(
        [
            -cO * sw - sO * cw * ci,
            -sO * sw + cO * cw * ci,
            cw * si,
        ]
    )
    return r_mep, v_mep


# ---------------------------------------------------------------------------
# Forward propagation and q computation
# ---------------------------------------------------------------------------
@dataclass
class PropagationResult:
    """Result of a single forward propagation of the transfer arc."""

    q: np.ndarray  # (q1, q2, q3)
    perigee_lon_mep: float
    perigee_lat_mep: float
    perilune_lon_mep: float
    perilune_lat_mep: float
    v_perilune: float
    i_m: float
    C3: float
    r1_norm: float
    converged: bool = True
    p: tuple[float, float, float] = (0.0, 0.0, 0.0)  # solved (v_p, Ω, ω)


def propagate(
    p: tuple[float, float, float],
    h_p: float,
    i_deg: float,
    t_arrive: float,
    T: float,
    frame: mf.MEPFrame,
    moon_model: str = "spice",
    polar_ref: str = "mep",
    use_sun: bool = False,
) -> PropagationResult:
    """Forward-propagate from perigee to arrival and compute ``q`` and outputs.

    Parameters
    ----------
    p : (v_p, Ω, ω)
        Perigee velocity (km/s, signed), ascending-node longitude (deg),
        argument of perigee (deg).
    h_p : float
        Perigee altitude (km).
    i_deg : float
        Initial orbit inclination (deg) relative to the MEP fundamental plane.
    t_arrive : float
        Arrival epoch (ET).
    T : float
        Transfer time (seconds).
    frame : mf.MEPFrame
        The frozen MEP frame at ``t_arrive``.
    moon_model : str
        ``"spice"`` (default) uses the SPICE truth (non-circular orbit);
        ``"circular"`` uses the paper's circular-orbit assumption at the mean
        Earth-Moon distance 384,400 km (for the R-1V golden regression).
    polar_ref : {'mep', 'moon_equator'}
        Reference axis for q3 (and hence ``i_m``):
        ``'mep'`` (default, paper convention) uses the MEP Z axis (the white
        pole); ``'moon_equator'`` uses the IAU_MOON Z axis at ``t_arrive``
        expressed in J2000, i.e. the true lunar polar direction.  The switch
        makes the solver target the actual polar-lunar orbit directly
        (PRD v7 §四.1).
    use_sun : bool
        If True, include solar third-body perturbation in the dynamics
        (default False, backward-compatible).

    Returns
    -------
    PropagationResult
        The computed q vector and all downstream outputs.
    """
    r_mep, v_mep = initial_state_mep(p, h_p, i_deg)
    r0 = mf.from_mep(r_mep, frame)
    v0 = mf.from_mep(v_mep, frame)
    t0 = t_arrive - T
    y0 = np.concatenate([r0, v0])

    sol = solve_ivp(
        _rhs,
        (t0, t_arrive),
        y0,
        method="DOP853",
        rtol=RTOL,
        atol=ATOL,
        dense_output=True,
        args=(t_arrive, frame, moon_model, use_sun),
    )
    y1 = sol.y[:, -1]
    r1_j2000 = y1[:3]
    v1_j2000 = y1[3:6]

    # Moon-centred state at arrival.
    if moon_model == "circular":
        r_moon, v_moon = _moon_state_circular(t_arrive, t_arrive, frame)
    else:
        r_moon, v_moon = _moon_state_j2000(t_arrive)
    r1 = r1_j2000 - r_moon
    v1 = v1_j2000 - v_moon

    r1_norm = np.linalg.norm(r1)
    v1_norm = np.linalg.norm(v1)
    q1 = r1_norm - R_MOON  # perilune altitude (km)
    q2 = np.dot(r1, v1) / (r1_norm * v1_norm)
    h_vec = np.cross(r1, v1)
    h_norm = np.linalg.norm(h_vec)
    # q3 reference axis (PRD v7 §四.1): 'mep' = frozen white-pole normal,
    # 'moon_equator' = IAU_MOON Z axis at t_arrive in J2000 (true lunar pole).
    ref_axis = _polar_ref_axis(polar_ref, t_arrive, frame)
    q3 = np.dot(ref_axis, h_vec) / h_norm

    # Downstream outputs.
    r1_mep = mf.to_mep(r1, frame)
    perilune_lon, perilune_lat = mf.to_mep_lunar_lonlat(r1_mep)
    r_p_mep = mf.to_mep(r0, frame)
    perigee_lon, perigee_lat = mf.to_mep_lonlat(r_p_mep)
    i_m = math.degrees(math.acos(np.clip(q3, -1.0, 1.0)))
    C3 = v1_norm**2 - 2.0 * MU_MOON / r1_norm

    return PropagationResult(
        q=np.array([q1, q2, q3]),
        perigee_lon_mep=perigee_lon,
        perigee_lat_mep=perigee_lat,
        perilune_lon_mep=perilune_lon,
        perilune_lat_mep=perilune_lat,
        v_perilune=v1_norm,
        i_m=i_m,
        C3=C3,
        r1_norm=r1_norm,
    )


# ---------------------------------------------------------------------------
# Newton solver
# ---------------------------------------------------------------------------
def _residual(
    p: np.ndarray,
    h_p: float,
    i_deg: float,
    t_arrive: float,
    T: float,
    frame: mf.MEPFrame,
    q_target: np.ndarray,
    moon_model: str = "spice",
    polar_ref: str = "mep",
    use_sun: bool = False,
) -> np.ndarray:
    """Normalised residual f(p) - q_target.

    q1 is normalised by Q1_SCALE_KM (1000 km) to keep the Jacobian
    well-conditioned (PRD R-1).
    """
    res = propagate(
        tuple(p), h_p, i_deg, t_arrive, T, frame,
        moon_model=moon_model, polar_ref=polar_ref, use_sun=use_sun,
    )
    q = res.q
    q_norm = np.array([q[0] / Q1_SCALE_KM, q[1], q[2]])
    q_target_norm = np.array(
        [q_target[0] / Q1_SCALE_KM, q_target[1], q_target[2]]
    )
    return q_norm - q_target_norm


def solve_r1(
    h_p: float,
    T: float,
    i_deg: float,
    t_arrive: float,
    q_target: np.ndarray,
    p0: tuple[float, float, float],
    frame: Optional[mf.MEPFrame] = None,
    max_iter: int = 200,
    moon_model: str = "spice",
    polar_ref: str = "mep",
    use_sun: bool = False,
) -> PropagationResult:
    """Solve the R-1 lunar-side boundary value problem.

    Parameters
    ----------
    h_p : float
        Perigee altitude (km).
    T : float
        Transfer time (seconds).
    i_deg : float
        Initial orbit inclination (deg) relative to the MEP fundamental plane.
    t_arrive : float
        Arrival epoch (ET).
    q_target : (3,) array
        Target q vector: (H_m, 0, ±1) coplanar or (H_m, 0, 0) polar.
    p0 : (v_p, Ω, ω)
        Initial guess for the unknowns.
    frame : mf.MEPFrame, optional
        The frozen MEP frame at ``t_arrive``.  Built if not given.
    max_iter : int
        Maximum Newton iterations.
    moon_model : str
        ``"spice"`` (default) or ``"circular"`` (paper's circular-orbit
        assumption for the R-1V golden regression).
    polar_ref : {'mep', 'moon_equator'}
        Reference axis for q3 (see :func:`propagate`).

    Returns
    -------
    PropagationResult
        The converged solution with all downstream outputs.

    Raises
    ------
    RuntimeError
        If the Newton iteration fails to converge.
    """
    if frame is None:
        frame = mf.mep_at_epoch(t_arrive)

    def fun(p):
        return _residual(
            p, h_p, i_deg, t_arrive, T, frame, q_target,
            moon_model=moon_model, polar_ref=polar_ref, use_sun=use_sun,
        )

    # Use a tight xtol so the Newton iteration drives q2/q3 to ~1e-10.
    # (hybr does not accept ftol; xtol is the step tolerance.)
    sol = root(
        fun,
        np.asarray(p0, dtype=float),
        method="hybr",
        options={"xtol": 1e-14, "maxfev": max_iter * 3},
    )

    p_sol = tuple(sol.x)
    result = propagate(
        p_sol, h_p, i_deg, t_arrive, T, frame,
        moon_model=moon_model, polar_ref=polar_ref, use_sun=use_sun,
    )
    result.p = p_sol

    # Convergence check on the actual (unnormalised) q values.  This is the
    # authoritative criterion (PRD R-1: |Δq| < (0.1 km, 1e-8, 1e-8)); we do
    # not rely on sol.success alone because hybr can report "not making good
    # progress" even when the residual is already below tolerance.
    q = result.q
    q_target_arr = np.asarray(q_target, dtype=float)
    tol = np.array([0.1, 1e-8, 1e-8])  # (0.1 km, 1e-8, 1e-8)
    if not np.all(np.abs(q - q_target_arr) < tol):
        raise RuntimeError(
            f"R-1 solution did not meet q tolerance: q={q}, target={q_target_arr} "
            f"(root success={sol.success}: {sol.message})"
        )

    return result


# ---------------------------------------------------------------------------
# Coplanar 2-D solver (fixes ω = 0)
# ---------------------------------------------------------------------------
def solve_coplanar_2d(
    h_p: float,
    T: float,
    i_deg: float,
    t_arrive: float,
    p0: tuple[float, float],
    H_m: float = DEFAULT_H_M,
    prograde: bool = True,
    frame: Optional[mf.MEPFrame] = None,
    moon_model: str = "spice",
    max_iter: int = 200,
    polar_ref: str = "mep",
    use_sun: bool = False,
) -> PropagationResult:
    """Solve a coplanar R-1 case with the argument of perigee fixed at 0.

    For a coplanar orbit the Moon-centred orbit lies in the MEP fundamental
    plane, so ``q3 = ±1`` is automatically satisfied and the system reduces to
    two constraints (q1 = H_m, q2 = 0) in two unknowns (v_p, Ω) with ω = 0.
    Fixing ω removes the degeneracy of the full 3-unknown parameterisation and
    makes the Newton iteration converge reliably (used by the R-1V golden
    regression).

    Parameters
    ----------
    h_p : float
        Perigee altitude (km).
    T : float
        Transfer time (seconds).
    i_deg : float
        Initial orbit inclination (deg) relative to the MEP fundamental plane
        (0 for prograde, 180 for retrograde).
    t_arrive : float
        Arrival epoch (ET).
    p0 : (v_p, Ω)
        Initial guess for the perigee velocity (km/s) and ascending-node
        longitude (deg).
    H_m : float
        Target perilune altitude (km).
    prograde : bool
        True for q3 = +1 (prograde), False for q3 = −1 (retrograde).
    frame : mf.MEPFrame, optional
        The frozen MEP frame at ``t_arrive``.  Built if not given.
    moon_model : str
        ``"spice"`` (default) or ``"circular"``.
    max_iter : int
        Maximum Newton iterations.
    polar_ref : {'mep', 'moon_equator'}
        Reference axis for q3 (see :func:`propagate`); for a true coplanar
        orbit q3 is not part of the constraint set, so the value of
        ``polar_ref`` only affects the reported ``i_m`` and downstream
        diagnostics.

    Returns
    -------
    PropagationResult
        The converged solution with all downstream outputs.

    Raises
    ------
    RuntimeError
        If the Newton iteration fails to converge.
    """
    if frame is None:
        frame = mf.mep_at_epoch(t_arrive)

    q3_target = 1.0 if prograde else -1.0

    def fun(x):
        vp, Omega = x
        res = propagate(
            (vp, Omega, 0.0), h_p, i_deg, t_arrive, T, frame,
            moon_model=moon_model, polar_ref=polar_ref, use_sun=use_sun,
        )
        q = res.q
        return np.array(
            [
                q[0] / Q1_SCALE_KM - H_m / Q1_SCALE_KM,
                q[1] - 0.0,
            ]
        )

    # The coplanar 2-unknown system (q1, q2) does not constrain q3, so the
    # Newton iteration can land on either the prograde or retrograde branch.
    # We therefore multi-start from both the given v_p sign and its negation,
    # and keep the first solution whose q3 matches the requested branch.
    vp0, Om0 = p0
    starts = [(vp0, Om0), (-vp0, Om0)]
    last_err = None
    for vp_seed, Om_seed in starts:
        sol = root(
            fun,
            np.array([vp_seed, Om_seed], dtype=float),
            method="hybr",
            options={"xtol": 1e-14, "maxfev": max_iter * 3},
        )
        vp_sol, Omega_sol = sol.x
        p_sol = (vp_sol, Omega_sol, 0.0)
        result = propagate(
            p_sol,
            h_p,
            i_deg,
            t_arrive,
            T,
            frame,
            moon_model=moon_model,
            polar_ref=polar_ref,
            use_sun=use_sun,
        )
        result.p = p_sol

        # Convergence check on the actual (unnormalised) q values, including
        # q3 so we reject solutions on the wrong prograde/retrograde branch.
        q = result.q
        q_target_arr = np.array([H_m, 0.0, q3_target])
        tol = np.array([0.1, 1e-8, 1e-8])  # (0.1 km, 1e-8, 1e-8)
        if np.all(np.abs(q - q_target_arr) < tol):
            return result
        last_err = f"q={q}, target={q_target_arr} (root success={sol.success}: {sol.message})"

    raise RuntimeError(
        f"R-1 coplanar-2d solution did not meet q tolerance on either branch: "
        f"{last_err}"
    )


# ---------------------------------------------------------------------------
# Convenience wrappers for the two target families
# ---------------------------------------------------------------------------
def solve_coplanar(
    h_p: float,
    T: float,
    i_deg: float,
    t_arrive: float,
    p0: tuple[float, float, float],
    H_m: float = DEFAULT_H_M,
    prograde: bool = True,
    frame: Optional[mf.MEPFrame] = None,
    moon_model: str = "spice",
    polar_ref: str = "mep",
    use_sun: bool = False,
) -> PropagationResult:
    """Solve a coplanar R-1 case (q_target = (H_m, 0, ±1))."""
    q_target = np.array([H_m, 0.0, 1.0 if prograde else -1.0])
    return solve_r1(
        h_p, T, i_deg, t_arrive, q_target, p0, frame=frame,
        moon_model=moon_model, polar_ref=polar_ref, use_sun=use_sun,
    )


def solve_polar(
    h_p: float,
    T: float,
    i_deg: float,
    t_arrive: float,
    p0: tuple[float, float, float],
    H_m: float = DEFAULT_H_M,
    frame: Optional[mf.MEPFrame] = None,
    moon_model: str = "spice",
    polar_ref: str = "mep",
    use_sun: bool = False,
) -> PropagationResult:
    """Solve a polar-lunar R-1 case (q_target = (H_m, 0, 0)), for CE-7.

    Parameters
    ----------
    polar_ref : {'mep', 'moon_equator'}
        Reference axis for q3 (see :func:`propagate`).  Default ``'mep'``
        reproduces the paper's white-pole approximation; ``'moon_equator'``
        makes q3 = cos(angle between orbit normal and IAU_MOON Z axis), so
        a converged solution is polar w.r.t. the *true* lunar equator
        (PRD v7 §四.1).  Recommended for engineering work to avoid a
        post-solve plane-change dv (~190 m/s).
    """
    q_target = np.array([H_m, 0.0, 0.0])
    return solve_r1(
        h_p, T, i_deg, t_arrive, q_target, p0, frame=frame,
        moon_model=moon_model, polar_ref=polar_ref, use_sun=use_sun,
    )


def solve_non_coplanar(
    h_p: float,
    T: float,
    i_deg: float,
    t_arrive: float,
    p0: tuple[float, float, float],
    im_deg: float,
    H_m: float = DEFAULT_H_M,
    frame: Optional[mf.MEPFrame] = None,
    moon_model: str = "spice",
    polar_ref: str = "mep",
    use_sun: bool = False,
) -> PropagationResult:
    """Solve a non-coplanar R-1 case targeting lunar-orbit inclination ``im_deg``.

    For Yang paper Tables 2 & 3 (im = 85° or 95°, not 90°): the target is
    ``q_target = (H_m, 0, cos(im_deg))`` where ``cos(im_deg)`` is taken with
    the sign of the paper's reference perigee velocity (east-launch →
    positive q3 for im < 90°; the solver handles the sign via Newton).

    Parameters
    ----------
    im_deg : float
        Target lunar-orbit inclination (deg).  85° → q3 ≈ +0.0872;
        95° → q3 ≈ -0.0872.  90° reproduces the polar case (prefer
        :func:`solve_polar`).
    """
    q3_target = math.cos(math.radians(im_deg))
    q_target = np.array([H_m, 0.0, q3_target])
    return solve_r1(
        h_p, T, i_deg, t_arrive, q_target, p0, frame=frame,
        moon_model=moon_model, polar_ref=polar_ref, use_sun=use_sun,
    )


# ---------------------------------------------------------------------------
# Free-return BVP solver (R-1V-d)
# ---------------------------------------------------------------------------
def solve_free_return(
    h_p: float,
    T: float,
    i_deg: float,
    t_arrive: float,
    p0: tuple[float, float],
    frame: Optional[mf.MEPFrame] = None,
    moon_model: str = "circular",
    max_iter: int = 200,
    use_sun: bool = False,
) -> PropagationResult:
    """Solve the free-return BVP: perilune at 180° lunar longitude.

    R-1V-d is an independent functional check whose target is *not* the
    standard q_target of §3: the free-return orbit must reach perilune at
    lunar-surface longitude 180° (the far side, sub-anti-Earth point).  The
    unknowns are (v_p, λ_p) with ω = 0 (coplanar), and the two constraints
    are

        q2 = 0            (arrival state is the perilune)
        perilune_lon = 180°

    The paper's free-return orbits (east v_p=10.940288, λ_p=9.464357°,
    T=69.87 h; west v_p=10.943344, λ_p=-7.123997°, T=69.25 h) are computed
    with the circular lunar-orbit assumption, so ``moon_model`` defaults to
    ``"circular"``.

    Parameters
    ----------
    h_p : float
        Perigee altitude (km).
    T : float
        Transfer time (seconds).
    i_deg : float
        Initial orbit inclination (0 prograde / 180 retrograde).
    t_arrive : float
        Arrival epoch (ET).
    p0 : (v_p, λ_p)
        Initial guess for perigee velocity (km/s) and perigee longitude (deg).
    frame : mf.MEPFrame, optional
        The frozen MEP frame at ``t_arrive``.  Built if not given.
    moon_model : str
        ``"circular"`` (default, paper's assumption) or ``"spice"``.
    max_iter : int
        Maximum Newton iterations.

    Returns
    -------
    PropagationResult
        The converged solution with all downstream outputs.

    Raises
    ------
    RuntimeError
        If the Newton iteration fails to converge.
    """
    if frame is None:
        frame = mf.mep_at_epoch(t_arrive)

    def fun(x):
        vp, lp = x
        res = propagate(
            (vp, lp, 0.0), h_p, i_deg, t_arrive, T, frame,
            moon_model=moon_model, use_sun=use_sun,
        )
        q2 = res.q[1]
        dlon = ((res.perilune_lon_mep - 180.0 + 180.0) % 360.0) - 180.0
        return np.array([q2, dlon])

    sol = root(
        fun,
        np.asarray(p0, dtype=float),
        method="hybr",
        options={"xtol": 1e-14, "maxfev": max_iter * 3},
    )
    vp_sol, lp_sol = sol.x
    p_sol = (vp_sol, lp_sol, 0.0)
    result = propagate(
        p_sol, h_p, i_deg, t_arrive, T, frame,
        moon_model=moon_model, use_sun=use_sun,
    )
    result.p = p_sol

    # Convergence check on the actual constraints.
    q2 = result.q[1]
    dlon = ((result.perilune_lon_mep - 180.0 + 180.0) % 360.0) - 180.0
    if abs(q2) > 1e-8 or abs(dlon) > 1e-6:
        raise RuntimeError(
            f"R-1V-d free-return solution did not meet tolerance: "
            f"q2={q2:.3e}, dlon={dlon:.3e} deg "
            f"(root success={sol.success}: {sol.message})"
        )
    return result


def solve_free_return_3dof(
    h_p: float,
    i_deg: float,
    t_arrive: float,
    p0: tuple[float, float, float],
    H_m: float = DEFAULT_H_M,
    frame: Optional[mf.MEPFrame] = None,
    moon_model: str = "circular",
    max_iter: int = 300,
    use_sun: bool = False,
) -> PropagationResult:
    """Solve the free-return BVP with transfer time as an unknown (R-1V-d).

    The full free-return BVP has three unknowns (v_p, λ_p, T) and three
    constraints:

        q1 = H_m            (perilune altitude)
        q2 = 0              (arrival state is the perilune)
        perilune_lon = 180° (lunar far side, sub-anti-Earth point)

    This is the cleanest R-1V-d check: the solved (v_p, λ_p, T) are compared
    directly with the paper's free-return orbits (east v_p=10.940288,
    λ_p=9.464357°, T=69.87 h; west v_p=10.943344, λ_p=-7.123997°,
    T=69.25 h), with the PRD criteria |ΔT| < 0.5 h and perilune longitude
    180° ± 1°.

    Parameters
    ----------
    h_p : float
        Perigee altitude (km).
    i_deg : float
        Initial orbit inclination (0 prograde / 180 retrograde).
    t_arrive : float
        Arrival epoch (ET).
    p0 : (v_p, λ_p, T_h)
        Initial guess for perigee velocity (km/s), perigee longitude (deg)
        and transfer time (hours).
    H_m : float
        Target perilune altitude (km).
    frame : mf.MEPFrame, optional
        The frozen MEP frame at ``t_arrive``.  Built if not given.
    moon_model : str
        ``"circular"`` (default, paper's assumption) or ``"spice"``.
    max_iter : int
        Maximum Newton iterations.

    Returns
    -------
    PropagationResult
        The converged solution with all downstream outputs.  ``p`` holds
        ``(v_p, λ_p, T_h)``.

    Raises
    ------
    RuntimeError
        If the Newton iteration fails to converge.
    """
    if frame is None:
        frame = mf.mep_at_epoch(t_arrive)

    def fun(x):
        vp, lp, T_h = x
        T = T_h * 3600.0
        res = propagate(
            (vp, lp, 0.0), h_p, i_deg, t_arrive, T, frame,
            moon_model=moon_model, use_sun=use_sun,
        )
        q1, q2 = res.q[0], res.q[1]
        dlon = ((res.perilune_lon_mep - 180.0 + 180.0) % 360.0) - 180.0
        return np.array([q1 - H_m, q2, dlon])

    sol = root(
        fun,
        np.asarray(p0, dtype=float),
        method="hybr",
        options={"xtol": 1e-14, "maxfev": max_iter * 3},
    )
    vp_sol, lp_sol, T_h_sol = sol.x
    T_sol = T_h_sol * 3600.0
    p_sol = (vp_sol, lp_sol, 0.0)
    result = propagate(
        p_sol, h_p, i_deg, t_arrive, T_sol, frame,
        moon_model=moon_model, use_sun=use_sun,
    )
    result.p = (vp_sol, lp_sol, T_h_sol)

    # Convergence check on the actual constraints.
    q = result.q
    dlon = ((result.perilune_lon_mep - 180.0 + 180.0) % 360.0) - 180.0
    if abs(q[0] - H_m) > 0.1 or abs(q[1]) > 1e-8 or abs(dlon) > 1e-6:
        raise RuntimeError(
            f"R-1V-d free-return 3-DOF solution did not meet tolerance: "
            f"q={q}, dlon={dlon:.3e} deg "
            f"(root success={sol.success}: {sol.message})"
        )
    return result


# ---------------------------------------------------------------------------
# XYZ-trajectory sampling (R-1V-e: visualisation of the 16 paper orbits)
# ---------------------------------------------------------------------------
@dataclass
class TrajectorySamples:
    """Dense XYZ trajectory of a single transfer arc in MEP coordinates.

    All position / velocity arrays are expressed in the frozen MEP frame
    (``mep_frame.MEPFrame``) and use the Earth centre as origin, exactly as
    defined in the PRD §S-0.  The Moon reference trajectory (``r_moon_mep``)
    is provided at the *same* sample times so visualisations can render the
    Moon moving along its (sidereal mean motion) circular arc while the
    transfer arc unfolds.

    Attributes
    ----------
    p : (v_p, Ω, ω)
        Solved perigee state (km/s, deg, deg).
    t_s : (N,) ndarray
        Sample epochs in SPICE ET (seconds past J2000).
    r_mep : (N, 3) ndarray
        Spacecraft position in MEP XYZ (km), Earth-centre origin.
    v_mep : (N, 3) ndarray
        Spacecraft velocity in MEP XYZ (km/s).
    r_moon_mep : (N, 3) ndarray
        Moon-centre position in MEP XYZ (km) at the same sample times.
    result : PropagationResult
        The end-of-arc ``PropagationResult`` (carries q, lon/lat, C3, ...).
    """

    p: tuple[float, float, float]
    t_s: np.ndarray
    r_mep: np.ndarray
    v_mep: np.ndarray
    r_moon_mep: np.ndarray
    result: "PropagationResult"


def sample_trajectory(
    p: tuple[float, float, float],
    h_p: float,
    i_deg: float,
    t_arrive: float,
    T: float,
    frame: mf.MEPFrame,
    n_samples: int = 121,
    moon_model: str = "spice",
    use_sun: bool = False,
) -> TrajectorySamples:
    """Forward-propagate a transfer arc and return N samples in MEP XYZ.

    The integrator, tolerances and dynamics are identical to :func:`propagate`
    (DOP853, rtol=1e-12, atol=1e-9).  Only the additional step is asking
    ``scipy.integrate.solve_ivp`` for a dense output sampled at
    ``np.linspace(t0, t_arrive, n_samples)``.  All sample positions /
    velocities are rotated from J2000 into the *frozen* MEP frame, so the
    output is ready for plotting without any further coordinate work.

    Parameters
    ----------
    p : (v_p, Ω, ω)
        Perigee state (km/s, deg, deg).  Normally the converged solution
        from one of the ``solve_*`` functions, but the routine accepts any
        valid perigee state for ad-hoc visualisations.
    h_p : float
        Perigee altitude (km).
    i_deg : float
        Initial orbit inclination relative to the MEP fundamental plane
        (deg).
    t_arrive : float
        Arrival epoch (SPICE ET).
    T : float
        Transfer time (seconds).
    frame : mf.MEPFrame
        The frozen MEP frame at ``t_arrive``.
    n_samples : int
        Number of evenly spaced samples returned (default 121 ≈ one sample
        per 36 min on a 72 h arc).  Must be ≥ 2.
    moon_model : {"spice", "circular"}
        ``"spice"`` (default) uses the SPICE-truth non-circular lunar orbit;
        ``"circular"`` uses the paper's circular-orbit assumption (R-1V
        golden regression only).

    Returns
    -------
    TrajectorySamples
        Dataclass with N samples each of (t_s, r_mep, v_mep, r_moon_mep).
    """
    if n_samples < 2:
        raise ValueError(f"n_samples must be >= 2, got {n_samples}")

    # -- 1. Initial state in J2000, exactly as in propagate() --------------
    r_mep_p, v_mep_p = initial_state_mep(p, h_p, i_deg)
    r0 = mf.from_mep(r_mep_p, frame)
    v0 = mf.from_mep(v_mep_p, frame)
    t0 = t_arrive - T
    y0 = np.concatenate([r0, v0])
    t_eval = np.linspace(t0, t_arrive, n_samples)

    # -- 2. Forward integrate with dense_output + t_eval -------------------
    sol = solve_ivp(
        _rhs,
        (t0, t_arrive),
        y0,
        method="DOP853",
        rtol=RTOL,
        atol=ATOL,
        dense_output=True,
        t_eval=t_eval,
        args=(t_arrive, frame, moon_model, use_sun),
    )
    if sol.status != 0:
        raise RuntimeError(
            f"sample_trajectory: solve_ivp failed (status={sol.status}): "
            f"{sol.message}"
        )

    # sol.y has shape (6, N).  Convert each sample to MEP XYZ in one shot.
    r_j2000 = sol.y[:3, :]   # (3, N)
    v_j2000 = sol.y[3:6, :]  # (3, N)
    rot = frame.rot_j2000_to_mep  # rows = MEP basis vectors in J2000
    r_mep = (rot @ r_j2000).T    # (N, 3)
    v_mep = (rot @ v_j2000).T    # (N, 3)

    # -- 3. Moon position at the *same* sample times (MEP XYZ) -------------
    r_moon_mep = np.empty_like(r_mep)
    for k, t_k in enumerate(t_eval):
        if moon_model == "circular":
            r_m_j2000, _ = _moon_state_circular(t_k, t_arrive, frame)
        else:
            r_m_j2000, _ = _moon_state_j2000(t_k)
        r_moon_mep[k] = rot @ r_m_j2000

    # -- 4. End-of-arc metadata: piggyback on propagate() for q / lon / lat
    # This reuses the existing PropagationResult code path (1 extra DOP853
    # pass on the same initial condition, negligible cost) so the sampled
    # trajectory and the end-of-arc summary are guaranteed consistent.
    result = propagate(p, h_p, i_deg, t_arrive, T, frame, moon_model=moon_model, use_sun=use_sun)

    return TrajectorySamples(
        p=p,
        t_s=t_eval,
        r_mep=r_mep,
        v_mep=v_mep,
        r_moon_mep=r_moon_mep,
        result=result,
    )
