"""launch_window/lambert.py — Lambert's problem solver.

Implements the Lambert targeting used by the PRD v2 ``lambert_window.py``
engine.  Given two position vectors and a transfer time, it returns the
departure/arrival velocity vectors of the conic arc connecting them about a
central body.

The PRD v2 explicitly permits either ``poliastro`` or a hand-written Izzo
algorithm as the Lambert solver:

    "依赖:``spiceypy`` + ``numpy`` + 一个 Lambert 求解器(``poliastro`` 或
     自己写 Izzo 算法)。"

This module wraps poliastro's ``izzo`` implementation (Izzo 2015), which is
well-tested and validated against the standard Vallado test case.  It is a
pure numerical solver (no SPICE, no calibration constants); the only physical
constant is the gravitational parameter ``mu`` of the central body, passed by
the caller (Earth for the TLI leg, Moon for the arrival leg).

References
----------
- Izzo, D. (2015). Revisiting Lambert's problem. Celestial Mechanics and
  Dynamical Astronomy, 121(1), 1-15.
- Vallado, D. (2013). Fundamentals of Astrodynamics and Applications, 4th ed.
  Chapter 7 (Lambert's problem).
"""

from __future__ import annotations

import math

import numpy as np

_NITER = 1000
_RTOL = 1e-08


def lambert(
    r1: np.ndarray,
    r2: np.ndarray,
    tof: float,
    mu: float,
    prograde: bool = True,
    low_path: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    # Lazy import so the module loads without poliastro installed (only the
    # solver body actually needs it).
    from poliastro.core.iod import izzo

    r1 = np.asarray(r1, dtype=float)
    r2 = np.asarray(r2, dtype=float)

    if np.linalg.norm(r1) < 1e-09 or np.linalg.norm(r2) < 1e-09:
        raise ValueError("Position vectors must be non-zero")

    v1, v2 = izzo(
        mu,
        r1,
        r2,
        tof,
        0,
        bool(prograde),
        bool(low_path),
        _NITER,
        _RTOL,
    )

    return (np.asarray(v1, dtype=float), np.asarray(v2, dtype=float))


def lambert_initial_guess(
    r1_eci_km: np.ndarray,
    r2_eci_km: np.ndarray,
    dt_s: float,
    mu_earth_km3s2: float,
    prograde: bool = True,
    low_path: bool = True,
) -> np.ndarray:
    """Lambert initial-velocity generator: returns the TLI-side ECI velocity.

    Thin wrapper around :func:`lambert` that returns only ``v1`` (the
    departure velocity at ``r1``).  Implements the contract documented in
    REQ-LW-2026-LAMBERT-IG-001 §六.

    Parameters
    ----------
    r1_eci_km : (3,) array-like
        TLI-side position vector in J2000 (km).
    r2_eci_km : (3,) array-like
        Arrival-side (Moon) position vector in J2000 (km).
    dt_s : float
        Transfer time ``T_e - t_TLI`` (seconds).
    mu_earth_km3s2 : float
        Earth gravitational parameter (``tli_solver.MU_EARTH``).
    prograde : bool, default True
        Direction flag for the Earth→Moon transfer (independent of the
        Moon-centered segment inclination ``i_m``).
    low_path : bool, default True
        izzo ``low_path`` flag — selects the lower-energy branch.

    Returns
    -------
    (3,) ndarray
        Departure velocity vector at ``r1`` in J2000 (km/s).

    Raises
    ------
    RuntimeError
        If izzo fails to converge or the geometry is degenerate.  Callers
        are expected to catch and fall back to the family-table path.
    """
    try:
        v1, _v2 = lambert(
            r1_eci_km,
            r2_eci_km,
            float(dt_s),
            float(mu_earth_km3s2),
            prograde=bool(prograde),
            low_path=bool(low_path),
        )
    except RuntimeError:
        # izzo already raises RuntimeError on non-convergence / multi-rev
        # degeneracy — propagate verbatim per REQ contract.
        raise
    return np.asarray(v1, dtype=float)


def classical_from_state(
    r_eci_km: np.ndarray,
    v_eci_kms: np.ndarray,
    mu_km3s2: float,
) -> tuple[float, float, float, float]:
    """Classical Keplerian elements from a J2000 state vector.

    Implements the contract documented in REQ-LW-2026-LAMBERT-IG-001 §六
    (auxiliary).  Returns the inclination, RAAN, argument of perigee and
    semi-major axis — the three angles map directly onto the
    ``(inc_E, Omega_E, omega_E)`` consumed by ``_G_te``.

    Inclination uses ``acos(h_z / |h|)``; RAAN uses
    ``atan2(n_y, n_x)`` where ``n = k × h`` (line of nodes); the argument
    of perigee uses the standard two-argument form on the perifocal
    basis.  Angles are returned in degrees in ``[0, 360)``; ``a`` is in km.

    Parameters
    ----------
    r_eci_km : (3,) array-like
        Position vector in J2000 (km).
    v_eci_kms : (3,) array-like
        Velocity vector in J2000 (km/s).
    mu_km3s2 : float
        Gravitational parameter of the central body (km^3/s^2).

    Returns
    -------
    i_deg : float
        Inclination relative to the J2000 equator (deg, [0, 180]).
    raan_deg : float
        Right ascension of the ascending node (deg, [0, 360)).
    omega_deg : float
        Argument of perigee (deg, [0, 360)).
    a_km : float
        Semi-major axis (km; sign indicates bound vs unbound: negative
        for hyperbolic orbits).
    """
    r = np.asarray(r_eci_km, dtype=float)
    v = np.asarray(v_eci_kms, dtype=float)
    mu = float(mu_km3s2)

    r_norm = float(np.linalg.norm(r))
    v_norm2 = float(np.dot(v, v))
    if r_norm < 1e-09 or mu <= 0.0:
        raise ValueError("classical_from_state: r and mu must be non-zero")

    # Specific angular momentum h = r × v.
    h = np.cross(r, v)
    h_norm = float(np.linalg.norm(h))

    if h_norm < 1e-12:
        # Equatorial orbit: RAAN and omega are undefined.
        raise ValueError(
            "classical_from_state: angular momentum is zero (degenerate orbit)"
        )

    # Inclination: acos(h_z / |h|).
    inc = math.acos(max(-1.0, min(1.0, float(h[2]) / h_norm)))
    i_deg = math.degrees(inc)

    # Node vector n = k × h (line of nodes, points toward ascending node).
    k_hat = np.array([0.0, 0.0, 1.0])
    n_vec = np.cross(k_hat, h)
    n_norm = float(np.linalg.norm(n_vec))

    if n_norm < 1e-12:
        # Purely equatorial orbit: RAAN is undefined; default to 0.
        raan_deg = 0.0
        omega_deg = 0.0
    else:
        n_hat = n_vec / n_norm
        raan_rad = math.atan2(float(n_hat[1]), float(n_hat[0]))
        raan_deg = math.degrees(raan_rad) % 360.0

        # Argument of perigee: angle in the orbital plane from the
        # ascending node to the eccentricity vector, measured in the
        # direction of motion.
        # e = ((v^2 - mu/r) r - (r·v) v) / mu
        e_vec = ((v_norm2 - mu / r_norm) * r - float(np.dot(r, v)) * v) / mu
        e_norm = float(np.linalg.norm(e_vec))

        if e_norm < 1e-12:
            # Circular orbit: omega is undefined; default to 0.
            omega_deg = 0.0
        else:
            # Project e onto the perifocal basis (n̂, ĥ×n̂).  Use
            # atan2(e·(ĥ×n̂), e·n̂).
            e_n = float(np.dot(e_vec, n_hat))
            h_hat = h / h_norm
            m_hat = np.cross(h_hat, n_hat)
            e_m = float(np.dot(e_vec, m_hat))
            omega_rad = math.atan2(e_m, e_n)
            omega_deg = math.degrees(omega_rad) % 360.0

    # Semi-major axis via vis-viva: 1/a = 2/r - v^2/mu.
    inv_a = 2.0 / r_norm - v_norm2 / mu
    if abs(inv_a) < 1e-15:
        a_km = math.inf
    else:
        a_km = 1.0 / inv_a

    return float(i_deg), float(raan_deg), float(omega_deg), float(a_km)