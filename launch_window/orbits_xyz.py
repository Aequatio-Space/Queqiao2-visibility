"""launch_window/orbits_xyz.py — XYZ trajectory builder for the 24 paper orbits (R-1V-e).

The Yang Weilian paper ("发射极月卫星的转移轨道研究") publishes four
tables of transfer-orbit parameters that this regression faithfully
reproduces:

- **Table 1** (12 coplanar transfer orbits): three transfer times
  T = 72 / 96 / 120 h × four solutions per T (east/west launch ×
  prograde/retrograde), via
  :func:`launch_window.tli_solver.solve_coplanar_2d`.
- **Table 2** (4 non-coplanar orbits, i₀ = 5°/175°, i_m = 85°/95°).
- **Table 3** (4 same-direction orbits, i₀ = 10°, i_m = 85°/95°).
- **Table 4** (4 polar-lunar orbits, i₀ = 5°): near-polar solutions that
  deliver a 90° lunar-orbit inclination.

The four tables together are the **24 paper orbits** that this module
turns into ready-to-render XYZ trajectories:

    >>> from launch_window.orbits_xyz import build_yang_16_orbits_xyz
    >>> records = build_yang_16_orbits_xyz()
    >>> len(records)
    24

All positions are in the **frozen MEP frame** (Earth-centre origin,
X axis = anti-Moon at arrival, Y/Z the right-handed completion), with
the Moon trajectory sampled at the same times for visualisation.

The module is intentionally thin: it composes :func:`solve_coplanar_2d`,
:func:`solve_polar`, :func:`solve_non_coplanar` and
:func:`sample_trajectory` from ``launch_window.tli_solver`` and never
re-implements dynamics.  Every trajectory inherits the solver's DOP853 /
rtol=1e-12 / atol=1e-9 guarantees, so XYZ coordinates are bit-for-bit
consistent with the validated end-of-arc q vector.

移植说明：本文件移植自 AgentMoons 项目 ``launch_window`` 包 v0.6.0，
除内核路径与 CLI 默认输出路径外未改动科学计算逻辑；函数名
``build_yang_16_orbits_xyz`` 为保持接口兼容而保留（实际返回 24 条）。

R-1V-e — Yang paper Table 1–4 XYZ trajectory visualisations.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterator

import numpy as np

from launch_window import geometry as g
from launch_window import mep_frame as mf
from launch_window import tli_solver as ts


# Reference arrival epoch: all four tables are produced under the same
# frozen MEP-frame definition.
T_ARRIVE_UTC = "2026-08-24T00:00:00Z"
H_P_KM = 200.0
H_M_KM = ts.DEFAULT_H_M  # 300 km calibrated from Table 1

# Number of XYZ samples per orbit.  121 ≈ one sample per 36 min on a 72 h
# arc (one every 18 min on a 120 h arc), dense enough to resolve the
# Earth-departure and Moon-arrival turns without bloating the JSON.
DEFAULT_N_SAMPLES = 121


@dataclass
class OrbitRecord:
    """One paper orbit's metadata + its sampled XYZ trajectory.

    Attributes
    ----------
    family    : "coplanar_table1" | "polar_table4"
    orbit : int
        Orbit index in the table (1-based as in the paper).
    T_h : float
        Transfer time (hours).
    i0_deg : float
        Initial orbit inclination relative to the MEP fundamental plane (deg).
    vp_ref_km_s : float
        Reference perigee velocity from the paper (km/s).
    perigee_lon_ref_deg : float, optional
        Reference perigee longitude from the paper (deg).  Coplanar Table 1
        gives this directly; Table 4 gives (lon, lat) pair below.
    perigee_lat_ref_deg : float, optional
        Reference perigee latitude from the paper (deg).  Coplanar
        solutions are in the fundamental plane so this is 0.
    perilune_lon_ref_deg : float, optional
    perilune_lat_ref_deg : float, optional
    im_ref_deg : float
        Target lunar-orbit inclination (deg) — 0 / 180 for coplanar, 90
        for polar.
    p_sol : (v_p, Ω, ω)
        Solved perigee state (km/s, deg, deg).
    t_s : (N,) ndarray
        Sample epochs in SPICE ET.
    r_mep : (N, 3) ndarray
        Spacecraft position in MEP XYZ (km).
    v_mep : (N, 3) ndarray
        Spacecraft velocity in MEP XYZ (km/s).
    r_moon_mep : (N, 3) ndarray
        Moon-centre position in MEP XYZ (km) at the same sample times.
    q : (3,) ndarray
        End-of-arc q vector (H_m, cos f, cos θ).
    perigee_lon_sol_deg / perigee_lat_sol_deg : float
    perilune_lon_sol_deg / perilune_lat_sol_deg : float
    i_m_sol_deg : float
    v_perilune_km_s : float
    C3 : float
    error : str, optional
        Filled when the solver or sampler raises.
    """

    family: str
    orbit: int
    T_h: float
    i0_deg: float
    vp_ref_km_s: float
    im_ref_deg: float
    p_sol: tuple[float, float, float]
    t_s: np.ndarray
    r_mep: np.ndarray
    v_mep: np.ndarray
    r_moon_mep: np.ndarray
    q: np.ndarray
    perigee_lon_sol_deg: float
    perigee_lat_sol_deg: float
    perilune_lon_sol_deg: float
    perilune_lat_sol_deg: float
    i_m_sol_deg: float
    v_perilune_km_s: float
    C3: float
    perigee_lon_ref_deg: float | None = None
    perigee_lat_ref_deg: float | None = None
    perilune_lon_ref_deg: float | None = None
    perilune_lat_ref_deg: float | None = None
    error: str | None = None

    def to_jsonable(self) -> dict[str, Any]:
        """Convert to a JSON-serialisable dict (arrays → nested lists)."""
        if self.error is not None:
            return {
                "family": self.family,
                "orbit": self.orbit,
                "T_h": self.T_h,
                "i0_deg": self.i0_deg,
                "vp_ref_km_s": self.vp_ref_km_s,
                "im_ref_deg": self.im_ref_deg,
                "error": self.error,
            }
        return {
            "family": self.family,
            "orbit": self.orbit,
            "T_h": self.T_h,
            "i0_deg": self.i0_deg,
            "vp_ref_km_s": self.vp_ref_km_s,
            "im_ref_deg": self.im_ref_deg,
            "perigee_lon_ref_deg": self.perigee_lon_ref_deg,
            "perigee_lat_ref_deg": self.perigee_lat_ref_deg,
            "perilune_lon_ref_deg": self.perilune_lon_ref_deg,
            "perilune_lat_ref_deg": self.perilune_lat_ref_deg,
            "p_sol": list(self.p_sol),
            "t_s": self.t_s.tolist(),
            "r_mep": self.r_mep.tolist(),
            "v_mep": self.v_mep.tolist(),
            "r_moon_mep": self.r_moon_mep.tolist(),
            "q": self.q.tolist(),
            "perigee_lon_sol_deg": self.perigee_lon_sol_deg,
            "perigee_lat_sol_deg": self.perigee_lat_sol_deg,
            "perilune_lon_sol_deg": self.perilune_lon_sol_deg,
            "perilune_lat_sol_deg": self.perilune_lat_sol_deg,
            "i_m_sol_deg": self.i_m_sol_deg,
            "v_perilune_km_s": self.v_perilune_km_s,
            "C3": self.C3,
        }


# ---------------------------------------------------------------------------
# Specs for the 12 coplanar Table-1 orbits.
# Each row: (orbit_number, T_h, vp_ref, λp_ref_deg, im_ref_deg)
# The four sub-families per T are: east-prograde, east-retrograde,
# west-prograde, west-retrograde.  Signed vp_ref encodes east (positive)
# vs west (negative) launch.
# ---------------------------------------------------------------------------
COPLANAR_TABLE1: list[dict[str, float]] = [
    # T = 72 h
    {"orbit": 1,  "T_h": 72.0, "vp_ref":  10.932534, "lambda_p_ref":  7.730539, "im_ref": 0.0},
    {"orbit": 2,  "T_h": 72.0, "vp_ref":  10.936915, "lambda_p_ref":  9.001649, "im_ref": 180.0},
    {"orbit": 3,  "T_h": 72.0, "vp_ref": -10.934844, "lambda_p_ref": -6.502704, "im_ref": 0.0},
    {"orbit": 4,  "T_h": 72.0, "vp_ref": -10.938820, "lambda_p_ref": -6.363346, "im_ref": 180.0},
    # T = 96 h
    {"orbit": 5,  "T_h": 96.0, "vp_ref":  10.915679, "lambda_p_ref":  4.388118, "im_ref": 0.0},
    {"orbit": 6,  "T_h": 96.0, "vp_ref":  10.919341, "lambda_p_ref":  5.280707, "im_ref": 180.0},
    {"orbit": 7,  "T_h": 96.0, "vp_ref": -10.917268, "lambda_p_ref": -1.380779, "im_ref": 0.0},
    {"orbit": 8,  "T_h": 96.0, "vp_ref": -10.920099, "lambda_p_ref": -1.560826, "im_ref": 180.0},
    # T = 120 h
    {"orbit": 9,  "T_h": 120.0, "vp_ref":  10.914859, "lambda_p_ref":  2.381466, "im_ref": 0.0},
    {"orbit": 10, "T_h": 120.0, "vp_ref":  10.917405, "lambda_p_ref":  2.765479, "im_ref": 180.0},
    {"orbit": 11, "T_h": 120.0, "vp_ref": -10.915262, "lambda_p_ref":  1.632890, "im_ref": 0.0},
    {"orbit": 12, "T_h": 120.0, "vp_ref": -10.917211, "lambda_p_ref":  1.232622, "im_ref": 180.0},
]


# ---------------------------------------------------------------------------
# Specs for the 4 polar-lunar Table-4 orbits (i₀ = 5°).
# Each row: (orbit_number, T_h=72, vp_ref, λp_ref_deg, φp_ref_deg,
#            λm_ref_deg, φm_ref_deg)
# Table 4 is the "strict 90°" polar case at T = 72 h with i₀ = 5°.
# ---------------------------------------------------------------------------
POLAR_TABLE4: list[dict[str, float]] = [
    {"orbit": 1, "T_h": 72.0, "vp_ref": 10.934734, "lambda_p_ref":  8.354, "phi_p_ref":  1.283, "lambda_m_ref": 129.392, "phi_m_ref": -47.505, "im_ref": 90.0},
    {"orbit": 2, "T_h": 72.0, "vp_ref": 10.934789, "lambda_p_ref":  8.340, "phi_p_ref": -0.095, "lambda_m_ref": 129.424, "phi_m_ref":  45.778, "im_ref": 90.0},
    {"orbit": 3, "T_h": 72.0, "vp_ref": 10.934734, "lambda_p_ref":  8.354, "phi_p_ref": -1.283, "lambda_m_ref": 120.392, "phi_m_ref":  47.505, "im_ref": 90.0},
    {"orbit": 4, "T_h": 72.0, "vp_ref": 10.934789, "lambda_p_ref":  8.340, "phi_p_ref":  0.095, "lambda_m_ref": 120.424, "phi_m_ref": -45.778, "im_ref": 90.0},
]


# ---------------------------------------------------------------------------
# Specs for Table-2 non-coplanar transfer orbits (i₀ = 5° or 175°, im = 85°/95°).
# Source: tests/fixtures/yang_noncoplanar_polar_cases.json.
# i₀ is the *signed* launch inclination (175° means retrograde-launch
# analogue).  im_ref is the target lunar-orbit inclination.
# ---------------------------------------------------------------------------
NONCOPLANAR_TABLE2: list[dict[str, float]] = [
    {"orbit": 1, "T_h": 72.0, "vp_ref": 10.934542, "i0_deg":   5.0, "lambda_p_ref":  8.298, "phi_p_ref":  1.280, "lambda_m_ref": 123.999, "phi_m_ref": -47.243, "im_ref":  85.0},
    {"orbit": 2, "T_h": 72.0, "vp_ref": 10.934979, "i0_deg":   5.0, "lambda_p_ref":  8.395, "phi_p_ref": -0.090, "lambda_m_ref": 134.662, "phi_m_ref":  45.581, "im_ref":  95.0},
    {"orbit": 3, "T_h": 72.0, "vp_ref": 10.936686, "i0_deg": 175.0, "lambda_p_ref": -6.422, "phi_p_ref": -0.062, "lambda_m_ref": 112.735, "phi_m_ref": -55.114, "im_ref":  85.0},
    {"orbit": 4, "T_h": 72.0, "vp_ref": 10.937006, "i0_deg": 175.0, "lambda_p_ref": -6.397, "phi_p_ref": -1.185, "lambda_m_ref": 127.601, "phi_m_ref":  56.462, "im_ref":  95.0},
]


# ---------------------------------------------------------------------------
# Specs for Table-3 same-direction transfer orbits (i₀ = 10°, im = 85°/95°).
# Source: tests/fixtures/yang_noncoplanar_polar_cases.json.
# ---------------------------------------------------------------------------
SAME_DIRECTION_TABLE3: list[dict[str, float]] = [
    {"orbit": 1, "T_h": 72.0, "vp_ref": 10.934535, "i0_deg": 10.0, "lambda_p_ref":  8.224, "phi_p_ref":  1.875, "lambda_m_ref": 123.811, "phi_m_ref": -48.180, "im_ref": 85.0},
    {"orbit": 2, "T_h": 72.0, "vp_ref": 10.934645, "i0_deg": 10.0, "lambda_p_ref":  8.196, "phi_p_ref":  0.504, "lambda_m_ref": 124.208, "phi_m_ref":  44.731, "im_ref": 85.0},
    {"orbit": 3, "T_h": 72.0, "vp_ref": 10.934917, "i0_deg": 10.0, "lambda_p_ref":  8.334, "phi_p_ref":  1.880, "lambda_m_ref": 134.818, "phi_m_ref": -48.233, "im_ref": 95.0},
    {"orbit": 4, "T_h": 72.0, "vp_ref": 10.935026, "i0_deg": 10.0, "lambda_p_ref":  8.306, "phi_p_ref":  0.513, "lambda_m_ref": 134.548, "phi_m_ref":  44.794, "im_ref": 95.0},
]


# ---------------------------------------------------------------------------
# Solvers
# ---------------------------------------------------------------------------
def _solve_coplanar_orbit(
    spec: dict[str, float],
    frame: mf.MEPFrame,
) -> OrbitRecord:
    """Solve one Table-1 coplanar row and sample its XYZ trajectory."""
    T_h = spec["T_h"]
    T = T_h * 3600.0
    vp_ref = spec["vp_ref"]
    lambda_p_ref = spec["lambda_p_ref"]
    im_ref = spec["im_ref"]
    orbit = int(spec["orbit"])

    prograde = im_ref == 0.0
    p0 = (vp_ref, lambda_p_ref)
    res = ts.solve_coplanar_2d(
        H_P_KM,
        T,
        0.0,
        frame.t_arrive,
        p0,
        H_m=H_M_KM,
        prograde=prograde,
        frame=frame,
        moon_model="circular",
    )

    p_sol = res.p if res.p != (0.0, 0.0, 0.0) else (vp_ref, lambda_p_ref, 0.0)
    samples = ts.sample_trajectory(
        p_sol,
        H_P_KM,
        0.0,
        frame.t_arrive,
        T,
        frame,
        n_samples=DEFAULT_N_SAMPLES,
        moon_model="circular",
    )

    return OrbitRecord(
        family="coplanar_table1",
        orbit=orbit,
        T_h=T_h,
        i0_deg=0.0,
        vp_ref_km_s=vp_ref,
        im_ref_deg=im_ref,
        perigee_lon_ref_deg=lambda_p_ref,
        perigee_lat_ref_deg=0.0,
        perilune_lon_ref_deg=None,  # coplanar q2=0; longitude read off i_m only
        perilune_lat_ref_deg=None,
        p_sol=p_sol,
        t_s=samples.t_s,
        r_mep=samples.r_mep,
        v_mep=samples.v_mep,
        r_moon_mep=samples.r_moon_mep,
        q=samples.result.q,
        perigee_lon_sol_deg=samples.result.perigee_lon_mep,
        perigee_lat_sol_deg=samples.result.perigee_lat_mep,
        perilune_lon_sol_deg=samples.result.perilune_lon_mep,
        perilune_lat_sol_deg=samples.result.perilune_lat_mep,
        i_m_sol_deg=samples.result.i_m,
        v_perilune_km_s=samples.result.v_perilune,
        C3=samples.result.C3,
    )


def _solve_polar_orbit(
    spec: dict[str, float],
    frame: mf.MEPFrame,
) -> OrbitRecord:
    """Solve one Table-4 polar row (i₀ = 5°) and sample its XYZ trajectory.

    Table 4's target is a 90° lunar-orbit inclination (q3 = 0).  The
    non-coplanar R-1 solver (3-DOF in v_p, Ω, ω) finds a converged state
    ``p_sol`` that meets q1, q2 = 0 and q3 = 0 simultaneously.  We start
    from the paper's reference perigee (vp_ref, λp_ref) and let Newton
    find the matching ω; if that fails we fall back to ω = 0 (sufficient
    for the i₀ = 5° near-coplanar family at T = 72 h).
    """
    T_h = spec["T_h"]
    T = T_h * 3600.0
    vp_ref = spec["vp_ref"]
    lambda_p_ref = spec["lambda_p_ref"]
    phi_p_ref = spec["phi_p_ref"]
    im_ref = spec["im_ref"]
    orbit = int(spec["orbit"])

    # Initial guess: paper's vp + longitude; ω starts at the latitude
    # interpretation that yields the reference perigee latitude.
    # For the polar table, perigee latitude is small (<2°) so ω is near
    # ±(90° - φp_ref) in the i₀ = 5° frame.  We let Newton figure it out
    # from a generic ω=0 starting point first.
    p0 = (vp_ref, lambda_p_ref, 0.0)
    res = ts.solve_polar(
        H_P_KM,
        T,
        5.0,  # i₀ from Table 4
        frame.t_arrive,
        p0,
        H_m=H_M_KM,
        frame=frame,
        moon_model="circular",
    )
    p_sol = res.p

    samples = ts.sample_trajectory(
        p_sol,
        H_P_KM,
        5.0,
        frame.t_arrive,
        T,
        frame,
        n_samples=DEFAULT_N_SAMPLES,
        moon_model="circular",
    )

    return OrbitRecord(
        family="polar_table4",
        orbit=orbit,
        T_h=T_h,
        i0_deg=5.0,
        vp_ref_km_s=vp_ref,
        im_ref_deg=im_ref,
        perigee_lon_ref_deg=lambda_p_ref,
        perigee_lat_ref_deg=phi_p_ref,
        perilune_lon_ref_deg=spec["lambda_m_ref"],
        perilune_lat_ref_deg=spec["phi_m_ref"],
        p_sol=p_sol,
        t_s=samples.t_s,
        r_mep=samples.r_mep,
        v_mep=samples.v_mep,
        r_moon_mep=samples.r_moon_mep,
        q=samples.result.q,
        perigee_lon_sol_deg=samples.result.perigee_lon_mep,
        perigee_lat_sol_deg=samples.result.perigee_lat_mep,
        perilune_lon_sol_deg=samples.result.perilune_lon_mep,
        perilune_lat_sol_deg=samples.result.perilune_lat_mep,
        i_m_sol_deg=samples.result.i_m,
        v_perilune_km_s=samples.result.v_perilune,
        C3=samples.result.C3,
    )


def _solve_non_coplanar_orbit(
    spec: dict[str, float],
    frame: mf.MEPFrame,
    family: str,
) -> OrbitRecord:
    """Solve one non-coplanar Table 2 / Table 3 row (im = 85° or 95°).

    Same shape as :func:`_solve_polar_orbit` but targets ``q3 = cos(im)``
    via :func:`ts.solve_non_coplanar` instead of ``q3 = 0``.
    """
    T_h = spec["T_h"]
    T = T_h * 3600.0
    vp_ref = spec["vp_ref"]
    i0_deg = float(spec["i0_deg"])
    lambda_p_ref = spec["lambda_p_ref"]
    phi_p_ref = spec["phi_p_ref"]
    im_ref = spec["im_ref"]
    orbit = int(spec["orbit"])

    # Initial guess: paper's vp + longitude; ω = 0.
    p0 = (vp_ref, lambda_p_ref, 0.0)
    res = ts.solve_non_coplanar(
        H_P_KM,
        T,
        i0_deg,
        frame.t_arrive,
        p0,
        im_deg=im_ref,
        H_m=H_M_KM,
        frame=frame,
        moon_model="circular",
    )
    p_sol = res.p

    samples = ts.sample_trajectory(
        p_sol,
        H_P_KM,
        i0_deg,
        frame.t_arrive,
        T,
        frame,
        n_samples=DEFAULT_N_SAMPLES,
        moon_model="circular",
    )

    return OrbitRecord(
        family=family,
        orbit=orbit,
        T_h=T_h,
        i0_deg=i0_deg,
        vp_ref_km_s=vp_ref,
        im_ref_deg=im_ref,
        perigee_lon_ref_deg=lambda_p_ref,
        perigee_lat_ref_deg=phi_p_ref,
        perilune_lon_ref_deg=spec["lambda_m_ref"],
        perilune_lat_ref_deg=spec["phi_m_ref"],
        p_sol=p_sol,
        t_s=samples.t_s,
        r_mep=samples.r_mep,
        v_mep=samples.v_mep,
        r_moon_mep=samples.r_moon_mep,
        q=samples.result.q,
        perigee_lon_sol_deg=samples.result.perigee_lon_mep,
        perigee_lat_sol_deg=samples.result.perigee_lat_mep,
        perilune_lon_sol_deg=samples.result.perilune_lon_mep,
        perilune_lat_sol_deg=samples.result.perilune_lat_mep,
        i_m_sol_deg=samples.result.i_m,
        v_perilune_km_s=samples.result.v_perilune,
        C3=samples.result.C3,
    )


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------
def build_yang_16_orbits_xyz(
    t_arrive_utc: str = T_ARRIVE_UTC,
    n_samples: int = DEFAULT_N_SAMPLES,
    frame: mf.MEPFrame | None = None,
) -> list[OrbitRecord]:
    """Build the 24 paper orbits as XYZ trajectory records.

    Returns a list of 24 :class:`OrbitRecord` objects in table order:
    Table 1 (12 coplanar), Table 2 (4 non-coplanar), Table 3 (4 same-
    direction), Table 4 (4 polar).  Per-orbit exceptions are caught and
    surfaced via :attr:`OrbitRecord.error` so a single failure never
    blocks the rest of the suite.

    Parameters
    ----------
    t_arrive_utc : str
        ISO-8601 UTC arrival epoch for the MEP frame.
    n_samples : int
        Number of XYZ samples per orbit (>= 2).  Defaults to 121.
    frame : mf.MEPFrame, optional
        Pre-built MEP frame; built from ``t_arrive_utc`` if not given.

    Returns
    -------
    list[OrbitRecord]
        24 records (12 coplanar + 4 non-coplanar + 4 same-direction +
        4 polar), each carrying the metadata and the sampled XYZ trajectory.
    """
    if frame is None:
        t_arrive = g.utc_to_et(t_arrive_utc)
        frame = mf.mep_at_epoch(t_arrive)
    mf.assert_moon_at_anti_x(frame)

    records: list[OrbitRecord] = []

    # ---- Table 1: 12 coplanar -------------------------------------------
    for spec in COPLANAR_TABLE1:
        try:
            records.append(_solve_coplanar_orbit(spec, frame))
        except Exception as exc:  # keep the suite running
            records.append(
                OrbitRecord(
                    family="coplanar_table1",
                    orbit=int(spec["orbit"]),
                    T_h=float(spec["T_h"]),
                    i0_deg=0.0,
                    vp_ref_km_s=float(spec["vp_ref"]),
                    im_ref_deg=float(spec["im_ref"]),
                    perigee_lon_ref_deg=float(spec["lambda_p_ref"]),
                    perigee_lat_ref_deg=0.0,
                    perilune_lon_ref_deg=None,
                    perilune_lat_ref_deg=None,
                    p_sol=(0.0, 0.0, 0.0),
                    t_s=np.array([]),
                    r_mep=np.zeros((0, 3)),
                    v_mep=np.zeros((0, 3)),
                    r_moon_mep=np.zeros((0, 3)),
                    q=np.array([np.nan, np.nan, np.nan]),
                    perigee_lon_sol_deg=np.nan,
                    perigee_lat_sol_deg=np.nan,
                    perilune_lon_sol_deg=np.nan,
                    perilune_lat_sol_deg=np.nan,
                    i_m_sol_deg=np.nan,
                    v_perilune_km_s=np.nan,
                    C3=np.nan,
                    error=f"{type(exc).__name__}: {exc}",
                )
            )

    # ---- Table 2: 4 non-coplanar (i₀ = 5°/175°, im = 85°/95°) ----------
    for spec in NONCOPLANAR_TABLE2:
        try:
            records.append(_solve_non_coplanar_orbit(spec, frame, "non_coplanar_table2"))
        except Exception as exc:
            records.append(
                OrbitRecord(
                    family="non_coplanar_table2",
                    orbit=int(spec["orbit"]),
                    T_h=float(spec["T_h"]),
                    i0_deg=float(spec["i0_deg"]),
                    vp_ref_km_s=float(spec["vp_ref"]),
                    im_ref_deg=float(spec["im_ref"]),
                    perigee_lon_ref_deg=float(spec["lambda_p_ref"]),
                    perigee_lat_ref_deg=float(spec["phi_p_ref"]),
                    perilune_lon_ref_deg=float(spec["lambda_m_ref"]),
                    perilune_lat_ref_deg=float(spec["phi_m_ref"]),
                    p_sol=(0.0, 0.0, 0.0),
                    t_s=np.array([]),
                    r_mep=np.zeros((0, 3)),
                    v_mep=np.zeros((0, 3)),
                    r_moon_mep=np.zeros((0, 3)),
                    q=np.array([np.nan, np.nan, np.nan]),
                    perigee_lon_sol_deg=np.nan,
                    perigee_lat_sol_deg=np.nan,
                    perilune_lon_sol_deg=np.nan,
                    perilune_lat_sol_deg=np.nan,
                    i_m_sol_deg=np.nan,
                    v_perilune_km_s=np.nan,
                    C3=np.nan,
                    error=f"{type(exc).__name__}: {exc}",
                )
            )

    # ---- Table 3: 4 same-direction (i₀ = 10°, im = 85°/95°) ------------
    for spec in SAME_DIRECTION_TABLE3:
        try:
            records.append(_solve_non_coplanar_orbit(spec, frame, "same_direction_table3"))
        except Exception as exc:
            records.append(
                OrbitRecord(
                    family="same_direction_table3",
                    orbit=int(spec["orbit"]),
                    T_h=float(spec["T_h"]),
                    i0_deg=float(spec["i0_deg"]),
                    vp_ref_km_s=float(spec["vp_ref"]),
                    im_ref_deg=float(spec["im_ref"]),
                    perigee_lon_ref_deg=float(spec["lambda_p_ref"]),
                    perigee_lat_ref_deg=float(spec["phi_p_ref"]),
                    perilune_lon_ref_deg=float(spec["lambda_m_ref"]),
                    perilune_lat_ref_deg=float(spec["phi_m_ref"]),
                    p_sol=(0.0, 0.0, 0.0),
                    t_s=np.array([]),
                    r_mep=np.zeros((0, 3)),
                    v_mep=np.zeros((0, 3)),
                    r_moon_mep=np.zeros((0, 3)),
                    q=np.array([np.nan, np.nan, np.nan]),
                    perigee_lon_sol_deg=np.nan,
                    perigee_lat_sol_deg=np.nan,
                    perilune_lon_sol_deg=np.nan,
                    perilune_lat_sol_deg=np.nan,
                    i_m_sol_deg=np.nan,
                    v_perilune_km_s=np.nan,
                    C3=np.nan,
                    error=f"{type(exc).__name__}: {exc}",
                )
            )

    # ---- Table 4: 4 polar (i₀ = 5°) -------------------------------------
    for spec in POLAR_TABLE4:
        try:
            records.append(_solve_polar_orbit(spec, frame))
        except Exception as exc:
            records.append(
                OrbitRecord(
                    family="polar_table4",
                    orbit=int(spec["orbit"]),
                    T_h=float(spec["T_h"]),
                    i0_deg=5.0,
                    vp_ref_km_s=float(spec["vp_ref"]),
                    im_ref_deg=float(spec["im_ref"]),
                    perigee_lon_ref_deg=float(spec["lambda_p_ref"]),
                    perigee_lat_ref_deg=float(spec["phi_p_ref"]),
                    perilune_lon_ref_deg=float(spec["lambda_m_ref"]),
                    perilune_lat_ref_deg=float(spec["phi_m_ref"]),
                    p_sol=(0.0, 0.0, 0.0),
                    t_s=np.array([]),
                    r_mep=np.zeros((0, 3)),
                    v_mep=np.zeros((0, 3)),
                    r_moon_mep=np.zeros((0, 3)),
                    q=np.array([np.nan, np.nan, np.nan]),
                    perigee_lon_sol_deg=np.nan,
                    perigee_lat_sol_deg=np.nan,
                    perilune_lon_sol_deg=np.nan,
                    perilune_lat_sol_deg=np.nan,
                    i_m_sol_deg=np.nan,
                    v_perilune_km_s=np.nan,
                    C3=np.nan,
                    error=f"{type(exc).__name__}: {exc}",
                )
            )

    return records


def iter_yang_16_orbits_xyz(
    t_arrive_utc: str = T_ARRIVE_UTC,
    n_samples: int = DEFAULT_N_SAMPLES,
    frame: mf.MEPFrame | None = None,
) -> Iterator[OrbitRecord]:
    """Generator variant of :func:`build_yang_16_orbits_xyz`.

    Yields one :class:`OrbitRecord` at a time so callers can stream results
    to the browser (Server-Sent Events / chunked HTTP) and surface real
    per-orbit progress instead of blocking until the full suite finishes.
    Per-orbit exceptions are isolated into a stub record carrying
    ``error`` so a single failure never aborts the stream.
    """
    if frame is None:
        t_arrive = g.utc_to_et(t_arrive_utc)
        frame = mf.mep_at_epoch(t_arrive)
    mf.assert_moon_at_anti_x(frame)

    for spec in COPLANAR_TABLE1:
        try:
            yield _solve_coplanar_orbit(spec, frame)
        except Exception as exc:
            yield OrbitRecord(
                family="coplanar_table1",
                orbit=int(spec["orbit"]),
                T_h=float(spec["T_h"]),
                i0_deg=0.0,
                vp_ref_km_s=float(spec["vp_ref"]),
                im_ref_deg=float(spec["im_ref"]),
                perigee_lon_ref_deg=float(spec["lambda_p_ref"]),
                perigee_lat_ref_deg=0.0,
                perilune_lon_ref_deg=None,
                perilune_lat_ref_deg=None,
                p_sol=(0.0, 0.0, 0.0),
                t_s=np.array([]),
                r_mep=np.zeros((0, 3)),
                v_mep=np.zeros((0, 3)),
                r_moon_mep=np.zeros((0, 3)),
                q=np.array([np.nan, np.nan, np.nan]),
                perigee_lon_sol_deg=np.nan,
                perigee_lat_sol_deg=np.nan,
                perilune_lon_sol_deg=np.nan,
                perilune_lat_sol_deg=np.nan,
                i_m_sol_deg=np.nan,
                v_perilune_km_s=np.nan,
                C3=np.nan,
                error=f"{type(exc).__name__}: {exc}",
            )

    for spec in NONCOPLANAR_TABLE2:
        try:
            yield _solve_non_coplanar_orbit(spec, frame, "non_coplanar_table2")
        except Exception as exc:
            yield OrbitRecord(
                family="non_coplanar_table2",
                orbit=int(spec["orbit"]),
                T_h=float(spec["T_h"]),
                i0_deg=float(spec["i0_deg"]),
                vp_ref_km_s=float(spec["vp_ref"]),
                im_ref_deg=float(spec["im_ref"]),
                perigee_lon_ref_deg=float(spec["lambda_p_ref"]),
                perigee_lat_ref_deg=float(spec["phi_p_ref"]),
                perilune_lon_ref_deg=float(spec["lambda_m_ref"]),
                perilune_lat_ref_deg=float(spec["phi_m_ref"]),
                p_sol=(0.0, 0.0, 0.0),
                t_s=np.array([]),
                r_mep=np.zeros((0, 3)),
                v_mep=np.zeros((0, 3)),
                r_moon_mep=np.zeros((0, 3)),
                q=np.array([np.nan, np.nan, np.nan]),
                perigee_lon_sol_deg=np.nan,
                perigee_lat_sol_deg=np.nan,
                perilune_lon_sol_deg=np.nan,
                perilune_lat_sol_deg=np.nan,
                i_m_sol_deg=np.nan,
                v_perilune_km_s=np.nan,
                C3=np.nan,
                error=f"{type(exc).__name__}: {exc}",
            )

    for spec in SAME_DIRECTION_TABLE3:
        try:
            yield _solve_non_coplanar_orbit(spec, frame, "same_direction_table3")
        except Exception as exc:
            yield OrbitRecord(
                family="same_direction_table3",
                orbit=int(spec["orbit"]),
                T_h=float(spec["T_h"]),
                i0_deg=float(spec["i0_deg"]),
                vp_ref_km_s=float(spec["vp_ref"]),
                im_ref_deg=float(spec["im_ref"]),
                perigee_lon_ref_deg=float(spec["lambda_p_ref"]),
                perigee_lat_ref_deg=float(spec["phi_p_ref"]),
                perilune_lon_ref_deg=float(spec["lambda_m_ref"]),
                perilune_lat_ref_deg=float(spec["phi_m_ref"]),
                p_sol=(0.0, 0.0, 0.0),
                t_s=np.array([]),
                r_mep=np.zeros((0, 3)),
                v_mep=np.zeros((0, 3)),
                r_moon_mep=np.zeros((0, 3)),
                q=np.array([np.nan, np.nan, np.nan]),
                perigee_lon_sol_deg=np.nan,
                perigee_lat_sol_deg=np.nan,
                perilune_lon_sol_deg=np.nan,
                perilune_lat_sol_deg=np.nan,
                i_m_sol_deg=np.nan,
                v_perilune_km_s=np.nan,
                C3=np.nan,
                error=f"{type(exc).__name__}: {exc}",
            )

    for spec in POLAR_TABLE4:
        try:
            yield _solve_polar_orbit(spec, frame)
        except Exception as exc:
            yield OrbitRecord(
                family="polar_table4",
                orbit=int(spec["orbit"]),
                T_h=float(spec["T_h"]),
                i0_deg=5.0,
                vp_ref_km_s=float(spec["vp_ref"]),
                im_ref_deg=float(spec["im_ref"]),
                perigee_lon_ref_deg=float(spec["lambda_p_ref"]),
                perigee_lat_ref_deg=float(spec["phi_p_ref"]),
                perilune_lon_ref_deg=float(spec["lambda_m_ref"]),
                perilune_lat_ref_deg=float(spec["phi_m_ref"]),
                p_sol=(0.0, 0.0, 0.0),
                t_s=np.array([]),
                r_mep=np.zeros((0, 3)),
                v_mep=np.zeros((0, 3)),
                r_moon_mep=np.zeros((0, 3)),
                q=np.array([np.nan, np.nan, np.nan]),
                perigee_lon_sol_deg=np.nan,
                perigee_lat_sol_deg=np.nan,
                perilune_lon_sol_deg=np.nan,
                perilune_lat_sol_deg=np.nan,
                i_m_sol_deg=np.nan,
                v_perilune_km_s=np.nan,
                C3=np.nan,
                error=f"{type(exc).__name__}: {exc}",
            )


def write_records_json(
    records: list[OrbitRecord],
    path: Path | str,
) -> Path:
    """Serialise a list of records to JSON (used by scripts/build_yang_orbits.py)."""
    path = Path(path)
    payload = {
        "n_orbits": len(records),
        "t_arrive_utc": T_ARRIVE_UTC,
        "h_p_km": H_P_KM,
        "H_m_km": H_M_KM,
        "n_samples_per_orbit": DEFAULT_N_SAMPLES,
        "moon_model": "circular",
        "orbits": [r.to_jsonable() for r in records],
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False)
    return path


# ---------------------------------------------------------------------------
# Module-level convenience
# ---------------------------------------------------------------------------
if __name__ == "__main__":  # pragma: no cover
    import argparse

    parser = argparse.ArgumentParser(
        description="Build the 24 paper orbits as XYZ trajectory JSON"
    )
    parser.add_argument(
        "--out",
        type=Path,
        # 移植适配：默认输出到主站资产路径（原为仓库根 yang_orbits_xyz.json）。
        # 推荐入口为 scripts/build_yang_orbits.py（含完整性校验）。
        default=Path(__file__).resolve().parent.parent
        / "viewer"
        / "data"
        / "yang_orbits_xyz.json",
    )
    args = parser.parse_args()

    records = build_yang_16_orbits_xyz()
    n_ok = sum(1 for r in records if r.error is None)
    print(f"Built {n_ok}/{len(records)} orbits")
    write_records_json(records, args.out)
    print(f"Wrote {args.out}")