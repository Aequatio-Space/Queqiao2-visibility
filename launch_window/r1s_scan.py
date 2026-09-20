"""launch_window/r1s_scan.py — R-1S initial-inclination scan (PRD v5 §5).

The paper's polar-lunar examples only use i0 = 5°/10° and claim "1° is
enough", but the engineering range is far wider: the ecliptic-vs-equator
angle swings 18.3°-28.6°, and Wenchang's latitude 19.614° forces a parking
orbit inclination >= 19.614° relative to the Earth equator, which converts
to i0 ~ 20°-45° relative to the MEP fundamental plane.

This module scans i0 in [1°, 45°] in 1° steps, for the four solution
families (ascending/descending x N/S perilune), and outputs the perigee and
perilune trajectory families.  The perigee-latitude spread vs i0 is the
reachable declination band that R-2 needs (the exclusion force of the whole
model).

Method
------
Each family is continued from the Table-4 seed at i0 = 5° (the four polar
solutions, which the R-1V-a regression reproduces to < 1 m/s).  Continuation
uses the previous converged (v_p, Omega, omega) as the next initial guess,
which is robust because the solution varies smoothly with i0.

Outputs
-------
- perigee_lon/lat vs i0 for each family  (R-2 input)
- perilune_lon/lat vs i0 for each family (R-3/R-6 input)
- reachable perigee-declination band: min/max perigee latitude over the
  scan, the quantity that bounds the model's exclusion force.

Usage::

    python -m launch_window.r1s_scan [--out r1s_scan.json]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from launch_window import geometry as g
from launch_window import mep_frame as mf
from launch_window import tli_solver as ts

_REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUT = _REPO_ROOT / "logs/launch_window/r1s_scan.json"

# Reference arrival epoch (same as the R-1V regression; the circular model is
# time-independent, so any fixed epoch reproduces the paper's geometry).
T_ARRIVE_UTC = "2026-09-01T00:00:00Z"

H_P = 200.0
H_M = ts.DEFAULT_H_M
T_H = 72.0  # paper's polar-lunar transfer time

I0_MIN = 1.0
I0_MAX = 114.0
I0_STEP = 1.0

# Table-4 seeds at i0 = 5° (from the R-1V-a regression, circular model).
# family key: (perilune hemisphere, ascending/descending at perigee)
FAMILY_SEEDS: dict[str, tuple[float, float, float]] = {
    "S_asc": (10.933858, 353.576, 14.819),   # perilune S, perigee ascending
    "N_desc": (10.933913, 9.485, 358.838),   # perilune N, perigee descending
    "N_asc": (10.933858, 173.576, 194.819),  # perilune N, perigee ascending
    "S_desc": (10.933913, 189.485, 178.838),  # perilune S, perigee descending
}


def find_seed_at_T(
    family: str,
    seed: tuple[float, float, float],
    t_arrive: float,
    frame: mf.MEPFrame,
    T_target_h: float,
    moon_model: str = "circular",
    T_start_h: float = T_H,
    dT_h: float = 1.0,
    H_m: float = H_M,
    polar_ref: str = "mep",
) -> tuple[float, float, float]:
    """Continue the family seed from ``T_start_h`` to ``T_target_h``.

    The Table-4 seeds are the T = 72 h polar solutions.  At other transfer
    times (e.g. CE-7's ~112 h) the 72 h seed does not converge directly, so
    the seed is continued in 1 h steps, each step using the previous
    converged (v_p, Omega, omega) as the initial guess.  The solution varies
    smoothly with T, so this is robust.

    ``polar_ref`` selects the q3 reference axis used by :func:`tli_solver.solve_polar`
    (default ``'mep'`` for backward compatibility).  ``'moon_equator'`` uses
    the IAU_MOON Z axis at ``t_arrive`` (true lunar polar direction); since
    that axis is rotated ~6.7° from the MEP Z, a previously converged
    ``'mep'`` seed is an excellent initial guess for ``'moon_equator'``.
    """
    guess = seed
    T_cur = T_start_h
    step = dT_h if T_target_h > T_start_h else -dT_h
    n_steps = int(abs(T_target_h - T_start_h) / dT_h + 0.5)
    for _ in range(n_steps):
        T_cur += step
        if abs(T_cur - T_target_h) < abs(step):
            T_cur = T_target_h
        res = ts.solve_polar(
            H_P,
            T_cur * 3600.0,
            5.0,  # seed is at i0 = 5 deg
            t_arrive,
            guess,
            H_m=H_m,
            frame=frame,
            moon_model=moon_model,
            polar_ref=polar_ref,
        )
        guess = res.p
    return guess


def scan_family(
    family: str,
    seed: tuple[float, float, float],
    t_arrive: float,
    frame: mf.MEPFrame,
    moon_model: str = "circular",
    T_h: float = T_H,
    H_m: float = H_M,
    polar_ref: str = "mep",
) -> list[dict]:
    """Continuation-scan one family over i0 in [1, 45] deg, step 1 deg.

    Returns a list of per-i0 records with the solved p and all downstream
    outputs.  The scan starts at the seed's i0 (5°) and continues both
    downward (5 -> 1) and upward (5 -> 45).  If ``T_h`` differs from the
    seed's transfer time (72 h), the seed is first continued in T to
    ``T_h`` (see :func:`find_seed_at_T`).

    ``polar_ref`` selects the q3 reference axis (default ``'mep'`` for
    backward compatibility; ``'moon_equator'`` uses the true IAU_MOON pole).
    """
    records: list[dict] = []
    T = T_h * 3600.0

    if abs(T_h - T_H) > 1e-9:
        seed = find_seed_at_T(
            family, seed, t_arrive, frame, T_h, moon_model, H_m=H_m, polar_ref=polar_ref
        )

    def solve_at(i0: float, guess: tuple[float, float, float]) -> tuple[dict, tuple[float, float, float]]:
        res = ts.solve_polar(
            H_P, T, i0, t_arrive, guess, H_m=H_m, frame=frame,
            moon_model=moon_model, polar_ref=polar_ref,
        )
        rec = {
            "i0_deg": i0,
            "family": family,
            "vp_km_s": float(res.p[0]),
            "Omega_deg": float(res.p[1]),
            "omega_deg": float(res.p[2]),
            "perigee_lon_deg": float(res.perigee_lon_mep),
            "perigee_lat_deg": float(res.perigee_lat_mep),
            "perilune_lon_deg": float(res.perilune_lon_mep),
            "perilune_lat_deg": float(res.perilune_lat_mep),
            "v_perilune_km_s": float(res.v_perilune),
            "i_m_deg": float(res.i_m),
            "C3": float(res.C3),
        }
        return rec, res.p

    # Downward continuation: 5 -> 1
    guess = seed
    for i0 in np.arange(5.0, I0_MIN - 1e-9, -I0_STEP):
        rec, guess = solve_at(float(i0), guess)
        records.append(rec)

    # Upward continuation: 6 -> 45
    guess = seed
    for i0 in np.arange(6.0, I0_MAX + 1e-9, I0_STEP):
        rec, guess = solve_at(float(i0), guess)
        records.append(rec)

    records.sort(key=lambda r: r["i0_deg"])
    return records


def run_scan(
    moon_model: str = "circular",
    T_h: float = T_H,
    family: str | None = None,
    H_m: float = H_M,
    polar_ref: str = "mep",
) -> dict:
    t_arrive = g.utc_to_et(T_ARRIVE_UTC)
    frame = mf.mep_at_epoch(t_arrive)
    mf.assert_moon_at_anti_x(frame)

    families: dict[str, list[dict]] = {}
    fam_items = [family] if family else list(FAMILY_SEEDS.keys())
    for fam in fam_items:
        families[fam] = scan_family(
            fam, FAMILY_SEEDS[fam], t_arrive, frame, moon_model,
            T_h=T_h, H_m=H_m, polar_ref=polar_ref,
        )

    # Reachable perigee-declination band: min/max perigee latitude over all
    # families and all i0.  This is the width of the R-2 reachable band.
    all_perigee_lat = [
        r["perigee_lat_deg"] for recs in families.values() for r in recs
    ]
    all_perilune_lat = [
        r["perilune_lat_deg"] for recs in families.values() for r in recs
    ]

    # Perigee-latitude spread at each i0 (the "lat_p ~ sin i0" law check).
    lat_spread = []
    for i0 in np.arange(I0_MIN, I0_MAX + 1e-9, I0_STEP):
        lats = [
            r["perigee_lat_deg"]
            for recs in families.values()
            for r in recs
            if abs(r["i0_deg"] - i0) < 1e-9
        ]
        if lats:
            lat_spread.append(
                {
                    "i0_deg": float(i0),
                    "min_perigee_lat_deg": float(min(lats)),
                    "max_perigee_lat_deg": float(max(lats)),
                    "spread_deg": float(max(lats) - min(lats)),
                }
            )

    return {
        "calibration": {
            "h_p_km": H_P,
            "H_m_km": H_m,
            "T_h": T_h,
            "moon_model": moon_model,
            "t_arrive_utc": T_ARRIVE_UTC,
            "i0_range_deg": [I0_MIN, I0_MAX],
            "i0_step_deg": I0_STEP,
        },
        "families": families,
        "reachable_band": {
            "perigee_lat_min_deg": float(min(all_perigee_lat)),
            "perigee_lat_max_deg": float(max(all_perigee_lat)),
            "perigee_lat_span_deg": float(max(all_perigee_lat) - min(all_perigee_lat)),
            "perilune_lat_min_deg": float(min(all_perilune_lat)),
            "perilune_lat_max_deg": float(max(all_perilune_lat)),
        },
        "perigee_lat_spread_vs_i0": lat_spread,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="R-1S initial-inclination scan")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--moon-model", default="circular", choices=["circular", "spice"])
    parser.add_argument(
        "--transfer-hours",
        type=float,
        default=T_H,
        help="transfer time (h) at which to scan the family (default 72)",
    )
    parser.add_argument(
        "--family",
        default=None,
        choices=list(FAMILY_SEEDS.keys()),
        help="scan only this family (default: all four)",
    )
    parser.add_argument(
        "--perilune-alt",
        type=float,
        default=H_M,
        help=f"perilune altitude H_m (km), default {H_M}",
    )
    parser.add_argument(
        "--polar-ref",
        default="mep",
        choices=("mep", "moon_equator"),
        help=(
            "q3 reference axis: 'mep' (default, frozen white-pole) or "
            "'moon_equator' (IAU_MOON Z, true lunar polar direction)."
        ),
    )
    args = parser.parse_args()

    result = run_scan(
        T_h=args.transfer_hours, moon_model=args.moon_model, family=args.family,
        H_m=args.perilune_alt, polar_ref=args.polar_ref,
    )
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)

    band = result["reachable_band"]
    first_fam = next(iter(result["families"]))
    print(f"R-1S: {len(result['families'])} families x "
          f"{len(result['families'][first_fam])} i0 steps")
    print(f"  reachable perigee-lat band: [{band['perigee_lat_min_deg']:.3f}, "
          f"{band['perigee_lat_max_deg']:.3f}] deg (span {band['perigee_lat_span_deg']:.3f})")
    print(f"  reachable perilune-lat band: [{band['perilune_lat_min_deg']:.3f}, "
          f"{band['perilune_lat_max_deg']:.3f}] deg")
    print(f"  wrote {args.out}")


if __name__ == "__main__":
    main()