"""launch_window/scan_im.py — arbitrary lunar-orbit-inclination family table builder.

Extends :mod:`launch_window.r1s_scan` (which is locked to the polar target
``i_m = 90°``) to build R-1S perigee-family tables for an *arbitrary* target
lunar-orbit inclination ``i_m``.

Method
------
For each of the four solution families (``S_asc``, ``N_desc``, ``N_asc``,
``S_desc``):

1. Start from the Table-4 polar seed at ``i0 = 5°`` (``i_m = 90°``).
2. Continue the seed in transfer time from 72 h to ``T_h`` (reusing
   :func:`r1s_scan.find_seed_at_T`).
3. Continue ``i_m`` from 90° to the target in ``im_step_deg`` steps, solving
   R-1 with ``q_target = (H_m, 0, cos i_m)`` and using the previous converged
   ``(v_p, Ω, ω)`` as the next initial guess.
4. At the target ``i_m``, scan ``i0`` in [1°, 45°] in 1° steps by
   continuation (same as :func:`r1s_scan.scan_family`).

Output
------
A JSON file with the same schema as ``r1s_scan.json`` (``calibration``,
``families``, ``reachable_band``) so that
:func:`perigee_match.load_families` can consume it directly.

Usage::

    python -m launch_window.scan_im --i-m 43.7 --transfer-hours 112.47 \\
        --perilune-alt 200 --out r1s_scan_im437_T11247_Hm200.json
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from launch_window import geometry as g
from launch_window import mep_frame as mf
from launch_window import tli_solver as ts
from launch_window.r1s_scan import (
    FAMILY_SEEDS,
    H_P,
    I0_MAX,
    I0_MIN,
    I0_STEP,
    T_H,
    find_seed_at_T,
)

_REPO_ROOT = Path(__file__).resolve().parent.parent

# Reference arrival epoch (same as the R-1V regression; the circular model is
# time-independent, so any fixed epoch reproduces the paper's geometry).
T_ARRIVE_UTC = "2026-08-24T00:00:00Z"

# Default perilune altitude (km).  Matches the CE-7 engineering tables
# (``r1s_scan_T*_Hm200_meq.json``) and the B-3 CE-6 nominal (200 km / 3 burns).
DEFAULT_H_M_KM = 200.0

# Continuation step in i_m (deg) when moving from the polar seed (i_m = 90°)
# to the target inclination.
IM_STEP_DEG = 1.0


def solve_im(
    h_p: float,
    T: float,
    i0_deg: float,
    t_arrive: float,
    p0: tuple[float, float, float],
    i_m_deg: float,
    H_m: float,
    frame: mf.MEPFrame,
    moon_model: str = "circular",
    polar_ref: str = "mep",
    use_sun: bool = False,
) -> ts.PropagationResult:
    """Solve the R-1 BVP with an arbitrary lunar-orbit inclination target.

    The target q vector is ``(H_m, 0, cos i_m)``: q3 = cos θ is the cosine
    of the angle between the Moon-centred orbital-plane normal and the
    polar-reference axis, so ``i_m = acos(q3)``.  This generalises
    ``solve_polar`` (q3 = 0 ⇒ i_m = 90°) and the coplanar targets
    (q3 = ±1 ⇒ i_m = 0°/180°).

    When ``use_sun`` is True, solar third-body perturbation is included
    in the dynamics (DE440 SPICE ephemeris).  Recommended for engineering
    work; increases solve time ~2×.
    """
    q3 = math.cos(math.radians(i_m_deg))
    q_target = np.array([H_m, 0.0, q3])
    return ts.solve_r1(
        h_p, T, i0_deg, t_arrive, q_target, p0, frame=frame,
        moon_model=moon_model, polar_ref=polar_ref, use_sun=use_sun,
    )


def continue_im(
    family: str,
    seed: tuple[float, float, float],
    t_arrive: float,
    frame: mf.MEPFrame,
    i_m_target_deg: float,
    T_h: float,
    H_m: float,
    moon_model: str = "circular",
    polar_ref: str = "mep",
    im_step_deg: float = IM_STEP_DEG,
    skip_t_continuation: bool = False,
    use_sun: bool = False,
) -> tuple[float, float, float]:
    """Continue the polar seed (i_m = 90°) to the target i_m at i0 = 5°.

    First continues the seed in transfer time from 72 h to ``T_h`` (at
    i_m = 90°, reusing :func:`r1s_scan.find_seed_at_T`), then continues
    i_m from 90° to ``i_m_target_deg`` in ``im_step_deg`` steps, each step
    using the previous converged ``(v_p, Ω, ω)`` as the initial guess.

    A *direct jump* from i_m = 90° to the target is attempted first: the
    R-1 Newton iteration is robust enough that a single solve at the target
    i_m converges from the polar seed in ~1 s (measured), which is far
    cheaper than the ~70-step 1° continuation.  Only if the direct jump
    fails do we fall back to the stepwise continuation.

    When ``skip_t_continuation`` is True the seed is assumed to already be
    at the target ``T_h`` (e.g. loaded from a polar table at that T_h), so
    the 72 h → T_h continuation is skipped.
    """
    if skip_t_continuation:
        guess = seed
    else:
        guess = find_seed_at_T(
            family, seed, t_arrive, frame, T_h, moon_model, H_m=H_m, polar_ref=polar_ref
        )
    T = T_h * 3600.0

    # When the target polar_ref is 'moon_equator', the seed returned by
    # find_seed_at_T is already at i_m = 90° w.r.t. the IAU_MOON Z axis
    # (solve_polar targets q3 = 0 against that axis), so no re-anchoring is
    # needed — the continuation below starts from i_m = 90° in the target
    # frame and steps to i_m_target_deg.

    # Direct jump first.
    try:
        res = solve_im(
            H_P, T, 5.0, t_arrive, guess, i_m_target_deg, H_m, frame,
            moon_model=moon_model, polar_ref=polar_ref, use_sun=use_sun,
        )
        return res.p
    except RuntimeError:
        pass  # fall through to stepwise continuation

    i_m_cur = 90.0
    step = im_step_deg if i_m_target_deg > 90.0 else -im_step_deg
    n_steps = int(abs(i_m_target_deg - 90.0) / im_step_deg + 0.5)
    for _ in range(n_steps):
        i_m_cur += step
        if abs(i_m_cur - i_m_target_deg) < abs(step):
            i_m_cur = i_m_target_deg
        res = solve_im(
            H_P, T, 5.0, t_arrive, guess, i_m_cur, H_m, frame,
            moon_model=moon_model, polar_ref=polar_ref, use_sun=use_sun,
        )
        guess = res.p
    return guess


def load_polar_seed(
    family: str,
    seed_table: Path,
    i0_deg: float = 5.0,
) -> tuple[float, float, float]:
    """Load the i0 = 5° polar seed for one family from an existing table.

    The seed table must be a ``r1s_scan*.json``-schema file (with
    ``families[family]`` records carrying ``vp_km_s``, ``Omega_deg``,
    ``omega_deg``) at the target transfer time and i_m = 90°.  This avoids
    the slow 72 h → T_h continuation when a suitable polar table already
    exists (e.g. ``r1s_scan_T1125_Hm200_meq.json``).
    """
    import json as _json

    with open(seed_table, encoding="utf-8") as f:
        data = _json.load(f)
    recs = data["families"][family]
    rec = next(r for r in recs if abs(r["i0_deg"] - i0_deg) < 1e-9)
    return (float(rec["vp_km_s"]), float(rec["Omega_deg"]), float(rec["omega_deg"]))


def scan_family_im(
    family: str,
    seed: tuple[float, float, float],
    t_arrive: float,
    frame: mf.MEPFrame,
    i_m_deg: float,
    moon_model: str = "circular",
    T_h: float = T_H,
    H_m: float = DEFAULT_H_M_KM,
    polar_ref: str = "mep",
    im_step_deg: float = IM_STEP_DEG,
    seed_table: Path | None = None,
    use_sun: bool = False,
    i0_min: float = I0_MIN,
    i0_max: float = I0_MAX,
    i0_step: float = I0_STEP,
) -> list[dict]:
    """Continuation-scan one family over i0 in [i0_min, i0_max] at a fixed i_m.

    Returns a list of per-i0 records with the solved p and all downstream
    outputs, in the same schema as :func:`r1s_scan.scan_family`.

    When ``seed_table`` is given, the i0 = 5° polar seed is loaded from that
    table (which must be at the target ``T_h`` and i_m = 90°) instead of
    continuing the Table-4 seed from 72 h — much faster.

    Notes
    -----
    The continuation seed always sits at i0 = 5° (the Table-4 polar seeds /
    the i0 = 5° record of ``seed_table``); the grid is walked downward from
    5° to ``i0_min`` and upward from ``5° + i0_step`` to ``i0_max``.
    """
    records: list[dict] = []
    T = T_h * 3600.0

    if seed_table is not None:
        seed = load_polar_seed(family, seed_table, i0_deg=5.0)
        # The seed is already at i_m = 90° and the target T_h; jump to i_m.
        seed = continue_im(
            family, seed, t_arrive, frame, i_m_deg, T_h, H_m,
            moon_model=moon_model, polar_ref=polar_ref, im_step_deg=im_step_deg,
            skip_t_continuation=True, use_sun=use_sun,
        )
    else:
        # Continue the seed to (T_h, i_m) at i0 = 5°.
        seed = continue_im(
            family, seed, t_arrive, frame, i_m_deg, T_h, H_m,
            moon_model=moon_model, polar_ref=polar_ref, im_step_deg=im_step_deg,
            use_sun=use_sun,
        )

    def solve_at(
        i0: float, guess: tuple[float, float, float]
    ) -> tuple[dict, tuple[float, float, float]]:
        res = solve_im(
            H_P, T, i0, t_arrive, guess, i_m_deg, H_m, frame,
            moon_model=moon_model, polar_ref=polar_ref, use_sun=use_sun,
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

    # Downward continuation: 5 -> i0_min
    guess = seed
    for i0 in np.arange(5.0, i0_min - 1e-9, -i0_step):
        rec, guess = solve_at(float(i0), guess)
        records.append(rec)

    # Upward continuation: 5 + step -> i0_max
    guess = seed
    for i0 in np.arange(5.0 + i0_step, i0_max + 1e-9, i0_step):
        rec, guess = solve_at(float(i0), guess)
        records.append(rec)

    records.sort(key=lambda r: r["i0_deg"])
    return records


# ---------------------------------------------------------------------------
# MEP -> IAU_MOON point clouds and family deduplication
# ---------------------------------------------------------------------------
# Behavioural tolerance for "complete duplicate" detection.  Two records are
# treated as the same solved family member when every downstream field
# agrees within this tolerance (degrees for angles, km/s for velocities).
# The scan's smallest structural separations at the CE-7 operating point
# (T = 112.61 h, SPICE + Sun, i_m = 90°) are O(3e-2) deg between adjacent
# i0 samples and O(2e-1) deg across families, so 1e-6 has a 4-decade safety
# margin while still absorbing the ~1e-8-level solver noise between two
# converged determinations of the same physical orbit.
DUP_TOL = 1e-6


def mep_to_iau_moon_lonlat(
    lon_deg: float,
    lat_deg: float,
    frame: mf.MEPFrame,
) -> tuple[float, float]:
    """Rotate a MEP-frame direction (lon/lat, deg) into IAU_MOON lon/lat.

    The perigee and perilune records of the R-1S family are directions in
    the frozen MEP frame (Earth-side and Moon-side respectively).  The
    conversion is a pure rotation of the direction vector:

        MEP -> J2000 (frozen ``MEPFrame``) -> IAU_MOON (SPICE ``pxform``)

    evaluated at the frame epoch ``frame.t_arrive``.  The Moon-centred
    origin is irrelevant for a direction, so the same formula serves both
    the Earth-side perigee and the Moon-side perilune.
    """
    import spiceypy

    v_mep = mf.from_mep_lonlat(lon_deg, lat_deg, 1.0)
    v_j2000 = mf.from_mep(v_mep, frame)
    rot = spiceypy.pxform("J2000", "IAU_MOON", float(frame.t_arrive))
    v_iau = rot @ v_j2000
    lon = math.degrees(math.atan2(float(v_iau[1]), float(v_iau[0]))) % 360.0
    lat = math.degrees(math.asin(np.clip(float(v_iau[2]), -1.0, 1.0)))
    return lon, lat


def build_point_clouds(
    families: dict[str, list[dict]],
    frame: mf.MEPFrame,
    *,
    coord_frame: str = "iau_moon",
) -> dict[str, list[dict]]:
    """Build the four lon/lat point clouds for the scan families.

    Returns a dict with the four intermediate products requested for the
    i0 scan:

    - ``S_perigee``  : S-hemisphere families, Earth-side perigee points
    - ``S_perilune`` : S-hemisphere families, Moon-side perilune points
    - ``N_perigee``  : N-hemisphere families, Earth-side perigee points
    - ``N_perilune`` : N-hemisphere families, Moon-side perilune points

    Each cloud is a list of ``{family, i0_deg, lon_deg, lat_deg}`` records.
    ``coord_frame`` selects the output frame:

    - ``"iau_moon"`` (default): MEP directions rotated to the true lunar
      body frame via ``pxform('J2000','IAU_MOON',t_arrive)`` (what the
      CE-7 products use);
    - ``"mep"``: the raw MEP lunar-surface lon/lat of the family records —
      the paper's convention (needed to compare against the published
      figures 18-22 bounds).

    The S/N grouping follows the family name prefix (``S_asc``/``S_desc`` ->
    ``S``, ``N_asc``/``N_desc`` -> ``N``).
    """
    if coord_frame not in ("iau_moon", "mep"):
        raise ValueError(
            f"coord_frame must be 'iau_moon' or 'mep', got {coord_frame!r}"
        )
    clouds: dict[str, list[dict]] = {
        "S_perigee": [], "S_perilune": [],
        "N_perigee": [], "N_perilune": [],
    }
    for fam, recs in families.items():
        hemi = "S" if fam.upper().startswith("S") else "N"
        for r in recs:
            if coord_frame == "mep":
                p_lon = float(r["perigee_lon_deg"])
                p_lat = float(r["perigee_lat_deg"])
                m_lon = float(r["perilune_lon_deg"])
                m_lat = float(r["perilune_lat_deg"])
            else:
                p_lon, p_lat = mep_to_iau_moon_lonlat(
                    float(r["perigee_lon_deg"]), float(r["perigee_lat_deg"]),
                    frame,
                )
                m_lon, m_lat = mep_to_iau_moon_lonlat(
                    float(r["perilune_lon_deg"]), float(r["perilune_lat_deg"]),
                    frame,
                )
            clouds[f"{hemi}_perigee"].append({
                "family": fam, "i0_deg": float(r["i0_deg"]),
                "lon_deg": p_lon, "lat_deg": p_lat,
            })
            clouds[f"{hemi}_perilune"].append({
                "family": fam, "i0_deg": float(r["i0_deg"]),
                "lon_deg": m_lon, "lat_deg": m_lat,
            })
    return clouds


def build_iau_moon_point_clouds(
    families: dict[str, list[dict]],
    frame: mf.MEPFrame,
) -> dict[str, list[dict]]:
    """Backward-compatible alias for ``coord_frame='iau_moon'``."""
    return build_point_clouds(families, frame, coord_frame="iau_moon")


POINT_CLOUD_NAMES: tuple[str, ...] = (
    "S_perigee", "S_perilune", "N_perigee", "N_perilune",
)


def _circular_lon_limits(
    lon_deg: list[float],
    pad_frac: float = 0.06,
    min_pad_deg: float = 0.2,
    full_span_deg: float = 350.0,
) -> tuple[float, float]:
    """Smallest longitude window [lo, hi] covering the points (deg).

    Points are treated on the 360° circle: the window is the complement of
    the largest angular gap, so a cloud that straddles the 0°/360° seam
    (the perigee arcs) is zoomed to its real arc instead of the full
    0..360 range.  A nearly full circle (span >= ``full_span_deg``) falls
    back to [0, 360].  Seam-crossing windows are recentred into
    (-180, 180] so longitudes read naturally (negative values near 0°);
    the plotting routine draws ±360° copies of every point so nothing is
    clipped at the seam.
    """
    if not lon_deg:
        return 0.0, 360.0
    vals = sorted(float(v) % 360.0 for v in lon_deg)
    n = len(vals)
    if n == 1:
        return vals[0] - 5.0, vals[0] + 5.0
    gaps = [(vals[(i + 1) % n] - vals[i]) % 360.0 for i in range(n)]
    k = max(range(n), key=lambda i: gaps[i])
    start = vals[(k + 1) % n]
    span = 360.0 - gaps[k]
    if span >= full_span_deg:
        return 0.0, 360.0
    pad = max(min_pad_deg, span * pad_frac)
    window = (start - pad, start + span + pad)
    # Recentre a seam-crossing window into (-180, 180] so longitudes read
    # naturally (paper-convention MEP figures use small negative values);
    # the plotting routine draws ±360° copies so nothing is clipped.
    if 0.5 * (window[0] + window[1]) > 180.0:
        window = (window[0] - 360.0, window[1] - 360.0)
    return window


def _padded_limits(
    values: list[float],
    frac: float = 0.12,
    min_pad: float = 0.05,
    fallback: tuple[float, float] = (-90.0, 90.0),
) -> tuple[float, float]:
    """[min, max] of ``values`` padded by ``frac`` of the span."""
    if not values:
        return fallback
    lo, hi = min(values), max(values)
    pad = max(min_pad, (hi - lo) * frac)
    return lo - pad, hi + pad


def plot_iau_moon_point_clouds(
    clouds: dict,
    out_path: Path | str,
    *,
    title: str | None = None,
    subtitle: str | None = None,
    frame_label: str = "IAU_MOON",
    panel_boxes: dict[str, tuple[float, float, float, float]] | None = None,
    axis_limits: dict[str, tuple[float, float, float, float]] | None = None,
    box_pad_frac: float = 0.12,
    dpi: int = 200,
) -> Path:
    """Render the four point clouds to a 2x2 lon/lat PNG.

    Optional companion to :func:`build_point_clouds` (also called from
    :func:`assemble_scan_result` when ``point_cloud_plot`` is set).
    Each panel is one cloud (S/N hemisphere x perigee/perilune); the
    markers are coloured by the parking-orbit inclination ``i0`` on a
    shared colourbar and shaped by the ascending/descending family branch.
    By default the lon/lat windows adapt per panel (circular-arc fit for
    longitude, padded min/max for latitude) so each cloud fills its panel
    instead of being a speck on a fixed 0-360 / -90-90 canvas.

    ``frame_label`` names the coordinate frame in the axis labels/title
    (``"IAU_MOON"`` for :func:`build_point_clouds` output, ``"MEP"`` for
    the paper-convention clouds).

    ``panel_boxes`` optionally pins each panel's window to a reference
    bounding box ``(x_lo, x_hi, y_lo, y_hi)`` keyed by cloud name (e.g.
    the published figures 18-22 bounds).  The axes then span the box
    padded by ``box_pad_frac`` and the box itself is drawn dashed, so
    points spilling outside the reference range remain visible.

    ``axis_limits`` pins each panel's window to an explicit
    ``(x_lo, x_hi, y_lo, y_hi)`` without drawing any reference box; use it
    to keep the axes fixed across the frames of an animation (see
    ``viz_pointcloud_drift``).  It is mutually exclusive with
    ``panel_boxes``.

    ``clouds`` is the payload produced by :func:`build_point_clouds` —
    metadata keys such as ``frame`` / ``epoch_utc`` / ``note`` (present in
    the ``run_scan`` result) are ignored.  Matplotlib is imported lazily
    inside this function so the scan itself does not require it; the Agg
    backend is forced for headless use.

    Returns the written PNG path.
    """
    if panel_boxes is not None and axis_limits is not None:
        raise ValueError(
            "pass either panel_boxes or axis_limits to "
            "plot_iau_moon_point_clouds, not both"
        )
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import font_manager as _fm
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.ticker import FuncFormatter, MaxNLocator

    # CJK font: best-effort (falls back silently), same convention as
    # polar_window / stage_calendars.  Only available families are put into
    # rcParams so matplotlib does not warn for every missing candidate.
    for _fname in ("AlibabaPuHuiTi-3-55-Regular.ttf",
                   "AlibabaPuHuiTi-3-85-Bold.ttf"):
        try:
            p = Path.home() / "Library" / "Fonts" / _fname
            if p.exists():
                _fm.fontManager.addfont(str(p))
        except Exception:
            pass
    _available = {f.name for f in _fm.fontManager.ttflist}
    _cjk = next(
        (name for name in (
            "Alibaba PuHuiTi 3.0", "PingFang SC", "PingFang HK",
            "Hiragino Sans GB", "Heiti TC", "Arial Unicode MS",
            "STHeiti", "SimHei",
        ) if name in _available),
        None,
    )
    if _cjk is not None:
        plt.rcParams["font.family"] = [_cjk, "sans-serif"]
    plt.rcParams["axes.unicode_minus"] = False

    panels = {
        "S_perigee": "S族 — 近地点",
        "S_perilune": "S族 — 近月点",
        "N_perigee": "N族 — 近地点",
        "N_perilune": "N族 — 近月点",
    }
    points = {
        name: clouds.get(name) or [] for name in POINT_CLOUD_NAMES
    }
    i0_all = [p["i0_deg"] for pts in points.values() for p in pts]
    if not i0_all:
        raise ValueError("plot_iau_moon_point_clouds: no point-cloud records")
    norm = Normalize(vmin=min(i0_all), vmax=max(i0_all))
    cmap = plt.get_cmap("viridis")

    def _fmt_lon(v: float, _pos: float) -> str:
        """Tick label: wrap across the 0/360 seam, keep negatives negative."""
        w = float(v)
        if w >= 360.0 - 1e-6:
            w -= 360.0
            if abs(w) < 1e-6:
                w = 0.0
        return f"{w:g}"

    fig, axes = plt.subplots(2, 2, figsize=(13.0, 8.6), dpi=dpi,
                             layout="constrained")
    # Reserve the top band for the (CJK) suptitle + subtitle so the
    # constrained layout does not push the panel titles into them.
    fig.get_layout_engine().set(rect=(0.0, 0.0, 1.0, 0.90))
    axes = np.asarray(axes).ravel()
    for ax, name in zip(axes, POINT_CLOUD_NAMES):
        pts = points[name]
        for branch, marker in (("asc", "o"), ("desc", "s")):
            sub = [p for p in pts if p["family"].lower().endswith(branch)]
            if not sub:
                continue
            # Draw the ±360° wraparound copies so families crossing the
            # 0°/360° seam (the perigee arcs) read as continuous.
            lon = []
            lat = []
            col = []
            for p in sub:
                for shift in (-360.0, 0.0, 360.0):
                    lon.append(p["lon_deg"] + shift)
                    lat.append(p["lat_deg"])
                    col.append(p["i0_deg"])
            ax.scatter(
                lon, lat, c=col, cmap=cmap, norm=norm,
                marker=marker, s=14, alpha=0.85, linewidths=0.2,
                edgecolors="white", label=f"{branch}（{marker}）",
            )
        n_asc = sum(1 for p in pts if p["family"].lower().endswith("asc"))
        n_desc = len(pts) - n_asc
        ax.set_title(
            f"{panels[name]}   ({n_asc} asc + {n_desc} desc)", fontsize=12,
        )
        box = (panel_boxes or {}).get(name)
        lim = (axis_limits or {}).get(name)
        if box is not None:
            # Reference bounding box (e.g. paper figures 18-22): pad the
            # window so overflow stays visible and draw the box dashed.
            x_lo, x_hi, y_lo, y_hi = (float(v) for v in box)
            x_pad = max(1e-6, (x_hi - x_lo) * box_pad_frac)
            y_pad = max(1e-6, (y_hi - y_lo) * box_pad_frac)
            ax.set_xlim(x_lo - x_pad, x_hi + x_pad)
            ax.set_ylim(y_lo - y_pad, y_hi + y_pad)
            ax.add_patch(plt.Rectangle(
                (x_lo, y_lo), x_hi - x_lo, y_hi - y_lo,
                fill=False, edgecolor="#d62728", linewidth=1.2,
                linestyle="--", zorder=1, label="参考边界 (paper)",
            ))
        elif lim is not None:
            # Caller-pinned window (e.g. fixed axes for an animation).
            x_lo, x_hi, y_lo, y_hi = (float(v) for v in lim)
            ax.set_xlim(x_lo, x_hi)
            ax.set_ylim(y_lo, y_hi)
        else:
            # Adaptive per-panel window: circular-arc longitude (small arc
            # even when the points straddle the seam) + padded latitude.
            x_lo, x_hi = _circular_lon_limits([p["lon_deg"] for p in pts])
            y_lo, y_hi = _padded_limits([p["lat_deg"] for p in pts])
            ax.set_xlim(x_lo, x_hi)
            ax.set_ylim(y_lo, y_hi)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
        ax.xaxis.set_major_formatter(FuncFormatter(_fmt_lon))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
        ax.grid(alpha=0.25, linewidth=0.5)
        if y_lo <= 0.0 <= y_hi:
            ax.axhline(0.0, color="#888888", linewidth=0.7, linestyle="--")
        ax.tick_params(labelsize=9)
        ax.legend(loc="best", fontsize=8, framealpha=0.85)

    for ax in axes[2:]:
        ax.set_xlabel(f"{frame_label} 经度 (deg)", fontsize=10)
    for ax in axes[::2]:
        ax.set_ylabel(f"{frame_label} 纬度 (deg)", fontsize=10)

    fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=list(axes),
        shrink=0.85, pad=0.01, label="停车轨道倾角 i0 (deg)",
    )
    fig.suptitle(
        title or f"R-1S 族表点云 — {frame_label}", fontsize=14,
        fontweight="bold", y=0.985,
    )
    if subtitle:
        fig.text(0.5, 0.945, subtitle, ha="center", fontsize=10, color="#444444")

    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out


def _dup_field_diffs(r1: dict, r2: dict) -> dict[str, float]:
    """Component-wise differences between two family records.

    Angles (perigee / perilune lon/lat) are compared as great-circle
    separations so 360° wraparound and poles are handled correctly.
    """

    def sph_sep(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
        v1 = np.array([
            math.cos(math.radians(lat1)) * math.cos(math.radians(lon1)),
            math.cos(math.radians(lat1)) * math.sin(math.radians(lon1)),
            math.sin(math.radians(lat1)),
        ])
        v2 = np.array([
            math.cos(math.radians(lat2)) * math.cos(math.radians(lon2)),
            math.cos(math.radians(lat2)) * math.sin(math.radians(lon2)),
            math.sin(math.radians(lat2)),
        ])
        # atan2(|v1 x v2|, v1 . v2) is stable at zero separation, unlike
        # acos(v1 . v2) whose float noise floor is ~1e-6 deg.
        return math.degrees(
            math.atan2(float(np.linalg.norm(np.cross(v1, v2))), float(np.dot(v1, v2)))
        )

    return {
        "perigee_deg": sph_sep(
            float(r1["perigee_lat_deg"]), float(r1["perigee_lon_deg"]),
            float(r2["perigee_lat_deg"]), float(r2["perigee_lon_deg"]),
        ),
        "perilune_deg": sph_sep(
            float(r1["perilune_lat_deg"]), float(r1["perilune_lon_deg"]),
            float(r2["perilune_lat_deg"]), float(r2["perilune_lon_deg"]),
        ),
        "vp_km_s": abs(float(r1["vp_km_s"]) - float(r2["vp_km_s"])),
        "v_perilune_km_s": abs(
            float(r1["v_perilune_km_s"]) - float(r2["v_perilune_km_s"])
        ),
        "C3": abs(float(r1["C3"]) - float(r2["C3"])),
        "i_m_deg": abs(float(r1["i_m_deg"]) - float(r2["i_m_deg"])),
    }


def find_duplicate_pairs(
    families: dict[str, list[dict]],
    tol: float = DUP_TOL,
) -> list[dict]:
    """Find complete duplicates among all family records.

    Two records are a complete duplicate pair when every downstream field
    agrees within ``tol`` (deg / km/s).  The search is global: it covers
    records within one family and across families, including pairs with
    different ``i0`` (the i0 <-> 180° − i0 mirror is only an approximate
    symmetry, so it does not match this criterion).
    """
    flat: list[tuple[str, dict]] = [
        (fam, r) for fam, recs in families.items() for r in recs
    ]
    order = {name: k for k, name in enumerate(FAMILY_SEEDS.keys())}
    flat.sort(key=lambda fr: (fr[1]["i0_deg"], order.get(fr[0], 99)))

    pairs: list[dict] = []
    n = len(flat)
    for i in range(n):
        fam1, r1 = flat[i]
        for j in range(i + 1, n):
            fam2, r2 = flat[j]
            diffs = _dup_field_diffs(r1, r2)
            if all(v <= tol for v in diffs.values()):
                pairs.append({
                    "a": {"family": fam1, "i0_deg": float(r1["i0_deg"])},
                    "b": {"family": fam2, "i0_deg": float(r2["i0_deg"])},
                    "diffs": diffs,
                })
    return pairs


def dedupe_families(
    families: dict[str, list[dict]],
    tol: float = DUP_TOL,
) -> tuple[dict[str, list[dict]], dict]:
    """Remove complete duplicates from the family records.

    Deterministic policy: process records in ascending ``i0`` order
    (family order of ``FAMILY_SEEDS`` breaks ties) and drop a record only
    when it duplicates an earlier (therefore lower-``i0``) record.  This
    keeps the prograde-side representative that the launch-site search
    (i0 in [1°, 90°]) can actually select.  A record is never dropped if
    it would empty its family.

    Returns ``(deduped_families, report)`` where ``report`` carries the
    removed pairs and counts.
    """
    order = {name: k for k, name in enumerate(FAMILY_SEEDS.keys())}
    flat: list[tuple[str, dict]] = [
        (fam, r) for fam, recs in families.items() for r in recs
    ]
    flat.sort(key=lambda fr: (fr[1]["i0_deg"], order.get(fr[0], 99)))

    n = len(flat)
    dropped: set[int] = set()
    removed: list[dict] = []
    kept_count: dict[str, int] = {fam: 0 for fam in families}
    for _fam, _r in flat:
        kept_count[_fam] += 1
    for i in range(n):
        if i in dropped:
            continue
        fam1, r1 = flat[i]
        for j in range(i + 1, n):
            if j in dropped:
                continue
            fam2, r2 = flat[j]
            # Keep at least one record per family (never empty a family).
            if kept_count.get(fam2, 0) <= 1:
                continue
            diffs = _dup_field_diffs(r1, r2)
            if all(v <= tol for v in diffs.values()):
                dropped.add(j)
                kept_count[fam2] -= 1
                removed.append({
                    "kept": {"family": fam1, "i0_deg": float(r1["i0_deg"])},
                    "removed": {
                        "family": fam2, "i0_deg": float(r2["i0_deg"]),
                    },
                    "diffs": diffs,
                })

    out: dict[str, list[dict]] = {fam: [] for fam in families}
    for k, (fam, r) in enumerate(flat):
        if k not in dropped:
            out[fam].append(r)
    for fam in out:
        out[fam].sort(key=lambda r: r["i0_deg"])

    report = {
        "tol": tol,
        "n_removed": len(removed),
        "removed": removed,
    }
    return out, report


def assemble_scan_result(
    families: dict[str, list[dict]],
    frame: mf.MEPFrame,
    *,
    i_m_deg: float,
    moon_model: str,
    T_h: float,
    H_m: float,
    polar_ref: str,
    im_step_deg: float,
    i0_min: float,
    i0_max: float,
    i0_step: float,
    dedupe: bool = True,
    include_point_clouds: bool = True,
    point_cloud_plot: Path | None = None,
    point_cloud_plot_dpi: int = 200,
    tol: float = DUP_TOL,
    t_arrive_utc: str = T_ARRIVE_UTC,
) -> dict:
    """Assemble the scan JSON payload from solved family records.

    Applies the complete-duplicate dedup (if ``dedupe``), computes the
    reachable bands, and attaches the four IAU_MOON point clouds (if
    ``include_point_clouds``).  Split out of :func:`run_scan` so callers
    that solve the four families in parallel can share the same assembly.

    ``t_arrive_utc`` is the MEP-frame-freezing arrival epoch (ISO-8601,
    matching ``frame.t_arrive``); it is recorded in ``calibration`` and on
    the point clouds.

    When ``point_cloud_plot`` is given, the four clouds are additionally
    rendered to that PNG via :func:`plot_iau_moon_point_clouds` (requires
    ``include_point_clouds=True``).
    """
    if point_cloud_plot is not None and not include_point_clouds:
        raise ValueError(
            "point_cloud_plot requires include_point_clouds=True"
        )
    dedupe_report = {"enabled": bool(dedupe), "tol": tol, "n_removed": 0, "removed": []}
    if dedupe:
        families, dedupe_report = dedupe_families(families, tol=tol)
        dedupe_report["enabled"] = True

    all_perigee_lat = [
        r["perigee_lat_deg"] for recs in families.values() for r in recs
    ]
    all_perilune_lat = [
        r["perilune_lat_deg"] for recs in families.values() for r in recs
    ]

    doc = {
        "calibration": {
            "h_p_km": H_P,
            "H_m_km": H_m,
            "T_h": T_h,
            "i_m_deg": i_m_deg,
            "moon_model": moon_model,
            "t_arrive_utc": t_arrive_utc,
            "polar_ref": polar_ref,
            "i0_range_deg": [i0_min, i0_max],
            "i0_step_deg": i0_step,
            "im_step_deg": im_step_deg,
        },
        "families": families,
        "reachable_band": {
            "perigee_lat_min_deg": float(min(all_perigee_lat)),
            "perigee_lat_max_deg": float(max(all_perigee_lat)),
            "perigee_lat_span_deg": float(max(all_perigee_lat) - min(all_perigee_lat)),
            "perilune_lat_min_deg": float(min(all_perilune_lat)),
            "perilune_lat_max_deg": float(max(all_perilune_lat)),
        },
        "dedupe": dedupe_report,
    }
    if include_point_clouds:
        doc["iau_moon_point_clouds"] = {
            "frame": "IAU_MOON",
            "epoch_utc": t_arrive_utc,
            "note": (
                "MEP lon/lat directions rotated to IAU_MOON via the frozen "
                "MEP frame and SPICE pxform('J2000','IAU_MOON',t_arrive); "
                "perigee = Earth-side direction, perilune = Moon-side."
            ),
            **build_iau_moon_point_clouds(families, frame),
        }
        if point_cloud_plot is not None:
            plot_iau_moon_point_clouds(
                doc["iau_moon_point_clouds"], point_cloud_plot,
                subtitle=(
                    f"T = {T_h:g} h, i_m = {i_m_deg:g}°, H_m = {H_m:g} km, "
                    f"i0 = {i0_min:g}–{i0_max:g}° (step {i0_step:g}°), "
                    f"moon={moon_model}, sun perturbation, "
                    f"epoch {t_arrive_utc}"
                ),
                dpi=point_cloud_plot_dpi,
            )
    return doc


def run_scan(
    i_m_deg: float,
    moon_model: str = "circular",
    T_h: float = T_H,
    family: str | None = None,
    H_m: float = DEFAULT_H_M_KM,
    polar_ref: str = "mep",
    im_step_deg: float = IM_STEP_DEG,
    seed_table: Path | None = None,
    use_sun: bool = False,
    i0_min: float = I0_MIN,
    i0_max: float = I0_MAX,
    i0_step: float = I0_STEP,
    dedupe: bool = True,
    include_point_clouds: bool = True,
    point_cloud_plot: Path | None = None,
    point_cloud_plot_dpi: int = 200,
    t_arrive_utc: str = T_ARRIVE_UTC,
) -> dict:
    """Build the arbitrary-i_m family table (same schema as r1s_scan.json).

    ``t_arrive_utc`` freezes the MEP frame (arrival at perilune, ISO-8601
    UTC).  It defaults to the module constant :data:`T_ARRIVE_UTC`; pass a
    different epoch to build the family table for that arrival geometry.

    When ``seed_table`` is given, the i0 = 5° polar seeds are loaded from
    that table (which must be at the target ``T_h`` and i_m = 90°) instead
    of continuing the Table-4 seeds from 72 h — much faster.

    ``i0_min`` / ``i0_max`` / ``i0_step`` set the scan range of the parking
    orbit inclination relative to the MEP fundamental plane.  The defaults
    are the legacy ``r1s_scan`` range [1°, 114°] in 1° steps; passing
    ``i0_max=179`` extends the table to the full inclination circle
    (retrograde members included).

    Complete duplicates (records whose full downstream geometry matches
    within ``scan_im.DUP_TOL``) are dropped before returning; the dropped
    pairs are listed in ``result["dedupe"]``.  The four MEP->IAU_MOON
    perigee/perilune lon/lat point clouds are attached under
    ``result["iau_moon_point_clouds"]``, and (optionally) rendered to
    ``point_cloud_plot`` as a 2x2 PNG.
    """
    t_arrive = g.utc_to_et(t_arrive_utc)
    frame = mf.mep_at_epoch(t_arrive)
    mf.assert_moon_at_anti_x(frame)

    families: dict[str, list[dict]] = {}
    fam_items = [family] if family else list(FAMILY_SEEDS.keys())
    for fam in fam_items:
        families[fam] = scan_family_im(
            fam, FAMILY_SEEDS[fam], t_arrive, frame, i_m_deg, moon_model,
            T_h=T_h, H_m=H_m, polar_ref=polar_ref, im_step_deg=im_step_deg,
            seed_table=seed_table, use_sun=use_sun,
            i0_min=i0_min, i0_max=i0_max, i0_step=i0_step,
        )

    return assemble_scan_result(
        families, frame,
        i_m_deg=i_m_deg, moon_model=moon_model, T_h=T_h, H_m=H_m,
        polar_ref=polar_ref, im_step_deg=im_step_deg,
        i0_min=i0_min, i0_max=i0_max, i0_step=i0_step,
        dedupe=dedupe, include_point_clouds=include_point_clouds,
        point_cloud_plot=point_cloud_plot,
        point_cloud_plot_dpi=point_cloud_plot_dpi,
        t_arrive_utc=t_arrive_utc,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Arbitrary-i_m R-1S family table builder"
    )
    parser.add_argument(
        "--i-m", type=float, required=True,
        help="target lunar-orbit inclination (deg)",
    )
    parser.add_argument("--out", type=Path, required=True, help="output JSON path")
    parser.add_argument(
        "--moon-model", default="circular", choices=["circular", "spice"]
    )
    parser.add_argument("--transfer-hours", type=float, default=T_H)
    parser.add_argument("--family", default=None, choices=list(FAMILY_SEEDS.keys()))
    parser.add_argument(
        "--perilune-alt", type=float, default=DEFAULT_H_M_KM,
        help=f"perilune altitude H_m (km), default {DEFAULT_H_M_KM}",
    )
    parser.add_argument(
        "--polar-ref", default="mep", choices=("mep", "moon_equator"),
        help=(
            "q3 reference axis: 'mep' (default, frozen white-pole) or "
            "'moon_equator' (IAU_MOON Z, true lunar polar direction)."
        ),
    )
    parser.add_argument("--im-step", type=float, default=IM_STEP_DEG)
    parser.add_argument(
        "--t-arrive", type=str, default=T_ARRIVE_UTC,
        help=(
            f"arrival epoch (ISO-8601 UTC) freezing the MEP frame "
            f"(default {T_ARRIVE_UTC})"
        ),
    )
    parser.add_argument("--i0-min", type=float, default=I0_MIN,
                        help="minimum parking-orbit MEP inclination i0 (deg)")
    parser.add_argument("--i0-max", type=float, default=I0_MAX,
                        help="maximum parking-orbit MEP inclination i0 (deg)")
    parser.add_argument("--i0-step", type=float, default=I0_STEP,
                        help="i0 scan step (deg)")
    parser.add_argument(
        "--seed-table", type=Path, default=None,
        help=(
            "existing polar (i_m=90°) family table at the target T_h to seed "
            "the i0=5° solutions (avoids the slow 72h->T_h continuation)"
        ),
    )
    parser.add_argument(
        "--use-sun", action="store_true",
        help="include solar third-body perturbation in the dynamics",
    )
    parser.add_argument(
        "--point-cloud-plot", type=Path, default=None,
        help=(
            "optional PNG path: render the four MEP->IAU_MOON perigee/"
            "perilune point clouds as a 2x2 lon/lat figure"
        ),
    )
    parser.add_argument(
        "--point-cloud-plot-dpi", type=int, default=200,
        help="DPI for --point-cloud-plot (default 200)",
    )
    args = parser.parse_args()

    result = run_scan(
        i_m_deg=args.i_m,
        moon_model=args.moon_model,
        T_h=args.transfer_hours,
        family=args.family,
        H_m=args.perilune_alt,
        polar_ref=args.polar_ref,
        im_step_deg=args.im_step,
        seed_table=args.seed_table,
        use_sun=args.use_sun,
        i0_min=args.i0_min,
        i0_max=args.i0_max,
        i0_step=args.i0_step,
        point_cloud_plot=args.point_cloud_plot,
        point_cloud_plot_dpi=args.point_cloud_plot_dpi,
        t_arrive_utc=args.t_arrive,
    )
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)

    band = result["reachable_band"]
    first_fam = next(iter(result["families"]))
    print(
        f"scan_im: i_m={args.i_m} deg, {len(result['families'])} families x "
        f"{len(result['families'][first_fam])} i0 steps"
    )
    print(
        f"  reachable perigee-lat band: [{band['perigee_lat_min_deg']:.3f}, "
        f"{band['perigee_lat_max_deg']:.3f}] deg"
    )
    print(f"  wrote {args.out}")
    if args.point_cloud_plot is not None:
        print(f"  wrote point-cloud figure {args.point_cloud_plot}")


if __name__ == "__main__":
    main()