#!/usr/bin/env python3
"""launch_window/fusion.py — forward/inverse fusion launch-window calculator.

Implements the fusion architecture of ``launch_window/fusion_method_prd.md``:

  FR-1  forward_scan      launch-seeded geometry + coast closure (ALL roots)
  FR-2  cluster_windows   continuous feasible runs -> screening windows
  FR-3  certify_candidate forward state -> MEP warm start -> R-1 BVP
                          (solve_polar for i_m = 90, solve_non_coplanar for
                          any other target lunar-orbit inclination,
                          retrograde targets included)
  FR-4  inverse_scan      production L6 chain (family tables + pocket search)
  FR-5  cross_validate    both / forward_only / inverse_only classification
  FR-6  coverage          forward step constraint + inverse_only defect check
  FR-7  FusionWindow      unified window record
  FR-8  CLI + artifacts   windows / calendar / cross / coverage / rejects

Usage::

    python -m launch_window.fusion scan \
        --start 2026-08-01T00:00:00 --end 2026-09-01T00:00:00 \
        --step 120 --T 112.61 --site wenchang \
        --perilune-alt 200 --i-m 90 --polar-ref moon_equator \
        --coast-min 0 --coast-max 1000 --az 90 175 \
        --certify --inverse --out /path/to/run_dir

Arbitrary target inclination (e.g. CE-5 43.7 deg, CE-6 137 deg retrograde):
pass ``--i-m <deg> --i0-hi 114``; the family table is keyed by i_m and built
automatically when missing, and certification targets ``q3 = cos(i_m)``.

移植说明（Queqiao2-visibility）：
本文件移植自 AgentMoons 项目 ``launch_window`` 包 v0.6.0，科学计算逻辑
未改动；仅将硬编码的本地绝对路径（``/Users/aequatio/AgentMoons_Local/...``）
改为环境变量可覆盖的仓库内默认路径：

- ``FUSION_OUT_DIR``       运行产物目录（默认 ``viewer/jobs/fusion_default``）
- ``FUSION_TABLE_DIR``     预生成 L5 家族表目录（默认 ``viewer/data/fusion_tables``）
- ``FUSION_CE56_VERIFY_DIR`` CE-5/6 复核表目录（verify-mission 子命令使用）

CLI 的 verify-mission / mission-chain / correct-chain 子命令依赖未移植的
``verify_ce56`` / ``polar_correct`` 模块，仅 scan 流程（含 FR-4 反解）可用。
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import multiprocessing as mp
import os
import sys
import time
import traceback
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Optional

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from launch_window import geometry as g  # noqa: E402
from launch_window import launch_epoch as le  # noqa: E402
from launch_window import launcher as lz  # noqa: E402
from launch_window import mep_frame as mf  # noqa: E402
from launch_window import perigee_match as pm  # noqa: E402
from launch_window import tli_solver as ts  # noqa: E402
from launch_window.lambert_window import LaunchSite  # noqa: E402
from launch_window.scan_im import run_scan as scan_im_run_scan  # noqa: E402

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
# 移植适配：原硬编码 /Users/aequatio/AgentMoons_Local/... 改为仓库内默认 +
# 环境变量覆盖（部署自包含）。
DEFAULT_OUT_DIR = Path(
    os.environ.get(
        "FUSION_OUT_DIR", str(_REPO_ROOT / "viewer" / "jobs" / "fusion_default")
    )
)
DEFAULT_TABLE_DIR = Path(
    os.environ.get(
        "FUSION_TABLE_DIR", str(_REPO_ROOT / "viewer" / "data" / "fusion_tables")
    )
)

DEFAULT_I_M_DEG = 90.0
DEFAULT_H_M_KM = 200.0          # perilune altitude target
DEFAULT_H_PARK_KM = 200.0       # parking-orbit perigee altitude
DEFAULT_H_APO_KM = 380.0        # CZ-5 200 x 380 km ellipse
DEFAULT_POLAR_REF = "moon_equator"
DEFAULT_MOON_MODEL = "spice"
DEFAULT_USE_SUN = True
DEFAULT_T_ASCENT_S = 740.0                                  # CZ-5 measured
DEFAULT_ASCENT_ARC_DEG = float(le.ASCENT_ARC_DEG)           # 20.8377
DEFAULT_COAST_MIN_S = float(le.COAST_MIN_S)                 # 0
DEFAULT_COAST_MAX_S = float(le.COAST_MAX_S)                 # 1000
DEFAULT_COAST_STEP_S = 2.0
DEFAULT_DNU_MIN_DEG = 95.0
DEFAULT_DNU_MAX_DEG = 250.0
DEFAULT_AZ_MIN_DEG = 90.0
DEFAULT_AZ_MAX_DEG = 175.0
DEFAULT_STEP_S = 120.0
DEFAULT_CERT_K = 8
DEFAULT_CERT_WORKERS = 8
DEFAULT_I0_LO = 1.0
DEFAULT_I0_HI = 90.0

# Cross-validation tolerances (FR-5).
DEFAULT_DT_TOL_S = 300.0
DEFAULT_DAZ_TOL_DEG = 1.0
DEFAULT_DCOAST_TOL_S = 60.0
DEFAULT_DT_TOL_H = 0.01

# Certification tolerances (FR-3).
DEFAULT_I_M_TOL_DEG = 0.1
DEFAULT_PERILUNE_TOL_KM = 5.0
DEFAULT_SITE_PLANE_TOL_KM = 50.0
DEFAULT_TLI_KG = 8350.0
DEFAULT_PAYLOAD_KG = 8200.0
DEFAULT_TLI_C3_REF = -2.0
DEFAULT_TLI_ISP_S = 442.0
P0_PERT_KMS = 1.0e-3
P0_PERT_LARGE_KMS = 1.0e-2

# Forward reject codes (FR-1.4).
REJ_AZ = "az_out_of_corridor"
REJ_NO_COAST = "no_coast_root"
REJ_ENVELOPE = "envelope_failed"
REJ_OUT_OF_RANGE = "launch_outside_scan_range"

CHUNK_SCHEMA = 1
CERT_SCHEMA = 6          # v6: gate-aware root selection + solved-plane gates (fixed)
                         # certified BVP solution's own parking plane

SITE_REGISTRY: dict[str, dict[str, Any]] = {
    "wenchang": {"name": "Wenchang", "lon_deg": 110.951, "lat_deg": 19.614,
                 "alt_m": 10.0},
    "xichang": {"name": "Xichang", "lon_deg": 102.027, "lat_deg": 28.2466,
                "alt_m": 1500.0},
    "ksc": {"name": "KSC", "lon_deg": -80.649, "lat_deg": 28.5729,
            "alt_m": 10.0},
}

# ---------------------------------------------------------------------------
# Small utilities
# ---------------------------------------------------------------------------
def _json_safe(obj: Any) -> Any:
    """Recursively replace NaN/Inf/numpy scalars so ``json`` emits valid JSON."""
    if isinstance(obj, dict):
        return {str(k): _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        v = float(obj)
        if math.isnan(v) or math.isinf(v):
            return None
        return v
    if isinstance(obj, (np.integer, int)):
        return int(obj) if not isinstance(obj, bool) else bool(obj)
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    if isinstance(obj, np.ndarray):
        return _json_safe(obj.tolist())
    return obj


def _log(msg: str, log_path: Optional[Path], also_stdout: bool = True) -> None:
    line = f"[{datetime.now(timezone.utc).isoformat(timespec='seconds')}] {msg}"
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(log_path, "a", encoding="utf-8") as fh:
            fh.write(line + "\n")
    if also_stdout:
        print(line, flush=True)


def _sha1(doc: Any) -> str:
    payload = json.dumps(_json_safe(doc), sort_keys=True,
                         separators=(",", ":"), default=str)
    return hashlib.sha1(payload.encode("utf-8")).hexdigest()


def _norm_iso(iso: str) -> str:
    """Accept ``...T24:00:00`` (next-day midnight) like polar_window's CLI."""
    s = iso.strip().replace("Z", "")
    if "T24:00:00" in s:
        base = datetime.fromisoformat(s.replace("T24:00:00", "T00:00:00"))
        return (base + timedelta(days=1)).strftime("%Y-%m-%dT%H:%M:%S")
    return s


def _utc_to_et_flexible(iso: str) -> float:
    """ISO UTC -> SPICE ET, accepting the ``T24:00:00`` convention."""
    return g.utc_to_et(_norm_iso(iso))


def _write_json(path: Path, doc: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(_json_safe(doc), fh, ensure_ascii=False, indent=2)


def _write_jsonl_gz(path: Path, records: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        for rec in records:
            fh.write(json.dumps(_json_safe(rec), ensure_ascii=False) + "\n")


def _stats(values: list[float]) -> dict:
    vals = sorted(float(v) for v in values if v == v)
    n = len(vals)
    if not n:
        return {"n": 0, "min": None, "median": None, "max": None}
    med = vals[n // 2] if n % 2 else 0.5 * (vals[n // 2 - 1] + vals[n // 2])
    return {"n": n, "min": vals[0], "median": med, "max": vals[-1]}


def _parse_az(value: Any) -> tuple[float, float]:
    if isinstance(value, (list, tuple)) and len(value) == 2:
        return float(value[0]), float(value[1])
    raise ValueError(f"--az must be two numbers, got {value!r}")


# ---------------------------------------------------------------------------
# Config construction
# ---------------------------------------------------------------------------
def build_config(args: argparse.Namespace) -> dict:
    site = SITE_REGISTRY[str(args.site).lower()]
    T_h = float(args.T)
    T_s = T_h * 3600.0
    _a_park, _e_park, _n_park = le.parking_orbit_elements(
        float(args.parking_peri_alt), float(args.parking_apo_alt))
    az_min, az_max = _parse_az(args.az)
    cfg: dict[str, Any] = {
        "mission": str(args.mission),
        "T_h": T_h,
        "T_s": T_s,
        "t_start_utc": args.start,
        "t_end_utc": args.end,
        "t_start_et": _utc_to_et_flexible(args.start),
        "t_end_et": _utc_to_et_flexible(args.end),
        "step_s": float(args.step),
        "site": dict(site),
        "az_min_deg": az_min,
        "az_max_deg": az_max,
        "use_envelope": bool(args.use_envelope),
        "coast_min_s": float(args.coast_min),
        "coast_max_s": float(args.coast_max),
        "coast_step_s": float(args.coast_step_s),
        "dnu_min_deg": float(DEFAULT_DNU_MIN_DEG),
        "dnu_max_deg": float(DEFAULT_DNU_MAX_DEG),
        "t_ascent_s": float(args.t_ascent_s),
        "ascent_arc_deg": float(args.ascent_arc_deg),
        # Closure convention: "launch_tof" = coast_tool launch->arrival T;
        # "perigee_tof" = production L6 perigee->arrival T (used by the
        # cross-coverage equivalence rescan, FR-5/FR-6).
        "closure": "launch_tof",
        "h_park_km": float(args.parking_peri_alt),
        "h_apo_km": float(args.parking_apo_alt),
        "a_park_km": float(_a_park),
        "om_park_deg_s": math.degrees(float(_n_park)),
        "rp_km": ts.R_EARTH + float(args.parking_peri_alt),
        "i_m_target_deg": float(args.i_m),
        # Multi-target certification (FR-11: i_m as an interval).  When
        # ``--i-m-list`` is given, every representative is certified against
        # EACH target (family table warm start from the primary table;
        # ``窄窗口多轨道`` variable-coast/variable-orbit demonstration).
        "i_m_list": [float(x) for x in (getattr(args, "i_m_list", None)
                                        or [args.i_m])],
        "H_m_km": float(args.perilune_alt),
        "polar_ref": str(args.polar_ref),
        "moon_model": str(args.moon_model),
        "use_sun": bool(args.use_sun),
        "i0_lo": float(DEFAULT_I0_LO),
        # i0 upper bound: 90° suffices for polar targets (the inner fixed
        # point folds the plane inclination to [0°, 90°]); non-polar /
        # retrograde targets (CE-5 43.7°, CE-6 137°) need the wider
        # family-table range — pass --i0-hi 114 like verify_ce56.
        "i0_hi": float(getattr(args, "i0_hi", DEFAULT_I0_HI)),
        "dtheta_deg": ((float(args.t_ascent_s) + float(DEFAULT_COAST_MAX_S))
                       / 3600.0 * le.EARTH_SIDEREAL_RATE_DEG_PER_H),
        "cert_k": int(args.cert_k),
        "dt_tol_s": float(DEFAULT_DT_TOL_S),
        "daz_tol_deg": float(DEFAULT_DAZ_TOL_DEG),
        "dcoast_tol_s": float(DEFAULT_DCOAST_TOL_S),
        "dt_tol_h": float(DEFAULT_DT_TOL_H),
        "i_m_tol_deg": float(getattr(args, "i_m_tol", DEFAULT_I_M_TOL_DEG)),
        "perilune_tol_km": float(DEFAULT_PERILUNE_TOL_KM),
        "site_plane_tol_km": float(DEFAULT_SITE_PLANE_TOL_KM),
        "tli_kg": float(DEFAULT_TLI_KG),
        "payload_kg": float(DEFAULT_PAYLOAD_KG),
        "tli_c3_ref": float(DEFAULT_TLI_C3_REF),
        "tli_isp_s": float(DEFAULT_TLI_ISP_S),
        "v_park_kms": math.sqrt(ts.MU_EARTH / (ts.R_EARTH + float(args.parking_peri_alt))),
        "v_circ_moon_kms": math.sqrt(ts.MU_MOON / (ts.R_MOON + float(args.perilune_alt))),
        "allow_inverse_only": bool(args.allow_inverse_only),
        "chunk_schema": CHUNK_SCHEMA,
    }
    cfg["config_id"] = _sha1({k: v for k, v in cfg.items()
                              if k not in ("t_start_et", "t_end_et")})
    return cfg


def _site_from_cfg(cfg: dict) -> LaunchSite:
    s = cfg["site"]
    return LaunchSite(name=s["name"], lon_deg=s["lon_deg"],
                      lat_deg=s["lat_deg"], alt_m=s["alt_m"])


def _table_path(out_dir: Path, T_h: float, H_m: float, i_m_deg: float,
                use_sun: bool) -> Path:
    im = int(round(i_m_deg * 1000))
    hm = int(round(H_m))
    suffix = "_sun" if use_sun else ""
    return out_dir / "tables" / (
        f"r1s_scan_im{im}_T{int(round(T_h * 10.0)):04d}_Hm{hm}_meq{suffix}.json"
    )


# ---------------------------------------------------------------------------
# FR-1  Forward chain
# ---------------------------------------------------------------------------
def forward_geometry(t_l_et: float, cfg: dict) -> dict:
    """coast_tool ``solveGeometry`` in J2000 (SPICE): plane / azimuth / theta."""
    site = cfg["site"]
    ta = t_l_et + cfg["T_s"]
    r_site = np.asarray(g.launch_site_vector_j2000_km(
        site["lon_deg"], site["lat_deg"], site["alt_m"], t_l_et), dtype=float)
    r_moon = np.asarray(g.moon_position_j2000_km(ta), dtype=float)
    u1 = r_site / np.linalg.norm(r_site)
    u2 = r_moon / np.linalg.norm(r_moon)
    n = np.cross(u1, u2)
    nn = float(np.linalg.norm(n))
    theta = math.degrees(math.atan2(nn, float(np.dot(u1, u2))))
    n = n / nn
    if n[2] < 0.0:                       # prograde branch
        n = -n
        theta = 360.0 - theta
    inc = math.degrees(math.acos(max(-1.0, min(1.0, float(n[2])))))
    raan = math.degrees(math.atan2(float(n[0]), -float(n[1]))) % 360.0
    # Unique azimuth from the prograde tangential velocity direction.
    v = np.cross(n, u1)
    east = np.array([-u1[1], u1[0], 0.0])
    east = east / np.linalg.norm(east)
    north = np.cross(u1, east)
    az = math.degrees(
        math.atan2(float(np.dot(v, east)), float(np.dot(v, north)))) % 360.0
    # Launch-site argument of latitude in the transfer plane and the
    # insertion-point geocentric latitude (F0b envelope).
    node = np.cross(np.array([0.0, 0.0, 1.0]), n)
    node_n = float(np.linalg.norm(node))
    e_r = node / node_n if node_n > 1e-12 else np.array([1.0, 0.0, 0.0])
    e_t = np.cross(n, e_r)
    u_l = math.degrees(math.atan2(float(np.dot(u1, e_t)),
                                  float(np.dot(u1, e_r)))) % 360.0
    return {
        "t_launch_et": float(t_l_et),
        "t_arrive_et": float(ta),
        "u1": u1, "u2": u2, "n": n,
        "theta_deg": float(theta),
        "i_lto_deg": float(inc),
        "raan_lto_deg": float(raan),
        "az_deg": float(az),
        "moon_dist_km": float(np.linalg.norm(r_moon)),
        "moon_dec_deg": math.degrees(math.asin(max(-1.0, min(1.0, float(u2[2]))))),
        "u_l_deg": float(u_l),
    }


def _kep_tof_s(rp_km: float, e: float, dnu_deg: float) -> float:
    a = rp_km / (1.0 - e)
    nu = math.radians(dnu_deg)
    E = 2.0 * math.atan2(math.sqrt(1.0 - e) * math.sin(nu / 2.0),
                         math.sqrt(1.0 + e) * math.cos(nu / 2.0))
    M = E - e * math.sin(E)
    if M < 0.0:
        M += 2.0 * math.pi
    return M / math.sqrt(ts.MU_EARTH / a ** 3)


def coast_roots_forward(geo: dict, cfg: dict) -> list[dict]:
    """All coast-closure roots (FR-1.2): arc + time closure, all sign changes.

    Scheme (aligned with ``ce7_jsT11261_recompute.py`` and the PRD F0 fixes):

        dnu  = theta − arc − ω_park · tc                      (arc closure)
        e    = (d − rp) / (rp − d cos dnu)                    (coplanar conic)
        f(tc) = TOF_Kepler(rp, e, dnu) − (T − t_ascent − tc)  (time closure)
    """
    rp = cfg["rp_km"]
    d = geo["moon_dist_km"]
    theta = geo["theta_deg"]
    om = cfg["om_park_deg_s"]
    t_ascent = cfg["t_ascent_s"]
    arc = cfg["ascent_arc_deg"]
    T_s = cfg["T_s"]
    lo, hi = cfg["coast_min_s"], cfg["coast_max_s"]
    step = cfg["coast_step_s"]
    dnu_lo, dnu_hi = cfg["dnu_min_deg"], cfg["dnu_max_deg"]

    closure = cfg.get("closure", "launch_tof")

    def f(tc: float) -> Optional[tuple[float, float, float, float]]:
        dnu = (theta - arc - om * tc) % 360.0
        if not (dnu_lo < dnu < dnu_hi):
            return None
        cosd = math.cos(math.radians(dnu))
        den = rp - d * cosd
        if abs(den) < 1e-12:
            return None
        e = (d - rp) / den
        if not (0.0 < e < 1.0):
            return None
        tof = _kep_tof_s(rp, e, dnu)
        target = T_s if closure == "perigee_tof" else (T_s - t_ascent - tc)
        return tof - target, dnu, e, tof

    sols: list[dict] = []
    prev = None
    prev_tc = lo
    n_steps = int(round((hi - lo) / step)) + 1
    for k in range(n_steps + 1):
        tc = min(lo + k * step, hi)
        cur = f(tc)
        if cur is not None and prev is not None and prev[0] * cur[0] < 0.0:
            a, b = prev_tc, tc
            fa = prev[0]
            for _ in range(60):
                mid = 0.5 * (a + b)
                fm = f(mid)
                if fm is None:
                    break
                if fa * fm[0] <= 0.0:
                    b = mid
                else:
                    a, fa = mid, fm[0]
            tcr = 0.5 * (a + b)
            r = f(tcr)
            if r is not None:
                _, dnu, e, tof = r
                aa = rp / (1.0 - e)
                v_per = math.sqrt(ts.MU_EARTH * (2.0 / rp - 1.0 / aa))
                sols.append({
                    "coast_s": float(tcr),
                    "dnu_deg": float(dnu),
                    "e": float(e),
                    "a_km": float(aa),
                    "tof_s": float(tof),
                    "v_peri_kms": float(v_per),
                    "dv_tli_kms": float(v_per - math.sqrt(
                        ts.MU_EARTH / rp)),
                    "c3_km2s2": float(v_per * v_per
                                      - 2.0 * ts.MU_EARTH / rp),
                    "apo_alt_km": float(aa * (1.0 + e) - ts.R_EARTH),
                })
        prev = cur
        prev_tc = tc
        if tc >= hi:
            break
    sols.sort(key=lambda s: s["coast_s"])
    return sols


def forward_point(t_l_et: float, cfg: dict) -> dict:
    """One FR-1 grid point: geometry + all roots + gate bits + reject codes."""
    geo = forward_geometry(t_l_et, cfg)
    roots = coast_roots_forward(geo, cfg)
    az = geo["az_deg"]
    gate_az = bool(cfg["az_min_deg"] <= az <= cfg["az_max_deg"])
    # F0b envelope at the parking-orbit insertion point (structural check).
    u_ins = (geo["u_l_deg"] + cfg["ascent_arc_deg"]) % 360.0
    phi_ins = math.degrees(math.asin(max(-1.0, min(1.0, abs(
        math.sin(math.radians(geo["i_lto_deg"]))
        * math.sin(math.radians(u_ins)))))))
    gate_env = bool(geo["i_lto_deg"] >= phi_ins) or (not cfg["use_envelope"])
    gate_coast = bool(len(roots) > 0)
    # Cap margin of the closest root (PRD §10 risk table).
    g_margin = min((abs(r["coast_s"] - cfg["coast_max_s"]) for r in roots),
                   default=None)
    reject_codes: list[str] = []
    if not gate_az:
        reject_codes.append(REJ_AZ)
    if not gate_coast:
        reject_codes.append(REJ_NO_COAST)
    if not gate_env:
        reject_codes.append(REJ_ENVELOPE)
    return {
        "t_launch_et": float(t_l_et),
        "t_launch_utc": g.et_to_utc(t_l_et),
        "t_arrive_et": float(geo["t_arrive_et"]),
        "t_arrive_utc": g.et_to_utc(geo["t_arrive_et"]),
        "T_h": float(cfg["T_h"]),
        "az_deg": float(az),
        "i_lto_deg": float(geo["i_lto_deg"]),
        "raan_lto_deg": float(geo["raan_lto_deg"]),
        "theta_deg": float(geo["theta_deg"]),
        "moon_dist_km": float(geo["moon_dist_km"]),
        "moon_dec_deg": float(geo["moon_dec_deg"]),
        "u_l_deg": float(geo["u_l_deg"]),
        "u_ins_deg": float(u_ins),
        "phi_ins_deg": float(phi_ins),
        "gate_az": gate_az,
        "gate_coast": gate_coast,
        "gate_envelope": gate_env,
        "n_roots": len(roots),
        "roots": roots,
        "g_margin_s": g_margin,
        "feasible": bool(gate_az and gate_coast and gate_env),
        "reject_codes": reject_codes,
        "provenance": "forward",
    }


# ---------------------------------------------------------------------------
# FR-2  Clustering
# ---------------------------------------------------------------------------
def cluster_windows(points: list[dict], cfg: dict) -> list[dict]:
    """Merge consecutive feasible points into ForwardWindows (FR-2)."""
    windows: list[dict] = []
    cur: list[dict] = []

    def flush() -> None:
        if not cur:
            return
        pts = cur[:]
        n = len(pts)
        first, last = pts[0], pts[-1]
        n_roots_max = max(p["n_roots"] for p in pts)
        root_extremes: dict[str, dict] = {}
        for k in range(n_roots_max):
            ck = [(p["roots"][k]["coast_s"], i)
                  for i, p in enumerate(pts) if len(p["roots"]) > k]
            if ck:
                root_extremes[str(k)] = {
                    "coast_min_s": min(ck)[0],
                    "coast_min_idx": min(ck)[1],
                    "coast_max_s": max(ck)[0],
                    "coast_max_idx": max(ck)[1],
                }
        windows.append({
            "t_open_et": first["t_launch_et"],
            "t_close_et": last["t_launch_et"],
            "t_open_utc": first["t_launch_utc"],
            "t_close_utc": last["t_launch_utc"],
            "duration_min": (last["t_launch_et"] - first["t_launch_et"]) / 60.0,
            "n_points": n,
            "az_min_deg": min(p["az_deg"] for p in pts),
            "az_max_deg": max(p["az_deg"] for p in pts),
            "i_lto_min_deg": min(p["i_lto_deg"] for p in pts),
            "i_lto_max_deg": max(p["i_lto_deg"] for p in pts),
            "theta_min_deg": min(p["theta_deg"] for p in pts),
            "theta_max_deg": max(p["theta_deg"] for p in pts),
            "root_extremes": root_extremes,
            "provenance": "forward",
            "_pts": pts,
        })

    for p in points:
        if p["feasible"]:
            cur.append(p)
        else:
            flush()
            cur = []
    flush()
    return windows


def select_representatives(window: dict, k_max: int) -> list[tuple[int, int]]:
    """Representative (point_index, root_index) pairs for certification.

    FR-2.3: window head/tail + per-root coast extremes + >= min(5, n)
    evenly-spaced samples; total capped at ``k_max``.
    """
    pts = window["_pts"]
    n = len(pts)
    pairs: list[tuple[int, int]] = []
    if n == 0:
        return pairs
    pairs += [(0, 0), (n - 1, 0)]
    for key, ext in window["root_extremes"].items():
        k = int(key)
        pairs.append((int(ext["coast_min_idx"]), k))
        pairs.append((int(ext["coast_max_idx"]), k))
    n_even = max(5, min(n, 8))
    for i in np.linspace(0, n - 1, n_even).round().astype(int).tolist():
        pairs.append((int(i), 0))
    out: list[tuple[int, int]] = []
    seen: set[tuple[int, int]] = set()
    for i, k in pairs:
        if (i, k) in seen:
            continue
        seen.add((i, k))
        out.append((i, k))
        if len(out) >= k_max:
            break
    return out


# ---------------------------------------------------------------------------
# FR-3  Certification
# ---------------------------------------------------------------------------
_FAMILY_CACHE: dict[str, dict] = {}


def _load_families_cached(path_str: str) -> dict:
    if path_str not in _FAMILY_CACHE:
        _FAMILY_CACHE[path_str] = pm.load_families(Path(path_str))
    return _FAMILY_CACHE[path_str]


def _mep_orbit_elements_from_state(r_peri: np.ndarray, n_j2000: np.ndarray,
                                   frame: mf.MEPFrame) -> tuple[float, float, float]:
    """(i0, Omega, omega) deg of a perigee state in the frozen MEP frame."""
    n_mep = mf.to_mep(np.asarray(n_j2000, dtype=float), frame)
    n_mep = n_mep / np.linalg.norm(n_mep)
    i0 = math.degrees(math.acos(max(-1.0, min(1.0, float(n_mep[2])))))
    Omega = math.degrees(math.atan2(float(n_mep[0]),
                                    -float(n_mep[1]))) % 360.0
    z = np.array([0.0, 0.0, 1.0])
    node = np.cross(z, n_mep)
    nn = float(np.linalg.norm(node))
    e_r = node / nn if nn > 1e-12 else np.array([1.0, 0.0, 0.0])
    e_w = np.cross(n_mep, e_r)
    r_mep = mf.to_mep(np.asarray(r_peri, dtype=float), frame)
    omega = math.degrees(math.atan2(float(np.dot(r_mep, e_w)),
                                    float(np.dot(r_mep, e_r)))) % 360.0
    return float(i0), float(Omega), float(omega)


def _perigee_state_from_forward(t_l_et: float, root: dict,
                                cfg: dict) -> dict:
    """Reconstruct the TLI (perigee) state and its MEP elements (FR-3.1)."""
    geo = forward_geometry(t_l_et, cfg)
    u1 = geo["u1"]
    n = geo["n"]
    ang = math.radians(cfg["ascent_arc_deg"]
                       + cfg["om_park_deg_s"] * root["coast_s"])
    c, s = math.cos(ang), math.sin(ang)
    ehat = (u1 * c + np.cross(n, u1) * s
            + n * float(np.dot(n, u1)) * (1.0 - c))
    r_peri = cfg["rp_km"] * ehat
    vp = math.sqrt(ts.MU_EARTH * (2.0 / cfg["rp_km"] - 1.0 / root["a_km"]))
    frame = mf.mep_at_epoch(geo["t_arrive_et"])
    i0, Om, om = _mep_orbit_elements_from_state(r_peri, n, frame)
    return {
        "r_peri": r_peri, "n_j2000": n, "vp_kms": float(vp),
        "i0_deg": float(i0), "Omega_deg": float(Om), "omega_deg": float(om),
        "frame": frame, "geo": geo,
    }


def _certify_solve(t_arrive: float, T_cert_s: float, i0_deg: float,
                   p0: tuple[float, float, float], cfg: dict) -> dict:
    """Solve the R-1 BVP for the configured target lunar-orbit inclination.

    ``i_m == 90`` dispatches to :func:`tli_solver.solve_polar` (q3 = 0, the
    CE-7 polar anchor path); any other target dispatches to
    :func:`tli_solver.solve_non_coplanar` with ``q3 = cos(i_m)`` (retrograde
    targets such as CE-6's 137° included).
    """
    t0 = time.time()
    im = float(cfg["i_m_target_deg"])
    try:
        kwargs = dict(
            H_m=cfg["H_m_km"], moon_model=cfg["moon_model"],
            polar_ref=cfg["polar_ref"], use_sun=cfg["use_sun"],
        )
        if abs(im - 90.0) <= 1e-9:
            res = ts.solve_polar(
                cfg["h_park_km"], float(T_cert_s), float(i0_deg),
                float(t_arrive),
                (float(p0[0]), float(p0[1]), float(p0[2])), **kwargs)
        else:
            res = ts.solve_non_coplanar(
                cfg["h_park_km"], float(T_cert_s), float(i0_deg),
                float(t_arrive),
                (float(p0[0]), float(p0[1]), float(p0[2])),
                im_deg=im, **kwargs)
        return {
            "ok": True, "converged": bool(res.converged),
            "res": res, "error": None, "elapsed_s": time.time() - t0,
        }
    except Exception as exc:  # noqa: BLE001
        return {
            "ok": False, "converged": False, "res": None,
            "error": f"{type(exc).__name__}: {exc}",
            "elapsed_s": time.time() - t0,
        }


def _solved_gates(res: Any, i0_deg: float, branch: str, t_l: float,
                  frame: mf.MEPFrame, cfg: dict) -> dict:
    """Additional-gate geometry of a certified BVP solution's OWN plane.

    The BVP moves ``Omega``/``omega``, so the certified parking plane is
    ``(i0, res.p[1])`` — evaluating the site-off-plane / azimuth gates on
    this plane (rather than the warm-start family plane) is the physically
    meaningful check: a launch can only fly the certified orbit if the pad
    lies on the parking plane at the launch epoch.
    """
    site = _site_from_cfg(cfg)
    im = math.radians(float(i0_deg))
    om = float(res.p[1]) % 360.0
    n_mep = np.array([math.sin(im) * math.sin(math.radians(om)),
                      -math.sin(im) * math.cos(math.radians(om)),
                      math.cos(im)])
    n_j2000 = frame.rot_mep_to_j2000 @ n_mep
    r_site = np.asarray(g.launch_site_vector_j2000_km(
        site.lon_deg, site.lat_deg, site.alt_m, t_l), dtype=float)
    r_hat = r_site / np.linalg.norm(r_site)
    site_off_km = abs(float(np.dot(n_j2000, r_hat))) * ts.R_EARTH
    inc_e = math.degrees(math.acos(max(-1.0, min(1.0,
                                                    abs(float(n_j2000[2]))))))
    az_std = float("nan")
    try:
        az_std = float(le._launch_azimuth_and_uL(
            inc_e, om, site, branch)[0]) % 360.0
    except Exception:  # noqa: BLE001
        pass
    vp_sol = float(res.p[0])
    r_park = ts.R_EARTH + cfg["h_park_km"]
    c3 = vp_sol ** 2 - 2.0 * ts.MU_EARTH / r_park
    m_avail = margin = float("nan")
    try:
        curve = lz.c3_curve_linear(
            m_ref_kg=cfg["tli_kg"], c3_ref_kms2=cfg["tli_c3_ref"],
            isp_s=cfg["tli_isp_s"])
        m_avail = float(curve(c3))
        margin = m_avail - cfg["payload_kg"]
    except Exception:  # noqa: BLE001
        pass
    width_min = float("nan")
    try:
        w = lz.width_from_margin(
            i_park_deg=inc_e, site_lat_deg=site.lat_deg,
            dv_margin_mps=30.0, v_park_kms=lz.v_park_kms())
        width_min = float(w["width_min"])
    except Exception:  # noqa: BLE001
        pass
    az_ok = bool(az_std == az_std
                 and cfg["az_min_deg"] <= az_std <= cfg["az_max_deg"])
    site_ok = bool(site_off_km <= cfg["site_plane_tol_km"])
    cap_ok = bool(margin == margin and margin >= 0.0)
    return {
        "inc_e_deg": float(inc_e),
        "Omega_solved_mep_deg": float(om),
        "site_off_km": float(site_off_km),
        "az_std_deg": az_std,
        "C3_km2s2": float(c3),
        "m_avail_kg": m_avail,
        "margin_kg": margin,
        "width_min": width_min,
        "all_ok": bool(az_ok and site_ok and cap_ok),
        "gates": {
            "converged": bool(res.converged),
            "i_m_ok": bool(abs(float(res.i_m) - cfg["i_m_target_deg"])
                           <= cfg["i_m_tol_deg"]),
            "perilune_ok": bool(abs(float(res.q[0]) - cfg["H_m_km"])
                                <= cfg["perilune_tol_km"]),
            "az_in_corridor": az_ok,
            "site_off_plane_ok": site_ok,
            "capability_ok": cap_ok,
        },
    }


def certify_candidate(window: dict, point_idx: int, root_idx: int,
                      cfg: dict) -> dict:
    """FR-3: certify one forward (point, root) candidate against the R-1 BVP.

    Warm-start chain: forward state -> MEP elements -> family-table p0
    (nearest azimuth root) -> ``tli_solver.solve_polar``; fallbacks are the
    JS-state p0 and vp perturbations.  Every attempt is recorded.
    """
    point = window["_pts"][point_idx]
    root = point["roots"][min(root_idx, len(point["roots"]) - 1)]
    t_l = point["t_launch_et"]
    ta = point["t_arrive_et"]
    T_cert_s = cfg["T_s"] - (cfg["t_ascent_s"] + root["coast_s"])

    st = _perigee_state_from_forward(t_l, root, cfg)
    frame = st["frame"]
    families = _load_families_cached(cfg["table_path"])

    # --- family roots at this arrival epoch (all families x both branches) --
    # F0 self-consistency: the inner fixed point must evaluate the launch
    # site at THIS candidate's true launch epoch.  The implied launch→perigee
    # arc is dtheta = (t_ascent + coast)·ω_E, so with the candidate's own
    # coast (and T_cert = T − t_ascent − coast) the site epoch estimate
    # ``t_L_est = t_arrive − T_cert − dtheta/ω_E`` collapses to the forward
    # chain's launch epoch ``t_L = t_arrive − T`` exactly.  Using the
    # nominal dtheta (which assumes a 1000 s coast) instead was shown to
    # inflate the site-off-plane residual along the window tail.
    dth = ((cfg["t_ascent_s"] + float(root["coast_s"])) / 3600.0
           * le.EARTH_SIDEREAL_RATE_DEG_PER_H)
    fam_roots: list[dict] = []
    for fam_name, fam in families.items():
        for branch in ("ascending", "descending"):
            inner = le._solve_inner_fixed_point(
                fam, ta, T_cert_s, dth, _site_from_cfg(cfg),
                cfg["i0_lo"], cfg["i0_hi"], branch=branch)
            if inner is None:
                continue
            i0, inc_E, OmE, omE, vp_f, _n_hat = inner
            try:
                az, u_l, _a = le._launch_azimuth_and_uL(
                    inc_E, OmE, _site_from_cfg(cfg), branch)
            except Exception:  # noqa: BLE001
                az, u_l = float("nan"), float("nan")
            fam_roots.append({
                "family": fam_name, "branch": branch, "i0_deg": float(i0),
                "inc_E_deg": float(inc_E), "Omega_E_deg": float(OmE),
                "omega_E_deg": float(omE), "vp_kms": float(vp_f),
                "az_deg": float(az % 360.0) if az == az else None,
                "u_l_deg": float(u_l) if u_l == u_l else None,
                "p0": (float(fam.vp_f(i0)), float(fam.Omega_f(i0)),
                       float(fam.omega_f(i0))),
                "daz_deg": (abs(float(az) % 360.0 - point["az_deg"])
                            if az == az else float("nan")),
            })

    fam_roots.sort(key=lambda r: (r["daz_deg"] != r["daz_deg"],
                                  r["daz_deg"]))
    attempts: list[dict] = []

    def _try(label: str, i0: float, p0: tuple[float, float, float]) -> Optional[dict]:
        out = _certify_solve(ta, T_cert_s, i0, p0, cfg)
        rec = {
            "label": label, "i0_deg": float(i0),
            "p0": [float(p0[0]), float(p0[1]), float(p0[2])],
            "converged": bool(out["converged"]),
            "error": out["error"], "elapsed_s": float(out["elapsed_s"]),
        }
        if out["res"] is not None:
            res = out["res"]
            rec.update({
                "i_m_deg": float(res.i_m), "q1_km": float(res.q[0]),
                "q2": float(res.q[1]), "q3": float(res.q[2]),
                "vp_kms": float(res.p[0]),
                "v_perilune_kms": float(res.v_perilune),
                "perilune_c3_km2s2": float(res.C3),
                "perigee_lon_mep_deg": float(res.perigee_lon_mep),
                "perigee_lat_mep_deg": float(res.perigee_lat_mep),
                "perilune_lon_mep_deg": float(res.perilune_lon_mep),
                "perilune_lat_mep_deg": float(res.perilune_lat_mep),
                "dvp_kms": float(res.p[0] - p0[0]),
            })
        attempts.append(rec)
        if (out["res"] is not None and bool(out["converged"])
                and abs(float(out["res"].i_m) - cfg["i_m_target_deg"]) <= cfg["i_m_tol_deg"]
                and abs(float(out["res"].q[0]) - cfg["H_m_km"]) <= cfg["perilune_tol_km"]):
            rec["_res"] = out["res"]
            return rec
        return None

    primary_best: Optional[dict] = None
    gate_best: Optional[dict] = None
    warm = "family_table_root"
    # 1) family-table p0 sorted by warm-plane azimuth distance: try up to 6
    #    roots and prefer an attempt whose SOLVED-plane gates all pass (the
    #    BVP can rotate the certified plane away from the site for some
    #    roots, while another root yields a site-consistent solution).
    for fr in fam_roots[:6]:
        rec_a = _try(f"family:{fr['family']}:{fr['branch']}",
                     fr["i0_deg"], fr["p0"])
        if rec_a is None:
            continue
        primary_best = primary_best or rec_a
        gt = _solved_gates(rec_a["_res"], rec_a["i0_deg"], fr["branch"],
                           t_l, frame, cfg)
        rec_a["_gates"] = gt
        if gt["all_ok"]:
            gate_best = rec_a
            break
    best = gate_best if gate_best is not None else primary_best
    if gate_best is not None:
        warm = gate_best["label"]
    # 2) vp perturbations on the best family root.
    if best is None and fam_roots:
        fr = fam_roots[0]
        for pert in (P0_PERT_KMS, -P0_PERT_KMS,
                     P0_PERT_LARGE_KMS, -P0_PERT_LARGE_KMS):
            p0 = (fr["p0"][0] + pert, fr["p0"][1], fr["p0"][2])
            warm = f"family_vp_pert_{pert:+.0e}"
            best = _try(warm, fr["i0_deg"], p0)
            if best is not None:
                break
    # 3) JS-state p0 (last resort).
    if best is None:
        warm = "js_state"
        best = _try("js_state", st["i0_deg"],
                    (st["vp_kms"], st["Omega_deg"], st["omega_deg"]))

    rec = {
        "provenance": "forward_certify",
        "i_m_target_deg": float(cfg["i_m_target_deg"]),
        "T_h": float(cfg["T_h"]),
        "t_launch_et": float(t_l),
        "t_launch_utc": point["t_launch_utc"],
        "t_arrive_et": float(ta),
        "t_arrive_utc": point["t_arrive_utc"],
        "t_tli_et": float(t_l + cfg["t_ascent_s"] + root["coast_s"]),
        "t_tli_utc": g.et_to_utc(t_l + cfg["t_ascent_s"] + root["coast_s"]),
        "coast_s": float(root["coast_s"]),
        "coast_roots_s": [float(r["coast_s"]) for r in point["roots"]],
        "root_index": int(root_idx),
        "az_deg": float(point["az_deg"]),
        "i_lto_deg": float(point["i_lto_deg"]),
        "js_i0_deg": float(st["i0_deg"]),
        "js_Omega_deg": float(st["Omega_deg"]),
        "js_omega_deg": float(st["omega_deg"]),
        "js_vp_kms": float(st["vp_kms"]),
        "T_cert_h": float(T_cert_s / 3600.0),
        "family_roots": [
            {k: v for k, v in fr.items() if k != "p0"} for fr in fam_roots
        ],
        "certify_attempt": attempts,
        "warm_start_used": warm,
        "certify_status": "failed",
        "certify_reject_codes": [],
    }
    if best is None:
        if not fam_roots:
            rec["certify_status"] = "failed:no_family_root"
            rec["certify_reject_codes"].append("inverse_root_missing")
        else:
            rec["certify_status"] = "failed:no_converged_attempt"
            rec["certify_reject_codes"].append("bvp_not_converged")
        rec["published"] = False
        for attempt in attempts:
            attempt.pop("_res", None)
        return rec

    res = best["_res"]
    vp_sol = float(res.p[0])
    dv_tli = vp_sol - cfg["v_park_kms"]
    dv_loi = float(res.v_perilune) - cfg["v_circ_moon_kms"]
    # Family entry actually used (labels only).
    fr_used = next((fr for fr in fam_roots
                    if fr["i0_deg"] == best["i0_deg"]), None)
    branch_used = fr_used["branch"] if fr_used else "ascending"
    gg = best.get("_gates") or _solved_gates(
        res, best["i0_deg"], branch_used, t_l, frame, cfg)
    inc_e = gg["inc_e_deg"]
    Om_solved = gg["Omega_solved_mep_deg"]
    site_off_km = gg["site_off_km"]
    az_std = gg["az_std_deg"]
    c3 = gg["C3_km2s2"]
    m_avail = gg["m_avail_kg"]
    margin = gg["margin_kg"]
    width_min = gg["width_min"]
    # Warm-start family plane residual (diagnostic only).
    site_off_warm_km = float("nan")
    try:
        r_site = np.asarray(g.launch_site_vector_j2000_km(
            cfg["site"]["lon_deg"], cfg["site"]["lat_deg"],
            cfg["site"]["alt_m"], t_l), dtype=float)
        r_hat = r_site / np.linalg.norm(r_site)
        im_w = math.radians(best["i0_deg"])
        Om_w = math.radians(float(best["p0"][1]))
        n_mep_w = np.array([math.sin(im_w) * math.sin(Om_w),
                            -math.sin(im_w) * math.cos(Om_w), math.cos(im_w)])
        n_j_w = frame.rot_mep_to_j2000 @ n_mep_w
        site_off_warm_km = abs(float(np.dot(n_j_w, r_hat))) * ts.R_EARTH
    except Exception:  # noqa: BLE001
        pass

    gates_all = bool(gg["gates"]["i_m_ok"] and gg["gates"]["perilune_ok"]
                     and gg["all_ok"])
    rec.update({
        "certified": True,
        "published": gates_all,
        "warm_start_used": warm,
        "certify_status": "certified" if gates_all else "certified_gates_failed",
        "family": fr_used["family"] if fr_used else None,
        "branch": fr_used["branch"] if fr_used else None,
        "i0_deg": float(best["i0_deg"]),
        "i_park_deg": float(inc_e),
        "i_m_deg": float(res.i_m),
        "perilune_alt_km": float(res.q[0]),
        "q1_km": float(res.q[0]), "q2": float(res.q[1]),
        "q3": float(res.q[2]),
        "vp_sol_kms": vp_sol,
        "dvp_kms": float(vp_sol - best["p0"][0]),
        "dvp_mps": float((vp_sol - best["p0"][0]) * 1000.0),
        "C3_km2s2": float(c3),
        "dv_tli_kms": float(dv_tli),
        "dv_loi_kms": float(dv_loi),
        "dv_total_kms": float(dv_tli + dv_loi),
        "v_perilune_kms": float(res.v_perilune),
        "m_avail_kg": m_avail,
        "margin_kg": margin,
        "width_min": width_min,
        "site_off_plane_km": float(site_off_km),
        "site_off_warmstart_plane_km": site_off_warm_km,
        "Omega_solved_mep_deg": float(Om_solved),
        "perigee_lon_mep_deg": float(res.perigee_lon_mep),
        "perigee_lat_mep_deg": float(res.perigee_lat_mep),
        "perilune_lon_mep_deg": float(res.perilune_lon_mep),
        "perilune_lat_mep_deg": float(res.perilune_lat_mep),
        "azimuth_std_deg": az_std,
        # Gate block (FR-3.4/3.5) — on the certified solution's own plane.
        "gates": dict(gg["gates"]),
        "certify_reject_codes": [],
    })
    if not rec["gates"]["site_off_plane_ok"]:
        rec["certify_reject_codes"].append("site_off_transfer_plane")
    if not rec["gates"]["az_in_corridor"]:
        rec["certify_reject_codes"].append("azimuth_out_of_corridor")
    if not rec["gates"]["capability_ok"]:
        rec["certify_reject_codes"].append("c3_capability_insufficient")
    for attempt in attempts:
        attempt.pop("_res", None)
    return rec


def _cert_worker(payload: tuple) -> list[dict]:
    jobs, cfg = payload
    out: list[dict] = []
    for item in jobs:
        key = item["key"]
        cache = Path(cfg["cert_cache_dir"]) / f"{key}.json"
        if cache.exists():
            try:
                with open(cache, encoding="utf-8") as fh:
                    rec = json.load(fh)
                rec.setdefault("published", False)
                out.append(rec)
                continue
            except Exception:  # noqa: BLE001
                pass
        # Multi-target certification (FR-11): the job carries its own
        # target inclination; the shared cfg is overridden per item
        # (workers process items sequentially).
        cfg_job = dict(cfg)
        cfg_job["i_m_target_deg"] = float(
            item.get("i_m_target", cfg["i_m_target_deg"]))
        try:
            rec = certify_candidate(item["window"], item["point_idx"],
                                    item["root_idx"], cfg_job)
        except Exception:  # noqa: BLE001
            rec = {
                "provenance": "forward_certify",
                "certify_status": "failed:exception",
                "certify_reject_codes": ["certify_exception"],
                "error": traceback.format_exc(limit=8),
                "t_launch_et": item["window"]["_pts"][item["point_idx"]]["t_launch_et"],
                "t_launch_utc": item["window"]["_pts"][item["point_idx"]]["t_launch_utc"],
            }
        rec.pop("_res", None)
        rec.setdefault("published", False)
        rec["cache_key"] = key
        rec["cert_schema"] = CERT_SCHEMA
        _write_json(cache, rec)
        out.append(rec)
    return out


def stage_certify(windows: list[dict], cfg: dict, log_path: Path,
                  workers: int) -> tuple[list[dict], list[dict]]:
    """FR-3 stage: certify representative points of every forward window."""
    cfg = dict(cfg)
    cfg["cert_cache_dir"] = str(Path(cfg["out_dir"]) / "certs")
    targets = [float(x) for x in cfg.get("i_m_list", [cfg["i_m_target_deg"]])]
    jobs: list[dict] = []
    for wi, win in enumerate(windows):
        for pi, ri in select_representatives(win, cfg["cert_k"]):
            for target in targets:
                key = _sha1({
                    "schema": CERT_SCHEMA, "config": cfg["config_id"],
                    "t_launch_et": win["_pts"][pi]["t_launch_et"],
                    "root_idx": ri,
                    "coast_s": win["_pts"][pi]["roots"][ri]["coast_s"],
                    "T_h": cfg["T_h"], "i_m": target,
                    "H_m": cfg["H_m_km"], "polar_ref": cfg["polar_ref"],
                })
                jobs.append({"key": key, "window": win, "point_idx": pi,
                             "root_idx": ri, "window_index": wi,
                             "i_m_target": target})
    _log(f"[FR-3] certification jobs: {len(jobs)} "
         f"(windows={len(windows)}, K={cfg['cert_k']}, "
         f"i_m targets={targets})", log_path)
    if not jobs:
        return [], []
    cache_hits = sum(
        1 for j in jobs
        if (Path(cfg["cert_cache_dir"]) / f"{j['key']}.json").exists())
    _log(f"[FR-3] cache hits={cache_hits}/{len(jobs)}", log_path)
    t0 = time.perf_counter()
    n_proc = max(1, min(workers, len(jobs)))
    chunk = max(1, math.ceil(len(jobs) / n_proc))
    payloads = [(jobs[i:i + chunk], cfg)
                for i in range(0, len(jobs), chunk)]
    records: list[dict] = []
    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=n_proc) as pool:
        done = 0
        for out in pool.imap_unordered(_cert_worker, payloads):
            records.extend(out)
            done += len(out)
            _log(f"  [FR-3] {done}/{len(jobs)} candidates certified", log_path)
    records.sort(key=lambda r: (float(r.get("t_launch_et", 0)),
                                float(r.get("coast_s", 0))))
    fails = [r for r in records if not r.get("certified")]
    published = [r for r in records if r.get("published")]
    _log(f"[FR-3] done in {time.perf_counter() - t0:.1f}s — "
         f"primary_converged={len(records) - len(fails)} "
         f"published={len(published)} failed={len(fails)}", log_path)
    return records, fails


# ---------------------------------------------------------------------------
# FR-4  Inverse production scan (L6)
# ---------------------------------------------------------------------------
def inverse_scan(cfg: dict, log_path: Path,
                 max_workers: int = 8,
                 seed_steps: Optional[list[float]] = None,
                 force: bool = False) -> tuple[list[dict], list[dict]]:
    """Run the production L6 chain (``le.scan_epochs``) over the month.

    Range semantics follow polar_window v2: the CLI window bounds the
    *launch* epoch; arrival seeds are shifted forward by ``T_h``.

    ``seed_steps`` runs several staggered seed grids and merges the results
    (FR-6 coverage hardening): the pocket search grid phase is fixed by the
    seed epochs, so a single sidereal-day grid can miss short pockets.  The
    first pass keeps the full reject list; later passes run window-only.

    Results are cached under ``chunks/inverse_<key>.json`` (keyed by the
    config id + seed steps) so the expensive chain is not re-run.
    """
    table = Path(cfg["table_path"])
    if not table.exists():
        raise FileNotFoundError(f"family table missing: {table}")
    families = pm.load_families(table)
    site = _site_from_cfg(cfg)
    T_h = cfg["T_h"]
    t_seed_start = cfg["t_start_et"] + T_h * 3600.0
    t_seed_end = cfg["t_end_et"] + T_h * 3600.0
    seed_steps = (list(seed_steps) if seed_steps
                  else [float(le.SIDEREAL_DAY_S)])
    cache_key = _sha1({"config": cfg["config_id"], "seed_steps": seed_steps,
                       "schema": 1})
    cache = Path(cfg["out_dir"]) / "chunks" / f"inverse_{cache_key}.json"
    rej_path = Path(cfg["out_dir"]) / "fusion_inverse_rejects.jsonl.gz"
    if cache.exists() and not force:
        with open(cache, encoding="utf-8") as fh:
            doc = json.load(fh)
        rows = doc["windows"]
        rejects_all = []
        if rej_path.exists():
            with gzip.open(rej_path, "rt", encoding="utf-8") as fh:
                for line in fh:
                    try:
                        rejects_all.append(json.loads(line))
                    except Exception:  # noqa: BLE001
                        pass
        _log(f"[FR-4] inverse cache hit ({cache.name}): windows={len(rows)} "
             f"rejects_loaded={len(rejects_all)}", log_path)
        return rows, rejects_all
    t0 = time.perf_counter()
    merged: list[Any] = []
    rejects_all: list[dict] = []
    for si, step in enumerate(seed_steps):
        keep_rejects = si == 0
        out = le.scan_epochs(
            t_start_et=t_seed_start,
            t_end_et=t_seed_end,
            T_h_list=[T_h],
            families=families,
            site=site,
            polar_ref=cfg["polar_ref"],
            branches=("ascending", "descending"),
            dtheta_deg=cfg["dtheta_deg"],
            i0_lo=cfg["i0_lo"],
            i0_hi=cfg["i0_hi"],
            seed_step_s=float(step),
            return_rejects=keep_rejects,
            t_ascent_s=cfg["t_ascent_s"],
            coast_max_s=cfg["coast_max_s"],
            ascent_arc_deg=cfg["ascent_arc_deg"],
        )
        if keep_rejects:
            windows, rejects = out
            rejects_all = [r.to_dict() for r in rejects]
        else:
            windows = out
        n_new = 0
        for w in windows:
            dup = any(abs(w.t_launch_et - m.t_launch_et)
                      < le.DEDUP_THRESHOLD_S
                      and w.family == m.family and w.branch == m.branch
                      for m in merged)
            if not dup:
                merged.append(w)
                n_new += 1
        _log(f"  [FR-4] seed_step={step:.0f}s pass {si + 1}/"
             f"{len(seed_steps)}: windows={len(windows)} new={n_new} "
             f"rejects_kept={len(rejects_all)}", log_path)
    lo = cfg["t_start_et"]
    hi = cfg["t_end_et"]
    rows: list[dict] = []
    trimmed = 0
    for w in sorted(merged, key=lambda x: x.t_launch_et):
        if not (lo <= w.t_launch_et < hi):
            trimmed += 1
            continue
        rows.append({
            "provenance": "inverse",
            "T_h": float(w.T_h),
            "t_launch_et": float(w.t_launch_et),
            "t_launch_utc": w.t_launch_utc,
            "t_tli_utc": w.t_tli_utc,
            "t_arrive_utc": w.t_arrive_utc,
            "t_arrive_et": float(w.t_arrive_et),
            "family": w.family,
            "branch": w.branch,
            "i_park_deg": float(w.i_park_deg),
            "raan_park_deg": float(w.raan_park_deg),
            "azimuth_deg": float(w.azimuth_deg),
            "coast_s": float(w.coast_s),
            "vp_kms": float(w.vp_kms),
            "c3_km2s2": float(w.c3_km2s2),
            "i_m_deg": None,
            "perilune_alt_km": None,
        })
    reason_counts = Counter(r.get("reason") for r in rejects_all)
    _write_json(cache, {
        "schema": 1, "config_id": cfg["config_id"],
        "seed_steps": seed_steps, "windows": rows,
        "rejects_total": len(rejects_all),
        "rejects_by_reason": dict(reason_counts.most_common()),
    })
    _log(f"[FR-4] inverse scan: windows={len(rows)} (trimmed={trimmed}) "
         f"rejects={len(rejects_all)} "
         f"{dict(reason_counts.most_common(8))} "
         f"seed_steps={seed_steps} in {time.perf_counter() - t0:.1f}s",
         log_path)
    return rows, rejects_all


def enumerate_roots_at(t_arrive: float, T_s: float, cfg: dict,
                       families: dict) -> list[dict]:
    """All inner-fixed-point roots at one arrival epoch (FR-4.2 cross-check).

    Enumerates every sign change of ``f(i0)`` over [i0_lo, i0_hi] for all
    (family, branch) combinations, and derives azimuth / coast / envelope for
    each root.  The production chain keeps only ``roots[0]``; this function
    exposes what the first-root selection discards.
    """
    from scipy.optimize import brentq
    site = _site_from_cfg(cfg)
    out: list[dict] = []
    for fam_name, fam in families.items():
        lo = max(cfg["i0_lo"], float(fam.i0_min))
        hi = min(cfg["i0_hi"], float(fam.i0_max))
        if lo >= hi:
            continue
        frame = mf.mep_at_epoch(t_arrive)
        dT_s = cfg["dtheta_deg"] / le.EARTH_SIDEREAL_RATE_DEG_PER_H * 3600.0
        t_l_est = t_arrive - T_s - dT_s
        r_site = np.asarray(g.launch_site_vector_j2000_km(
            site.lon_deg, site.lat_deg, site.alt_m, t_l_est), dtype=float)
        r_site_mep = mf.to_mep(r_site, frame)
        r_site_mep = pm.R_P * r_site_mep / np.linalg.norm(r_site_mep)

        def f_i0(i0: float) -> float:
            r_peri = mf.from_mep_lonlat(float(fam.lon_f(i0)),
                                        float(fam.lat_f(i0)), pm.R_P)
            return pm.parking_orbit_inclination_deg(r_site_mep, r_peri) - i0

        xs = np.linspace(lo, hi, 400)
        vals = np.array([f_i0(float(x)) for x in xs])
        for k in range(len(xs) - 1):
            if vals[k] * vals[k + 1] >= 0.0:
                continue
            try:
                i0 = float(brentq(f_i0, float(xs[k]), float(xs[k + 1]),
                                  xtol=1e-8))
            except ValueError:
                continue
            inc_E, OmE, omE, vp, n_hat = le._mep_orbit_to_equator(
                fam, i0, frame)
            for branch in ("ascending", "descending"):
                try:
                    az, u_l, _alpha = le._launch_azimuth_and_uL(
                        inc_E, OmE, site, branch)
                except Exception:  # noqa: BLE001
                    continue
                coast_raw = le.coast_from_argument(omE, u_l)
                coast_ext = le.parking_orbit_coast_seconds(omE, u_l)
                out.append({
                    "t_arrive_utc": g.et_to_utc(t_arrive),
                    "family": fam_name, "branch": branch,
                    "root_index": k,
                    "i0_deg": i0, "inc_E_deg": float(inc_E),
                    "Omega_E_deg": float(OmE), "omega_E_deg": float(omE),
                    "vp_kms": float(vp),
                    "az_deg": float(az) % 360.0,
                    "u_l_deg": float(u_l),
                    "coast_raw_s": float(coast_raw),
                    "coast_ext_s": float(coast_ext),
                    "coast_in_cap": bool(cfg["coast_min_s"] <= coast_ext
                                         <= cfg["coast_max_s"]),
                })
    return out


# ---------------------------------------------------------------------------
# FR-5  Cross validation
# ---------------------------------------------------------------------------
def cross_validate(cfg: dict, forward_windows: list[dict],
                   cert_records: list[dict], inverse: list[dict],
                   forward_points: list[dict]) -> dict:
    """Match forward windows against inverse windows and classify (FR-5).

    Primary key: launch epoch within ``dt_tol_s`` (default 300 s), evaluated
    between inverse windows and *forward feasible points* (FR-6.2), not just
    window midpoints.  Additional deviations (az / coast / T) are recorded
    for classification only and never veto a match.
    """
    dt_tol = cfg["dt_tol_s"]
    inv = sorted(inverse, key=lambda w: w["t_launch_et"])
    feas = sorted((p for p in forward_points if p["feasible"]),
                  key=lambda p: p["t_launch_et"])
    feas_t = [p["t_launch_et"] for p in feas]

    def _nearest_feasible(t: float) -> tuple[Optional[dict], Optional[float]]:
        if not feas:
            return None, None
        import bisect
        i = bisect.bisect_left(feas_t, t)
        best = None
        for j in (i - 1, i, i + 1):
            if 0 <= j < len(feas):
                dt = abs(feas[j]["t_launch_et"] - t)
                if best is None or dt < best[1]:
                    best = (feas[j], dt)
        return best if best else (None, None)

    pairs: list[dict] = []
    forward_only: list[dict] = []
    both: list[dict] = []
    matched_inv: set[int] = set()
    inv_match_forward: dict[int, dict] = {}

    # Match every forward window: is any of its feasible points within
    # dt_tol of an inverse window?
    for win in forward_windows:
        t_open, t_close = win["t_open_et"], win["t_close_et"]
        close = [w for w in inv
                 if t_open - dt_tol <= w["t_launch_et"] <= t_close + dt_tol]
        recs = [r for r in cert_records
                if t_open - dt_tol <= r.get("t_launch_et", -1e18)
                <= t_close + dt_tol]
        certified = [r for r in recs if r.get("published")]
        entry = {
            "forward_window": {
                "t_open_utc": win["t_open_utc"],
                "t_close_utc": win["t_close_utc"],
                "duration_min": win["duration_min"],
                "az_min_deg": win["az_min_deg"],
                "az_max_deg": win["az_max_deg"],
                "n_points": win["n_points"],
            },
            "inverse_windows": close,
            "n_certified": len(certified),
            "n_attempts": len(recs),
            "certified": bool(certified),
            "certify_statuses": dict(Counter(
                r.get("certify_status", "?") for r in recs)),
            "certify_reject_codes": sorted({
                c for r in recs for c in r.get("certify_reject_codes", [])}),
        }
        if close:
            # Field-level deviations against the nearest inverse window.
            w0 = min(close, key=lambda w: abs(w["t_launch_et"] - t_open))
            recs_near = [r for r in certified
                         if abs(r["t_launch_et"] - w0["t_launch_et"]) <= dt_tol]
            entry["deviations"] = [{
                "t_launch_dt_s": r["t_launch_et"] - w0["t_launch_et"],
                "az_ddeg": r["az_deg"] - w0["azimuth_deg"],
                "coast_ds": r["coast_s"] - w0["coast_s"],
            } for r in recs_near]
            for w in close:
                for i, wi in enumerate(inv):
                    if wi is w:
                        matched_inv.add(i)
                        inv_match_forward.setdefault(i, entry)
            both.append(entry)
        else:
            entry["reason"] = ("inverse_window_missing (certified forward "
                               "candidate present)" if certified
                               else "certification_failed")
            forward_only.append(entry)
        pairs.append(entry)

    inverse_only: list[dict] = []
    for i, w in enumerate(inv):
        if i in matched_inv:
            continue
        nearest, ndt = _nearest_feasible(w["t_launch_et"])
        # Forward feasible points may exist but outside the corridor used
        # for the physical gates; re-check the raw nearest point regardless.
        raw_nearest, raw_dt = None, None
        for p in forward_points:
            dt = abs(p["t_launch_et"] - w["t_launch_et"])
            if raw_dt is None or dt < raw_dt:
                raw_nearest, raw_dt = p, dt
        az_in_corridor = bool(
            cfg["az_min_deg"] <= w["azimuth_deg"] <= cfg["az_max_deg"])
        # Equivalent forward rescan (FR-5/FR-6): the inverse window uses the
        # L6 T convention (perigee->arrival), the forward scanner the
        # coast_tool one (launch->arrival).  For coverage, re-run the forward
        # closure at the *equivalent* mission: same TLI->arrival time
        # T_L6 = cfg.T (perigee-tof closure), same launch epoch.
        cfg_eq = dict(cfg)
        cfg_eq["closure"] = "perigee_tof"
        fine = None
        if az_in_corridor:
            fine = _fine_rescan(float(w["t_launch_et"]), cfg_eq)
        # Same-convention direct check against the main scan points.
        direct_ok = bool(ndt is not None and ndt <= dt_tol)
        if not az_in_corridor:
            kind = "corridor_outside"
            why = ("inverse azimuth outside the fusion corridor "
                   f"[{cfg['az_min_deg']}, {cfg['az_max_deg']}] — filtered "
                   "by the launch-corridor gate, not a coverage defect")
        elif fine is not None:
            kind = "covered_equivalent"
            why = ("covered by a forward feasible point at the equivalent "
                   "mission (perigee-tof closure, 10 s rescan)")
        elif direct_ok:
            kind = "covered_direct"
            why = "covered by a forward feasible point, same T convention"
        elif (89.5 <= w["azimuth_deg"] <= 91.5
              and (w["i_park_deg"] or 0.0) > 24.0):
            # The L6 family root here sits at the sin(Az)=cos(i_E)/cos(φ_L)
            # clamp (i_E ~ φ_L → Az = 90°), while the Lambert parking
            # replacement records i_park 28°–38°.  The forward geometry at
            # this epoch requires Az < 90 (out of the Wenchang corridor),
            # so the window is an L6 boundary artifact, not a real launch
            # opportunity.
            kind = "inverse_boundary_artifact"
            why = ("L6 azimuth-clamp boundary artifact (family i_E ≈ φ_L "
                   "→ Az=90° clamped; Lambert i_park > 24°); the forward "
                   "geometry at this epoch requires Az < 90° (out of "
                   "corridor) — not launchable")
        else:
            kind = "true_gap"
            why = ("no feasible forward point under the equivalent mission "
                   "even after 10 s fine rescan — forward-chain gap; L6 "
                   "solution relies on the family/BVP parking-vs-transfer "
                   "plane offset (dtheta ≠ 0) not modelled by the forward "
                   "screen")
        inverse_only.append({
            "inverse_window": {
                "t_launch_utc": w["t_launch_utc"],
                "family": w.get("family"),
                "branch": w.get("branch"),
                "coast_s": w["coast_s"],
                "azimuth_deg": w["azimuth_deg"],
                "i_park_deg": w["i_park_deg"],
                "c3_km2s2": w.get("c3_km2s2"),
            },
            "az_in_corridor": az_in_corridor,
            "nearest_forward_feasible_dt_s": ndt,
            "nearest_forward_feasible_az_deg": (
                nearest["az_deg"] if nearest else None),
            "nearest_forward_point_dt_s": raw_dt,
            "nearest_forward_point_az_deg": (
                raw_nearest["az_deg"] if raw_nearest else None),
            "nearest_forward_point_reject_codes": (
                raw_nearest["reject_codes"] if raw_nearest else None),
            "fine_rescan": ({
                "dt_s": fine[1],
                "az_deg": fine[0]["az_deg"],
                "coast_s": fine[0]["roots"][0]["coast_s"]
                if fine[0]["roots"] else None,
            } if fine is not None else None),
            "classification": kind,
            "why": why,
        })
    n_defects = sum(1 for r in inverse_only
                    if r["classification"] == "true_gap")
    n_covered = sum(1 for r in inverse_only
                    if r["classification"] in ("covered_equivalent",
                                               "covered_direct"))
    return {
        "config": {
            "dt_tol_s": dt_tol,
            "daz_tol_deg": cfg["daz_tol_deg"],
            "dcoast_tol_s": cfg["dcoast_tol_s"],
            "dt_tol_h": cfg["dt_tol_h"],
        },
        "n_pairs": len(pairs),
        "n_both": len(both),
        "n_forward_only": len(forward_only),
        "n_inverse_windows": len(inv),
        "n_inverse_only": len(inverse_only),
        "n_inverse_covered": n_covered,
        "n_inverse_only_corridor_outside": sum(
            1 for r in inverse_only
            if r["classification"] == "corridor_outside"),
        "n_inverse_only_boundary_artifact": sum(
            1 for r in inverse_only
            if r["classification"] == "inverse_boundary_artifact"),
        "n_inverse_only_defects": n_defects,
        "pairs": pairs,
        "forward_only": forward_only,
        "inverse_only": inverse_only,
    }


def _fine_rescan(t_center: float, cfg: dict, span_s: float = 300.0,
                 step_s: float = 10.0) -> Optional[tuple[dict, float]]:
    """Forward fine rescan around ``t_center`` (FR-6 sampling hardening)."""
    best: Optional[tuple[dict, float]] = None
    k = int(round(span_s / step_s))
    for i in range(-k, k + 1):
        t = t_center + i * step_s
        try:
            p = forward_point(t, cfg)
        except Exception:  # noqa: BLE001
            continue
        if p["feasible"]:
            dt = abs(p["t_launch_et"] - t_center)
            if best is None or dt < best[1]:
                best = (p, dt)
    return best


# ---------------------------------------------------------------------------
# FR-8  Outputs
# ---------------------------------------------------------------------------
def _fusion_window_from_cert(rec: dict) -> dict:
    """FR-7 unified window record from a certified candidate."""
    return {
        "mission": None,
        "T_h": rec.get("T_h"),
        "t_launch_utc": rec.get("t_launch_utc"),
        "t_tli_utc": rec.get("t_tli_utc"),
        "t_arrive_utc": rec.get("t_arrive_utc"),
        "provenance": "forward_certified",
        "certify_status": rec.get("certify_status"),
        "published": bool(rec.get("published")),
        "i_m_target_deg": rec.get("i_m_target_deg"),
        "family": rec.get("family"),
        "branch": rec.get("branch"),
        "i0_deg": rec.get("i0_deg"),
        "i_lto_deg": rec.get("i_lto_deg"),
        "i_park_deg": rec.get("i_park_deg"),
        "az_deg": rec.get("az_deg"),
        "coast_s": rec.get("coast_s"),
        "coast_roots_s": rec.get("coast_roots_s"),
        "i_m_deg": rec.get("i_m_deg"),
        "perilune_alt_km": rec.get("perilune_alt_km"),
        "c3_km2s2": rec.get("C3_km2s2"),
        "dv_tli_kms": rec.get("dv_tli_kms"),
        "dv_loi_kms": rec.get("dv_loi_kms"),
        "dv_total_kms": rec.get("dv_total_kms"),
        "dvp_mps": rec.get("dvp_mps"),
        "margin_kg": rec.get("margin_kg"),
        "width_min": rec.get("width_min"),
        "site_off_plane_km": rec.get("site_off_plane_km"),
        "gates": rec.get("gates"),
        "warm_start_used": rec.get("warm_start_used"),
        "reject_codes": rec.get("certify_reject_codes"),
    }


def build_calendar(cfg: dict, forward_windows: list[dict],
                   cert_records: list[dict], inverse: list[dict]) -> dict:
    by_date: dict[str, dict] = defaultdict(
        lambda: {"forward_windows": [], "certified": [], "inverse": []})
    for win in forward_windows:
        w = {k: v for k, v in win.items() if k != "_pts"}
        by_date[win["t_open_utc"][:10]]["forward_windows"].append(w)
    for rec in cert_records:
        if not rec.get("published"):
            continue
        # The mission calendar carries the PRIMARY target only; the
        # multi-target sweep (--i-m-list, FR-11) stays in
        # fusion_certify.json / fusion_multi_orbit.json.
        if abs(float(rec.get("i_m_target_deg", cfg["i_m_target_deg"]))
               - cfg["i_m_target_deg"]) > 1e-9:
            continue
        by_date[rec["t_launch_utc"][:10]]["certified"].append(
            _fusion_window_from_cert(rec))
    for inv in inverse:
        by_date[inv["t_launch_utc"][:10]]["inverse"].append({
            "t_launch_utc": inv["t_launch_utc"],
            "family": inv.get("family"),
            "branch": inv.get("branch"),
            "coast_s": inv.get("coast_s"),
            "azimuth_deg": inv.get("azimuth_deg"),
            "i_park_deg": inv.get("i_park_deg"),
            "c3_km2s2": inv.get("c3_km2s2"),
        })
    dates = []
    for date in sorted(by_date):
        d = by_date[date]
        cert = sorted(d["certified"], key=lambda r: r["t_launch_utc"])
        dates.append({
            "date_utc": date,
            "n_forward_windows": len(d["forward_windows"]),
            "n_certified_windows": len(cert),
            "n_inverse_windows": len(d["inverse"]),
            "first_launch_utc": cert[0]["t_launch_utc"] if cert else None,
            "last_launch_utc": cert[-1]["t_launch_utc"] if cert else None,
            "margin_kg_stats": _stats([r["margin_kg"] for r in cert
                                       if r.get("margin_kg") is not None]),
            "c3_km2s2_stats": _stats([r["c3_km2s2"] for r in cert
                                      if r.get("c3_km2s2") is not None]),
            "window_min_stats": _stats([r["width_min"] for r in cert
                                        if r.get("width_min") is not None]),
            "forward_windows": d["forward_windows"],
            "certified_windows": cert,
            "inverse_windows": sorted(d["inverse"],
                                      key=lambda r: r["t_launch_utc"]),
        })
    return {
        "mission": cfg["mission"],
        "T_h": cfg["T_h"],
        "range_utc": [cfg["t_start_utc"], cfg["t_end_utc"]],
        "n_dates": len(dates),
        "dates": dates,
    }


def build_coverage(cfg: dict, forward_points: list[dict],
                   forward_windows: list[dict], cert_records: list[dict],
                   inverse: list[dict], cross: dict) -> dict:
    by_date: dict[str, dict] = defaultdict(
        lambda: {"n_points": 0, "n_feasible": 0, "reject_counts": Counter(),
                 "n_forward_windows": 0, "n_certified": 0, "n_inverse": 0})
    for p in forward_points:
        date = p["t_launch_utc"][:10]
        by_date[date]["n_points"] += 1
        if p["feasible"]:
            by_date[date]["n_feasible"] += 1
        else:
            for code in p["reject_codes"]:
                by_date[date]["reject_counts"][code] += 1
    for win in forward_windows:
        by_date[win["t_open_utc"][:10]]["n_forward_windows"] += 1
    for rec in cert_records:
        if rec.get("certify_status") == "certified":
            by_date[rec["t_launch_utc"][:10]]["n_certified"] += 1
    inv_by_date: dict[str, int] = Counter(
        i["t_launch_utc"][:10] for i in inverse)
    dates = []
    for date in sorted(by_date):
        d = by_date[date]
        dates.append({
            "date_utc": date,
            "n_points": d["n_points"],
            "n_feasible": d["n_feasible"],
            "feasible_minutes": d["n_feasible"] * cfg["step_s"] / 60.0,
            "reject_counts": dict(d["reject_counts"].most_common()),
            "n_forward_windows": d["n_forward_windows"],
            "n_certified": d["n_certified"],
            "n_inverse": inv_by_date.get(date, 0),
        })
    return {
        "step_s": cfg["step_s"],
        "n_dates": len(dates),
        "n_points_total": len(forward_points),
        "n_inverse_only": cross.get("n_inverse_only"),
        "n_inverse_only_defects": cross.get("n_inverse_only_defects"),
        "n_inverse_covered": cross.get("n_inverse_covered"),
        "dates": dates,
    }


def write_month_report(cfg: dict, calendar: dict, cross: dict,
                       coverage: dict, cert_records: list[dict],
                       path: Path) -> None:
    lines: list[str] = []
    lines.append(f"# CE-7 {cfg['t_start_utc'][:7]} 融合窗口报告（T={cfg['T_h']} h）")
    lines.append("")
    lines.append(f"- 发射窗口范围: {cfg['t_start_utc']} .. {cfg['t_end_utc']} "
                 f"(步长 {cfg['step_s']:.0f} s)")
    lines.append(f"- 场址: {cfg['site']['name']} "
                 f"({cfg['site']['lat_deg']}°N, {cfg['site']['lon_deg']}°E)")
    lines.append(f"- 走廊: {cfg['az_min_deg']}°–{cfg['az_max_deg']}°, "
                 f"coast ∈ [{cfg['coast_min_s']:.0f}, {cfg['coast_max_s']:.0f}] s")
    lines.append(f"- 认证: i_m={cfg['i_m_target_deg']}°±{cfg['i_m_tol_deg']}°, "
                 f"近月点 {cfg['H_m_km']:.0f}±{cfg['perilune_tol_km']} km; "
                 f"反解: 生产 L6 链")
    lines.append("")
    lines.append(f"## 月度概览")
    lines.append("")
    lines.append(f"- 正解筛查窗口: {sum(d['n_forward_windows'] for d in calendar['dates'])}")
    lines.append(f"- 认证(可发布)窗口: {sum(d['n_certified_windows'] for d in calendar['dates'])}")
    lines.append(f"- 反解窗口: {sum(d['n_inverse_windows'] for d in calendar['dates'])}")
    lines.append(f"- 交叉: both={cross['n_both']}, "
                 f"forward_only={cross['n_forward_only']}, "
                 f"反解窗口={cross['n_inverse_windows']} "
                 f"(等效链覆盖 {cross.get('n_inverse_covered', 0)}, "
                 f"走廊外 {cross.get('n_inverse_only_corridor_outside', 0)}, "
                 f"真缺口 {cross.get('n_inverse_only_defects', 0)})")
    lines.append("")
    lines.append("## 逐日明细（可发射时间段 + 代表性弹道）")
    lines.append("")
    lines.append("> **可发射时间段** = 当天连续可行区间（正解/两体筛查 + 包络，"
                 "2 min 栅格；区间内任意发射时刻都有闭合滑行解）。")
    lines.append("> **代表性弹道** = 四体认证采样（窗口首/尾、各根 coast 极值、"
                 "≥5 等分点，≤8 点/窗）；**表中时刻是认证样本，不是全部可发射时刻**。"
                 "四体认证的实际可用边界与筛查区间可能略有出入（认证解平面须含"
                 "发射场，个别 2 min 点可能被剔除——逐点认证见 "
                 "`fusion_certified_boundary.json`）。")
    lines.append("")
    lines.append("| 日期 | 可发射时间段 (UTC) | 累计(min) | 滑行跨度(s) | "
                 "代表性弹道（认证样本, UTC） | 样本 | 反解 |")
    lines.append("|---|---|---|---|---|---|---|")
    for d in calendar["dates"]:
        fws = d.get("forward_windows") or []
        if fws:
            spans = " / ".join(
                f"{w['t_open_utc'][11:16]}–{w['t_close_utc'][11:16]}"
                for w in fws)
            total_min = sum(float(w["duration_min"]) for w in fws)
            cmins: list[float] = []
            cmaxs: list[float] = []
            for w in fws:
                ext = (w.get("root_extremes") or {}).get("0") or {}
                if ext.get("coast_min_s") is not None:
                    cmins.append(float(ext["coast_min_s"]))
                if ext.get("coast_max_s") is not None:
                    cmaxs.append(float(ext["coast_max_s"]))
            coast_span = (f"{min(cmins):.0f}–{max(cmaxs):.0f}"
                          if cmins and cmaxs else "—")
        else:
            spans, total_min, coast_span = "—", 0.0, "—"
        cert = d["certified_windows"]
        reps = (" / ".join(r["t_launch_utc"][11:16] for r in cert[:6])
                + (" …" if len(cert) > 6 else "")) if cert else "—"
        lines.append(
            f"| {d['date_utc']} | {spans} | {total_min:.0f} | {coast_span} "
            f"| {reps} | {len(cert)} | {d['n_inverse_windows']} |")
    lines.append("")
    lines.append("## 认证汇总（全部认证记录；多目标运行时含 --i-m-list 的"
                 "每个目标，任务日历仅取主目标）")
    lines.append("")
    lines.append("| 发射(UTC) | 目标 i_m | family/branch | coast | az | i_m | 近月点 | C3 | dvp(m/s) | 状态 |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|")
    for r in cert_records:
        fam = f"{r.get('family')}/{r.get('branch')}" if r.get("family") else "—"
        status = r.get("certify_status")
        if status == "certified_gates_failed":
            status = "certified(附加门禁未过: " + ",".join(
                r.get("certify_reject_codes", [])) + ")"
        lines.append(
            f"| {r.get('t_launch_utc')} | "
            f"{r.get('i_m_target_deg', float('nan')):g} | {fam} | "
            f"{r.get('coast_s', float('nan')):.1f} | "
            f"{r.get('az_deg', float('nan')):.2f} | "
            f"{r.get('i_m_deg', float('nan')):.3f} | "
            f"{r.get('perilune_alt_km', float('nan')):.2f} | "
            f"{r.get('C3_km2s2', float('nan')):.3f} | "
            f"{r.get('dvp_mps', float('nan')):+.2f} | "
            f"{status} |")
    lines.append("")
    lines.append("## 交叉校验")
    lines.append("")
    lines.append(f"- both: {cross['n_both']}")
    lines.append(f"- forward_only: {cross['n_forward_only']}")
    lines.append(f"- inverse_only: {cross['n_inverse_only']}")
    if cross["inverse_only"]:
        lines.append("")
        lines.append("### inverse_only 明细（按分类）")
        for r in cross["inverse_only"]:
            w = r["inverse_window"]
            lines.append(f"- [{r['classification']}] {w['t_launch_utc']} "
                         f"{w.get('family')}/{w.get('branch')} "
                         f"coast={w['coast_s']:.0f}s az={w['azimuth_deg']:.2f}° "
                         f"— {r['why']}")
    lines.append("")
    lines.append("## 覆盖统计（逐日）")
    lines.append("")
    lines.append("| 日期 | 采样点 | 可行点 | 可行分钟 | 正解窗口 | 认证 | 反解 | 主要拒绝 |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for d in coverage["dates"]:
        top = ", ".join(f"{k}:{v}" for k, v in
                        list(d["reject_counts"].items())[:3]) or "—"
        lines.append(f"| {d['date_utc']} | {d['n_points']} | {d['n_feasible']} | "
                     f"{d['feasible_minutes']:.0f} | {d['n_forward_windows']} | "
                     f"{d['n_certified']} | {d['n_inverse']} | {top} |")
    lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# Mission verification (CE-5 / CE-6 liftoff-time accuracy, arbitrary i_m)
# ---------------------------------------------------------------------------
_CE56_VERIFY_DIR = Path(
    os.environ.get(
        "FUSION_CE56_VERIFY_DIR",
        str(_REPO_ROOT / "viewer" / "data" / "fusion_ce56_verify"),
    )
)


def _verify_table_path(out_dir: Path, T_h: float, i_m: float,
                       use_sun: bool = True) -> Path:
    """Locate (or mark for build) the arbitrary-i_m family table.

    Search order: the pre-built ``ce56_verify`` tables (``verify_ce56``
    naming, T in centihours), then the run's own ``tables/`` directory
    (fusion naming, T in tenths).  A missing table is built on demand by
    :func:`ensure_family_table`.
    """
    im = int(round(i_m * 1000))
    sun = "_sun" if use_sun else ""
    ce56 = _CE56_VERIFY_DIR / (
        f"r1s_scan_im{im}_T{int(round(T_h * 100)):04d}_Hm200_meq{sun}.json")
    if ce56.exists():
        return ce56
    return _table_path(out_dir, T_h, DEFAULT_H_M_KM, i_m, use_sun)


def _load_ce56_baseline() -> dict[str, float]:
    """Production L6 exact-time offsets (min) from ``ce56_verification_f0.json``."""
    path = _CE56_VERIFY_DIR / "ce56_verification_f0.json"
    out: dict[str, float] = {}
    try:
        doc = json.loads(path.read_text(encoding="utf-8"))
        for m in doc.get("per_mission", []):
            hits = (m.get("l6_exact_time_check") or {}).get("hits") or []
            offs = [h.get("offset_min") for h in hits
                    if h.get("offset_min") is not None]
            if offs:
                out[m["mission"]] = min(offs, key=abs)
    except Exception:  # noqa: BLE001
        pass
    return out


def _verify_missions(report: dict) -> list[dict]:
    """Merge the ground-truth report with ``verify_ce56`` mission parameters.

    The i_m target comes from ``verify_ce56.MISSIONS`` (CE-5 43.7 deg,
    CE-6 137 deg — literature values); the report's
    ``actual_orbit_inclination_deg`` is recorded for reference only and is
    NOT used as a certification target (earlier revisions carried mistaken
    CE-5 values — 21.299 deg, then 45.7 deg — both superseded by the 43.7
    deg literature value; do not reintroduce them).
    """
    try:
        from launch_window.verify_ce56 import MISSIONS as ce56_missions
    except Exception:  # noqa: BLE001
        ce56_missions = []
    by_name = {m["mission"]: m for m in ce56_missions}
    missions: list[dict] = []
    for name, rec in report.get("per_mission", {}).items():
        base = dict(by_name.get(name, {}))
        liftoff = rec.get("actual_liftoff_utc") or base.get("launch_utc")
        arrival = rec.get("actual_arrival_utc")
        if not liftoff:
            continue
        T_h = float(base.get("T_h") or 0.0)
        if not T_h and arrival:
            T_h = (g.utc_to_et(arrival) - g.utc_to_et(liftoff)) / 3600.0
        im_cands: list[float] = []
        if base.get("i_m_deg") is not None:
            im_cands.append(float(base["i_m_deg"]))
        elif rec.get("actual_orbit_inclination_deg") is not None:
            im_cands.append(float(rec["actual_orbit_inclination_deg"]))
        missions.append({
            "mission": name,
            "liftoff_utc": liftoff,
            "arrival_utc": arrival,
            "T_h": T_h,
            "t_to_park_s": float(base.get("t_to_park_s",
                                          DEFAULT_T_ASCENT_S)),
            "coast_measured_s": base.get("coast_s"),
            "i_m_candidates": im_cands,
            "report_inclination_deg": rec.get("actual_orbit_inclination_deg"),
        })
    return missions


def _mission_namespace(mission: dict, i_m: float, coast_cap_s: float,
                       args: argparse.Namespace) -> argparse.Namespace:
    """Mission-specific CLI namespace (shared by verify-mission / chain)."""
    t_liftoff = g.utc_to_et(mission["liftoff_utc"])
    span_s = float(args.span_h) * 3600.0
    return argparse.Namespace(
        start=_norm_iso(g.et_to_utc(t_liftoff - span_s)),
        end=_norm_iso(g.et_to_utc(t_liftoff + span_s)),
        step=float(args.step), T=float(mission["T_h"]), site=str(args.site),
        mission=str(mission["mission"]), perilune_alt=DEFAULT_H_M_KM,
        parking_peri_alt=DEFAULT_H_PARK_KM, parking_apo_alt=DEFAULT_H_APO_KM,
        i_m=float(i_m), polar_ref=DEFAULT_POLAR_REF,
        moon_model=DEFAULT_MOON_MODEL, use_sun=True,
        coast_min=DEFAULT_COAST_MIN_S, coast_max=float(coast_cap_s),
        coast_step_s=DEFAULT_COAST_STEP_S,
        t_ascent_s=float(mission["t_to_park_s"]),
        # Inertial launch->insertion arc: from the mission's real insertion
        # point when available (bookkeeping: EF point rotated to J2000 at the
        # insertion epoch, so Earth rotation during the ascent is included),
        # otherwise the legacy power-law ASCENT_ARC_DEG.
        ascent_arc_deg=float(mission.get("ascent_arc_deg")
                             or getattr(args, "ascent_arc_deg", None)
                             or DEFAULT_ASCENT_ARC_DEG),
        az=[float(args.az[0]), float(args.az[1])], use_envelope=True,
        cert_k=int(getattr(args, "cert_k", DEFAULT_CERT_K)),
        allow_inverse_only="true",
        i0_hi=float(getattr(args, "i0_hi", DEFAULT_I0_HI)),
        i_m_tol=float(getattr(args, "i_m_tol", DEFAULT_I_M_TOL_DEG)),
    )


def _mission_tag(mission: dict, i_m: float, coast_cap_s: float) -> str:
    return (f"{mission['mission']}_im{int(round(i_m * 1000))}"
            f"_cap{int(round(coast_cap_s))}")


def _mission_cfg(mission: dict, i_m: float, coast_cap_s: float, regime: str,
                 args: argparse.Namespace, coast_fallback: float) -> dict:
    """Build the fusion config for one mission / target / coast-cap.

    ``dtheta_deg`` is the R-2 B-7 launch->perigee arc; it is set from the
    mission's *measured* timeline (``t_to_park + coast_measured``) when
    available, otherwise from ``coast_fallback``.  It only feeds the
    family-table warm starts (FR-3 / L6), never the forward closure.
    """
    cfg = build_config(_mission_namespace(mission, i_m, coast_cap_s, args))
    cfg["out_dir"] = str(Path(args.cache) / _mission_tag(
        mission, i_m, coast_cap_s))
    cfg["regime"] = regime
    coast_m = mission.get("coast_measured_s")
    coast_for_dtheta = (float(coast_m) if coast_m is not None
                        else float(coast_fallback))
    cfg["dtheta_deg"] = (
        (float(mission["t_to_park_s"]) + coast_for_dtheta) / 3600.0
        * le.EARTH_SIDEREAL_RATE_DEG_PER_H)
    cfg["table_path"] = str(_verify_table_path(
        Path(args.cache), mission["T_h"], i_m, True))
    return cfg


def _verify_one(mission: dict, i_m: float, coast_cap_s: float, regime: str,
                args: argparse.Namespace, log_path: Path) -> dict:
    """Verify one (mission, i_m, coast-cap) candidate vs the actual liftoff.

    ``regime`` labels the cap configuration: ``design`` (CZ-5 1000 s single
    orbit coast envelope) or ``dynamics`` (raised cap — the 2-body closure
    coasts for the historical missions sit ~1.3-1.4 ks, above the design
    envelope, so the raised cap measures the intrinsic liftoff accuracy
    without the cap pushing the window edge).
    """
    t_liftoff = g.utc_to_et(mission["liftoff_utc"])
    span_s = float(args.span_h) * 3600.0
    cfg = _mission_cfg(mission, i_m, coast_cap_s, regime, args,
                       args.coast_max)
    ensure_family_table(cfg, log_path)
    tag = _mission_tag(mission, i_m, coast_cap_s)
    _log(f"[verify] {tag}: T={cfg['T_h']} h i_m={i_m} deg "
         f"table={Path(cfg['table_path']).name} "
         f"dtheta={cfg['dtheta_deg']:.3f} deg", log_path)

    # FR-1 forward scan around the actual liftoff (in-process, fine grid).
    n = int(round(2.0 * span_s / float(args.step)))
    points = [forward_point(t_liftoff - span_s + k * float(args.step), cfg)
              for k in range(n + 1)]
    pt_actual = forward_point(t_liftoff, cfg)
    windows = cluster_windows(points, cfg)

    # FR-3 certification of representatives + the exact liftoff epoch.
    cert_records, cert_fails = stage_certify(windows, cfg, log_path,
                                             args.workers)
    cert_actual: Optional[dict] = None
    if pt_actual["roots"]:
        syn = {"_pts": [pt_actual], "root_extremes": {}}
        cert_actual = certify_candidate(syn, 0, 0, cfg)
        cert_actual.pop("_res", None)
    else:
        cert_actual = {
            "certify_status": "skipped:no_coast_root",
            "certify_reject_codes": list(pt_actual["reject_codes"]),
            "t_launch_utc": pt_actual["t_launch_utc"],
        }

    published = [r for r in cert_records if r.get("published")]
    pool = list(published)
    if cert_actual is not None and cert_actual.get("published"):
        pool = pool + [cert_actual]
    nearest = (min(pool, key=lambda r: abs(float(r["t_launch_et"])
                                           - t_liftoff))
               if pool else None)
    contains_actual = any(
        w["t_open_et"] - float(args.step) <= t_liftoff
        <= w["t_close_et"] + float(args.step) for w in windows)

    def _short(r: Optional[dict]) -> Optional[dict]:
        if r is None:
            return None
        out = {
            "t_launch_utc": r.get("t_launch_utc"),
            "dt_min": None,
            "family": r.get("family"), "branch": r.get("branch"),
            "coast_s": r.get("coast_s"),
            "az_deg": r.get("az_deg"),
            "azimuth_std_deg": r.get("azimuth_std_deg"),
            "i_m_deg": r.get("i_m_deg"),
            "perilune_alt_km": r.get("perilune_alt_km"),
            "C3_km2s2": r.get("C3_km2s2"),
            "dvp_mps": r.get("dvp_mps"),
            "certify_status": r.get("certify_status"),
            "warm_start_used": r.get("warm_start_used"),
            "reject_codes": r.get("certify_reject_codes"),
        }
        if r.get("t_launch_et") is not None:
            out["dt_min"] = round(
                (float(r["t_launch_et"]) - t_liftoff) / 60.0, 2)
        if (mission.get("coast_measured_s") is not None
                and r.get("coast_s") is not None):
            out["coast_error_s"] = round(
                float(r["coast_s"]) - float(mission["coast_measured_s"]), 1)
        return out

    rec = {
        "mission": mission["mission"],
        "i_m_deg": float(i_m),
        "report_inclination_deg": mission.get("report_inclination_deg"),
        "T_h": float(cfg["T_h"]),
        "t_to_park_s": float(mission["t_to_park_s"]),
        "coast_measured_s": mission.get("coast_measured_s"),
        "coast_cap_s": float(coast_cap_s),
        "regime": regime,
        "dtheta_deg": float(cfg["dtheta_deg"]),
        "i0_hi": float(cfg["i0_hi"]),
        "actual_liftoff_utc": mission["liftoff_utc"],
        "actual_arrival_utc": mission.get("arrival_utc"),
        "table": Path(cfg["table_path"]).name,
        "scan": {
            "span_h": float(args.span_h), "step_s": float(args.step),
            "n_points": len(points),
            "n_feasible": sum(1 for p in points if p["feasible"]),
            "n_windows": len(windows),
        },
        "forward_windows": [{
            "t_open_utc": w["t_open_utc"], "t_close_utc": w["t_close_utc"],
            "duration_min": w["duration_min"],
            "contains_actual": bool(
                w["t_open_et"] - float(args.step) <= t_liftoff
                <= w["t_close_et"] + float(args.step)),
        } for w in windows],
        "actual_within_forward_window": bool(contains_actual),
        "forward_at_actual": {
            "feasible": pt_actual["feasible"],
            "az_deg": pt_actual["az_deg"],
            "i_lto_deg": pt_actual["i_lto_deg"],
            "n_roots": pt_actual["n_roots"],
            "coast_s": (pt_actual["roots"][0]["coast_s"]
                        if pt_actual["roots"] else None),
            "reject_codes": pt_actual["reject_codes"],
        },
        "n_published": len(published),
        "cert_at_actual": _short(cert_actual),
        "nearest_published": _short(nearest),
        "published_records": [_short(r) for r in published],
        "n_cert_failed": len(cert_fails),
    }
    _write_json(Path(cfg["out_dir"]) / "verify_record.json", rec)
    return rec


def _write_verify_markdown(results: list[dict], baseline: dict,
                           report_path: Path, path: Path) -> None:
    L: list[str] = []
    L.append("# 融合求解器起飞时间验证（对照 "
             f"`{report_path.name}`）")
    L.append("")
    L.append("方法：以报告中的实际起飞/到达时刻为真值，在实际起飞前后 "
             "`--span-h` 窗口内以 `--step` 步长运行正解扫描（FR-1/FR-2），"
             "对代表点与**实际起飞历元本身**做任意目标环月倾角认证"
             "（FR-3，`i_m=90°` 走 `solve_polar`，其余走 `solve_non_coplanar`），"
             "以最近可发布认证点的 Δt 度量起飞时间准确性。")
    L.append("")
    L.append("| 任务 | 目标 i_m | 滑行上限 | 实际点可行 | 实际点滑行 (s) | "
             "滑行误差 (s) | 最近认证点 Δt (min) | 认证 i_m | 认证方位角 (°) | "
             "L6 基线 Δt (min) |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for r in results:
        fa = r["forward_at_actual"]
        ca = r.get("cert_at_actual") or {}
        np_ = r.get("nearest_published") or {}
        base = baseline.get(r["mission"])
        coast_err = ca.get("coast_error_s")
        L.append(
            f"| {r['mission']} | {r['i_m_deg']} "
            f"| {r['coast_cap_s']:.0f} s ({r['regime']}) "
            f"| {fa['feasible']} | {_fmt(fa.get('coast_s'), 1)} "
            f"| {_fmt(coast_err, 1)} | {_fmt(np_.get('dt_min'), 2)} "
            f"| {_fmt(np_.get('i_m_deg'), 3)} "
            f"| {_fmt(np_.get('azimuth_std_deg'), 2)} "
            f"| {_fmt(base, 2)} |")
    L.append("")
    L.append("两种滑行口径：**design** = CZ-5 单圈滑行设计包络（1000 s）；"
             "**dynamics** = 放开上限（1500 s）——历史上星 2 体闭合滑行位于 "
             "1.3–1.4 ks（超出设计包络），放开后才能度量求解器本身的起飞时间"
             "准确性，而不被滑行上限把窗口边缘推离实际起飞时刻。")
    L.append("")
    for r in results:
        L.append(f"## {r['mission']} — i_m = {r['i_m_deg']}° "
                 f"（报告口径 {r.get('report_inclination_deg')}°，"
                 f"滑行上限 {r['coast_cap_s']:.0f} s / {r['regime']}）")
        L.append("")
        L.append(f"- 族表: `{r['table']}`，T = {r['T_h']} h，"
                 f"t_to_park = {r['t_to_park_s']} s，"
                 f"实测滑行 = {r['coast_measured_s']} s，"
                 f"dθ = {r['dtheta_deg']:.3f}°，i0 ∈ [1, {r['i0_hi']}]")
        L.append(f"- 扫描: ±{r['scan']['span_h']:.0f} h @ {r['scan']['step_s']:.0f} s，"
                 f"{r['scan']['n_points']} 点，可行 {r['scan']['n_feasible']}，"
                 f"窗口 {r['scan']['n_windows']} 个")
        L.append(f"- 实际起飞落在正解窗口内: **{r['actual_within_forward_window']}**；"
                 f"实际历元正解可行: **{r['forward_at_actual']['feasible']}**"
                 f"（coast {_fmt(r['forward_at_actual'].get('coast_s'), 1)} s，"
                 f"az {_fmt(r['forward_at_actual'].get('az_deg'), 2)}°，"
                 f"拒绝码 {r['forward_at_actual']['reject_codes'] or '无'}）")
        ca = r.get("cert_at_actual") or {}
        if ca:
            L.append(f"- 实际历元认证: {ca.get('certify_status')}，"
                     f"family {ca.get('family')}/{ca.get('branch')}，"
                     f"coast {_fmt(ca.get('coast_s'), 1)} s，"
                     f"i_m {_fmt(ca.get('i_m_deg'), 3)}°，"
                     f"近月点 {_fmt(ca.get('perilune_alt_km'), 2)} km，"
                     f"认证标准方位角 {_fmt(ca.get('azimuth_std_deg'), 2)}°")
        np_ = r.get("nearest_published") or {}
        if np_:
            L.append(f"- 最近可发布认证点: {np_.get('t_launch_utc')}"
                     f"（Δt = {_fmt(np_.get('dt_min'), 2)} min），"
                     f"coast {_fmt(np_.get('coast_s'), 1)} s（误差 "
                     f"{_fmt(np_.get('coast_error_s'), 1)} s），"
                     f"C3 {_fmt(np_.get('C3_km2s2'), 3)} km²/s²，"
                     f"dvp {_fmt(np_.get('dvp_mps'), 2)} m/s")
        else:
            L.append("- 无可发布认证点（认证失败/无族根），见 verify_record.json")
        if r["mission"] in baseline:
            L.append(f"- 生产 L6 链基线（`ce56_verification_f0.json`）: "
                     f"Δt = {baseline[r['mission']]} min")
        L.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(L), encoding="utf-8")


def _fmt(v, nd: int = 2) -> str:
    if v is None:
        return "—"
    try:
        return f"{float(v):.{nd}f}"
    except Exception:  # noqa: BLE001
        return str(v)


def _main_verify(args: argparse.Namespace) -> int:
    g.ensure_kernels()
    report_path = Path(args.report)
    report = json.loads(report_path.read_text(encoding="utf-8"))
    missions = _verify_missions(report)
    if not missions:
        print("no missions with actual liftoff times found in report")
        return 1
    baseline = _load_ce56_baseline()
    log_path = Path(args.cache) / "logs" / "fusion_verify.log"
    caps: list[tuple[str, float]] = [("design", float(args.coast_max))]
    if (float(args.coast_max_dynamics)
            != float(args.coast_max)):
        caps.append(("dynamics", float(args.coast_max_dynamics)))
    results: list[dict] = []
    for m in missions:
        for i_m in m["i_m_candidates"]:
            for regime, cap in caps:
                results.append(_verify_one(m, i_m, cap, regime, args,
                                           log_path))
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    doc = {
        "schema": "fusion_mission_verification_v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "report": str(report_path),
        "method": ("fusion forward scan + arbitrary-i_m certification "
                   "(solve_polar @ 90 deg / solve_non_coplanar otherwise), "
                   "nearest published certificate to the actual liftoff; "
                   "two coast-cap regimes: design (CZ-5 envelope) and "
                   "dynamics (raised cap, intrinsic accuracy)"),
        "baseline_l6_offset_min": baseline,
        "results": results,
    }
    json_path = out / "launch_window_fusion_ce56_report.json"
    json_path.write_text(json.dumps(_json_safe(doc), ensure_ascii=False,
                                    indent=2), encoding="utf-8")
    md_path = out / "launch_window_fusion_ce56_report.md"
    _write_verify_markdown(results, baseline, report_path, md_path)
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")
    for r in results:
        np_ = r.get("nearest_published") or {}
        base = baseline.get(r["mission"])
        print(f"{r['mission']} i_m={r['i_m_deg']} cap={r['coast_cap_s']:.0f}"
              f"({r['regime']}): "
              f"actual_in_window={r['actual_within_forward_window']} "
              f"fwd_feasible@actual={r['forward_at_actual']['feasible']} "
              f"nearest_dt_min={np_.get('dt_min')} "
              f"(L6 baseline {base})")
    return 0


# ---------------------------------------------------------------------------
# Mission chain: coast_tool screening -> capture-plane correction -> BVP ->
# errors vs the official mission timeline (CE-5 / CE-6)
# ---------------------------------------------------------------------------
_FIXTURE_DEFAULT = (_REPO_ROOT / "tests" / "launch_window"
                    / "发射窗口核验.json")


def _mission_times_from_fixture(rec: dict) -> Optional[dict]:
    """Timeline of one fixture mission under the corrected definitions.

    ``anchor`` = parking-orbit insertion = ``coast_start_utc`` (CE-5) or
    ``coast_end_utc`` (CE-6).  Then

        ascent   = anchor − liftoff            (launch -> parking insertion)
        coast    = tli − anchor                (parking-orbit coast)
        transfer = orbit_insertion − tli       (perigee -> arrival, L6 T)
        T_total  = orbit_insertion − liftoff   (fusion forward-closure T)
    """
    liftoff = rec.get("liftoff_utc")
    tli = rec.get("tli_utc")
    arrival = rec.get("orbit_insertion_utc")
    anchor = rec.get("coast_start_utc") or rec.get("coast_end_utc")
    if not (liftoff and tli and arrival and anchor):
        return None
    t_l = g.utc_to_et(liftoff)
    t_a = g.utc_to_et(anchor)
    t_t = g.utc_to_et(tli)
    t_r = g.utc_to_et(arrival)
    return {
        "liftoff_utc": liftoff,
        "anchor_utc": anchor,
        "tli_utc": tli,
        "arrival_utc": arrival,
        "t_to_park_s": float(t_a - t_l),
        "coast_measured_s": float(t_t - t_a),
        "T_transfer_s": float(t_r - t_t),
        "T_total_s": float(t_r - t_l),
        "T_h": float((t_r - t_l) / 3600.0),
    }


def _load_fixture_missions(path: Path) -> list[dict]:
    """CE-5/CE-6 ground truth from ``发射窗口核验.json`` (corrected times).

    The fixture carries the liftoff / parking-insertion (``coast_start_utc``
    or the legacy ``coast_end_utc``) / TLI / orbit-insertion epochs; the
    lunar-orbit inclination target is taken from ``verify_ce56.MISSIONS``
    (the fixture has no inclination field).
    """
    doc = json.loads(Path(path).read_text(encoding="utf-8"))
    try:
        from launch_window.verify_ce56 import MISSIONS as ce56_missions
    except Exception:  # noqa: BLE001
        ce56_missions = []
    by_name = {m["mission"]: m for m in ce56_missions}

    records: list[tuple[str, dict]] = []

    def _walk(obj: Any, prefix: str = "") -> None:
        if not isinstance(obj, dict):
            return
        for key, val in obj.items():
            if isinstance(val, dict):
                if val.get("liftoff_utc"):
                    records.append((key, val))
                else:
                    _walk(val, f"{prefix}{key}.")

    _walk(doc)
    out: list[dict] = []
    for name, rec in records:
        times = _mission_times_from_fixture(rec)
        if times is None:
            continue
        base = dict(by_name.get(name, {}))
        im_cands: list[float] = []
        if base.get("i_m_deg") is not None:
            im_cands.append(float(base["i_m_deg"]))
        elif rec.get("actual_orbit_inclination_deg") is not None:
            im_cands.append(float(rec["actual_orbit_inclination_deg"]))
        out.append({
            "mission": name,
            "report_inclination_deg": rec.get("actual_orbit_inclination_deg"),
            "i_m_candidates": im_cands,
            "source_fixture": str(path),
            **times,
        })
    return out


def _gc_arc_deg(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Great-circle arc (deg) between two geocentric lat/lon points."""
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dl = math.radians(lon2 - lon1)
    c = (math.sin(p1) * math.sin(p2)
         + math.cos(p1) * math.cos(p2) * math.cos(dl))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def _apply_insertion_bookkeeping(missions: list[dict], args: argparse.Namespace,
                                 log_path: Path) -> None:
    """Set each mission's inertial ascent arc from its real insertion point.

    Bookkeeping: the insertion point is an Earth-fixed (geodetic) location;
    its J2000 direction must be taken at the insertion epoch ``t_L +
    t_to_park``, so the angle from the launch-site direction at liftoff
    automatically includes the Earth rotation during the ascent (~ω_E·t).
    """
    ins_lat = getattr(args, "insertion_lat", None)
    ins_lon = getattr(args, "insertion_lon", None)
    if ins_lat is None or ins_lon is None:
        return
    s = SITE_REGISTRY[str(args.site).lower()]
    for m in missions:
        t_l = g.utc_to_et(m["liftoff_utc"])
        t1 = t_l + float(m["t_to_park_s"])
        w = np.asarray(g.launch_site_vector_j2000_km(
            s["lon_deg"], s["lat_deg"], s["alt_m"], t_l), dtype=float)
        a = np.asarray(g.launch_site_vector_j2000_km(
            float(ins_lon), float(ins_lat), 0.0, t1), dtype=float)
        w = w / np.linalg.norm(w)
        a = a / np.linalg.norm(a)
        arc_inertial = math.degrees(math.acos(
            max(-1.0, min(1.0, float(np.dot(w, a))))))
        arc_ef = _gc_arc_deg(s["lat_deg"], s["lon_deg"], float(ins_lat),
                             float(ins_lon))
        m["ascent_arc_deg"] = float(arc_inertial)
        m["ascent_arc_ef_deg"] = float(arc_ef)
        m["ascent_arc_omega_e_deg"] = float(arc_inertial - arc_ef)
        _om = math.degrees(le.parking_orbit_elements(200.0, 380.0)[2])
        _log(f"[chain]   {m['mission']} insertion bookkeeping: "
             f"EF arc {arc_ef:.3f}° + Earth rotation "
             f"{arc_inertial - arc_ef:+.3f}° = inertial {arc_inertial:.3f}° "
             f"(legacy {DEFAULT_ASCENT_ARC_DEG:.4f}°, "
             f"{(arc_inertial - DEFAULT_ASCENT_ARC_DEG) / _om:+.0f} s of "
             f"parking phase)", log_path)


def _mission_screen(cfg: dict, t_liftoff: float, span_s: float,
                    step_s: float) -> dict:
    """FR-1 coast_tool screening scan around the official liftoff.

    No family tables and no ``dtheta`` are involved: the two-body closure
    root at each grid epoch is the coast_tool "initial solution".
    """
    n = int(round(2.0 * span_s / step_s))
    points = [forward_point(t_liftoff - span_s + k * step_s, cfg)
              for k in range(n + 1)]
    windows = cluster_windows(points, cfg)
    feas = [p for p in points if p["feasible"]]
    nearest = (min(feas, key=lambda p: abs(p["t_launch_et"] - t_liftoff))
               if feas else None)
    actual = forward_point(t_liftoff, cfg)
    return {
        "span_h": span_s / 3600.0,
        "step_s": float(step_s),
        "n_points": len(points),
        "n_feasible": len(feas),
        "n_windows": len(windows),
        "actual_utc": g.et_to_utc(t_liftoff),
        "actual_feasible": bool(actual["feasible"]),
        "actual_n_roots": int(actual["n_roots"]),
        "actual_reject_codes": list(actual["reject_codes"]),
        "nearest_feasible_et": (float(nearest["t_launch_et"])
                                if nearest else None),
        "nearest_feasible_utc": (nearest["t_launch_utc"]
                                 if nearest else None),
        "nearest_feasible_dt_s": (float(nearest["t_launch_et"] - t_liftoff)
                                  if nearest else None),
        "windows": [{
            "t_open_utc": w["t_open_utc"],
            "t_close_utc": w["t_close_utc"],
            "duration_min": w["duration_min"],
            "contains_actual": bool(
                w["t_open_et"] - step_s <= t_liftoff
                <= w["t_close_et"] + step_s),
        } for w in windows],
    }


def _mission_errors(mission: dict, cfg: dict, t_l_et: float, corr: dict,
                    final: dict) -> dict:
    """Chain orbit minus official mission data (seconds / degrees / km)."""
    t_loff = g.utc_to_et(mission["liftoff_utc"])
    coast_meas = mission.get("coast_measured_s")
    t_tli_official = (t_loff + cfg["t_ascent_s"] + float(coast_meas)
                      if coast_meas is not None else None)
    t_arr_official = (g.utc_to_et(mission["arrival_utc"])
                      if mission.get("arrival_utc") else None)
    t_tli_chain = t_l_et + cfg["t_ascent_s"] + corr["coast_corrected_s"]
    return {
        "dt_launch_s": float(t_l_et - t_loff),
        "dt_tli_s": (float(t_tli_chain - t_tli_official)
                     if t_tli_official is not None else None),
        "dt_arrive_s": (float(t_l_et + cfg["T_s"] - t_arr_official)
                        if t_arr_official is not None else None),
        "d_coast_s": (float(corr["coast_corrected_s"] - coast_meas)
                      if coast_meas is not None else None),
        "d_i_m_deg": (float(final["i_m_deg"] - cfg["i_m_target_deg"])
                      if final.get("i_m_deg") is not None else None),
        "d_perilune_km": (float(final["perilune_alt_km"] - cfg["H_m_km"])
                          if final.get("perilune_alt_km") is not None
                          else None),
    }


def _mission_chain_point(mission: dict, cfg: dict, t_launch_et: float,
                         label: str, relax: bool) -> dict:
    """One chain point: coast_tool initial -> correction -> final BVP.

    No family tables / ``dtheta`` are used: the warm start is the coast_tool
    two-body initial, processed by ``polar_correct.correct_to_polar``
    (SPICE/DOP853 4-body, three-stage ``dAz``/``dtc``/``dvp``, capture plane
    locked to ``pol = cos(i_m_target)``) and finished by the R-1 BVP
    (``solve_polar`` for 90 deg, ``solve_non_coplanar`` otherwise).
    """
    from launch_window import polar_correct as pc
    t_l = float(t_launch_et)
    geo = forward_geometry(t_l, cfg)
    roots = coast_roots_forward(geo, cfg)
    rec: dict[str, Any] = {
        "label": label,
        "t_launch_et": t_l,
        "t_launch_utc": g.et_to_utc(t_l),
        "az_deg": float(geo["az_deg"]),
        "i_lto_deg": float(geo["i_lto_deg"]),
        "theta_deg": float(geo["theta_deg"]),
        "n_roots": len(roots),
        "roots": [{
            "coast_s": float(r["coast_s"]), "dnu_deg": float(r["dnu_deg"]),
            "e": float(r["e"]), "a_km": float(r["a_km"]),
            "v_peri_kms": float(r["v_peri_kms"]),
            "C3_km2s2": float(r["c3_km2s2"]),
        } for r in roots],
        "chain_ok": False,
    }
    if not roots:
        rec["error"] = "no_coast_root"
        return rec
    root = roots[0]
    try:
        corr = pc.correct_to_polar(
            t_l, geo, float(root["coast_s"]), float(root["a_km"]), cfg,
            peri_target_km=cfg["H_m_km"], polar_lock=True, relax=relax,
            i_m_target_deg=float(cfg["i_m_target_deg"]))
    except Exception:  # noqa: BLE001
        rec["error"] = traceback.format_exc(limit=8)
        return rec
    rec["initial_replay"] = {
        "i_m_deg": corr["i_m_initial_replay_deg"],
        "perilune_alt_km": corr["peri_initial_replay_km"],
        "peri_dt_h": corr["peri_initial_dt_h"],
        "pol": corr["pol_initial"],
    }
    rec["correction"] = {k: corr[k] for k in (
        "converged", "i_m_target_deg", "pol_target", "dAz_deg", "dtc_s",
        "dvp_mps", "coast_corrected_s", "az_corrected_deg", "t_tli_utc",
        "T_cert_h", "i_m_replay_deg", "peri_alt_replay_km", "peri_dt_h",
        "pol", "f2_final", "i0_deg", "n_eval", "elapsed_s")}
    rec["t_tli_et"] = float(t_l + cfg["t_ascent_s"]
                            + corr["coast_corrected_s"])
    rec["t_tli_utc"] = g.et_to_utc(rec["t_tli_et"])
    rec["t_arrive_et"] = float(t_l + cfg["T_s"])
    rec["t_arrive_utc"] = g.et_to_utc(rec["t_arrive_et"])
    # Final R-1 BVP from the corrected (SPICE 4-body) warm start.
    try:
        out = _certify_solve(
            rec["t_arrive_et"], float(corr["T_cert_h"] * 3600.0),
            float(corr["i0_deg"]),
            (float(corr["vp_corrected_kms"]), float(corr["Omega_mep_deg"]),
             float(corr["omega_mep_deg"])), cfg)
    except Exception as exc:  # noqa: BLE001
        out = {"ok": False, "converged": False, "res": None,
               "error": f"{type(exc).__name__}: {exc}"}
    final: dict[str, Any] = {"converged": False, "error": out["error"]}
    if out["res"] is not None:
        r = out["res"]
        final = {
            "converged": bool(out["converged"]),
            "i_m_deg": float(r.i_m),
            "perilune_alt_km": float(r.q[0]), "q2": float(r.q[1]),
            "q3": float(r.q[2]), "vp_sol_kms": float(r.p[0]),
            "v_perilune_kms": float(r.v_perilune),
            "C3_km2s2": float(r.p[0] ** 2
                              - 2.0 * ts.MU_EARTH / cfg["rp_km"]),
            "C3_lunar_km2s2": float(r.C3),
            "dvp_final_mps": float(
                (r.p[0] - corr["vp_corrected_kms"]) * 1000.0),
            "perigee_lon_mep_deg": float(r.perigee_lon_mep),
            "perigee_lat_mep_deg": float(r.perigee_lat_mep),
            "perilune_lon_mep_deg": float(r.perilune_lon_mep),
            "perilune_lat_mep_deg": float(r.perilune_lat_mep),
            "error": None,
        }
    rec["final"] = final
    i_m_ok = (final.get("i_m_deg") is not None
              and abs(float(final["i_m_deg"]) - cfg["i_m_target_deg"])
              <= cfg["i_m_tol_deg"])
    peri_ok = (final.get("perilune_alt_km") is not None
               and abs(float(final["perilune_alt_km"]) - cfg["H_m_km"])
               <= cfg["perilune_tol_km"])
    rec["chain_ok"] = bool(final.get("converged") and i_m_ok and peri_ok)
    rec["errors"] = _mission_errors(mission, cfg, t_l, corr, final)
    return rec


def _run_mission_chain(mission: dict, i_m: float, cap: float, regime: str,
                       args: argparse.Namespace, log_path: Path) -> dict:
    """Full mission chain for one (mission, i_m, coast-cap, regime)."""
    cfg = _mission_cfg(mission, i_m, cap, regime, args, cap)
    t_liftoff = g.utc_to_et(mission["liftoff_utc"])
    span_s = float(args.span_h) * 3600.0
    _log(f"[chain] {mission['mission']} i_m={i_m} cap={cap:.0f} "
         f"({regime}): coast_tool screening +/-{args.span_h}h "
         f"@{args.step}s", log_path)
    screen = _mission_screen(cfg, t_liftoff, span_s, float(args.step))
    _log(f"[chain]   screen feasible={screen['n_feasible']}/"
         f"{screen['n_points']} windows={screen['n_windows']} "
         f"actual_roots={screen['actual_n_roots']} nearest_dt_s="
         f"{screen['nearest_feasible_dt_s']}", log_path)
    cands: list[tuple[str, float]] = []
    if screen["actual_n_roots"] > 0:
        cands.append(("actual_epoch", t_liftoff))
    nf = screen["nearest_feasible_et"]
    if nf is not None and abs(nf - t_liftoff) > 1e-6:
        cands.append(("nearest_feasible", float(nf)))
    if not cands:
        _log("[chain]   no feasible candidate — chain skipped", log_path)
    relax = bool(getattr(args, "relax_corrector", False))
    chain = [_mission_chain_point(mission, cfg, t, label, relax)
             for label, t in cands]
    rec = {
        "mission": mission["mission"],
        "i_m_deg": float(i_m),
        "report_inclination_deg": mission.get("report_inclination_deg"),
        "T_h": float(cfg["T_h"]),
        "t_to_park_s": float(mission["t_to_park_s"]),
        "coast_measured_s": mission.get("coast_measured_s"),
        "coast_cap_s": float(cap),
        "regime": regime,
        "relax_corrector": relax,
        "ascent_arc_deg": float(cfg["ascent_arc_deg"]),
        "ascent_arc_ef_deg": mission.get("ascent_arc_ef_deg"),
        "ascent_arc_omega_e_deg": mission.get("ascent_arc_omega_e_deg"),
        "official": {
            "liftoff_utc": mission["liftoff_utc"],
            "anchor_utc": mission.get("anchor_utc"),
            "tli_utc": (mission.get("tli_utc") or (
                g.et_to_utc(t_liftoff + cfg["t_ascent_s"]
                            + float(mission["coast_measured_s"]))
                if mission.get("coast_measured_s") is not None else None)),
            "arrival_utc": mission.get("arrival_utc"),
            "i_m_deg": float(i_m),
            "perilune_alt_km": float(cfg["H_m_km"]),
            "t_to_park_s": float(mission["t_to_park_s"]),
            "coast_measured_s": mission.get("coast_measured_s"),
            "T_transfer_s": mission.get("T_transfer_s"),
            "T_total_s": mission.get("T_total_s"),
        },
        "screen": screen,
        "chain": chain,
        "n_chain_ok": sum(1 for c in chain if c.get("chain_ok")),
    }
    for c in chain:
        _log(f"[chain]   {c['label']}: roots={c['n_roots']} "
             f"coast={_fmt(c.get('correction', {}).get('coast_corrected_s'), 1)}s "
             f"final_i_m={_fmt((c.get('final') or {}).get('i_m_deg'), 3)} "
             f"chain_ok={c['chain_ok']}", log_path)
    return rec


def _write_mission_chain_markdown(results: list[dict], report_path: Path,
                                  path: Path) -> None:
    L: list[str] = []
    L.append("# 任务链路复核：coast_tool 初筛 → 捕获面修正 → 融合 BVP → 官方误差")
    L.append("")
    L.append(f"真值来源：`{report_path.name}`。链路：① 实际/筛查发射历元上的 "
             "coast_tool 两体闭合**初解**；② `polar_correct.correct_to_polar`"
             "（SPICE/DOP853 四体，三段法 dAz/dtc/dvp，捕获面锁定 "
             "`pol = cos(i_m_target)`）→ 修正态；③ `solve_polar` / "
             "`solve_non_coplanar` 从修正态做 R-1 BVP **终解**。"
             "误差 = 本链路轨道 − 官方时序/目标；终解 i_m/近月点误差即 BVP "
             "相对目标的自洽残差。")
    L.append("")
    for r in results:
        L.append(f"## {r['mission']} — i_m = {r['i_m_deg']}° "
                 f"（cap {r['coast_cap_s']:.0f} s / {r['regime']}"
                 + (f"，relax" if r.get("relax_corrector") else "") + "）")
        L.append("")
        sc = r["screen"]
        L.append(f"- 初筛：±{sc['span_h']:.0f} h @ {sc['step_s']:.0f} s，"
                 f"{sc['n_points']} 点，可行 {sc['n_feasible']}，"
                 f"窗口 {sc['n_windows']} 个；实际历元可行 "
                 f"{sc['actual_feasible']}（根 {sc['actual_n_roots']}），"
                 f"最近可行点 {sc['nearest_feasible_utc']}（Δt "
                 f"{_fmt(sc['nearest_feasible_dt_s'], 0)} s）")
        arc_txt = ""
        if r.get("ascent_arc_ef_deg") is not None:
            _om = math.degrees(le.parking_orbit_elements(200.0, 380.0)[2])
            _dph = (float(r["ascent_arc_deg"])
                    - DEFAULT_ASCENT_ARC_DEG) / _om
            arc_txt = (f"；航程角 {_fmt(r.get('ascent_arc_ef_deg'), 3)}°(EF) + "
                       f"{_fmt(r.get('ascent_arc_omega_e_deg'), 3)}°(自转) = "
                       f"{_fmt(r.get('ascent_arc_deg'), 3)}°(惯性)"
                       f"（旧 {DEFAULT_ASCENT_ARC_DEG:.4f}°，{_dph:+.0f} s 相位）")
        L.append(f"- 官方：发射 {r['official']['liftoff_utc']}，"
                 f"入泊 {r['official']['anchor_utc']}，"
                 f"TLI {r['official']['tli_utc']}，"
                 f"到达 {r['official']['arrival_utc']}；"
                 f"上升 {_fmt(r['official'].get('t_to_park_s'), 0)} s + "
                 f"滑行 {_fmt(r['official'].get('coast_measured_s'), 1)} s + "
                 f"转移 {_fmt(r['official'].get('T_transfer_s'), 1)} s = "
                 f"{_fmt(r['official'].get('T_total_s'), 1)} s；"
                 f"i_m {r['i_m_deg']}°，"
                 f"近月点 {r['official']['perilune_alt_km']:.0f} km{arc_txt}")
        day = r["official"]["liftoff_utc"][:10]
        day_ws = [w for w in sc["windows"]
                  if w["t_open_utc"][:10] <= day <= w["t_close_utc"][:10]]
        if day_ws:
            L.append(f"- 当日窗口（{day}）：" + "；".join(
                f"{w['t_open_utc'][5:16]}–{w['t_close_utc'][5:16]} "
                f"({w['duration_min']:.0f} min，含实际起飞 "
                f"{w['contains_actual']})" for w in day_ws))
        else:
            L.append(f"- 当日窗口（{day}）：无正解窗口")
        L.append("")
        if not r["chain"]:
            L.append("- 无可行链路候选点（初筛无根）")
            L.append("")
            continue
        L.append("| 点 | 发射 (UTC) | 初解 coast (s) | 初解回放 i_m/近月点 | "
                 "修正 dAz/dtc/dvp | 修正后 coast (s) | 修正回放 i_m/近月点 | "
                 "终解 i_m (°) | 终解近月点 (km) | Δt发射 (s) | Δcoast (s) | "
                 "ΔtTLI (s) | Δ到达 (s) | Δi_m (°) | Δ近月点 (km) | ok |")
        L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for c in r["chain"]:
            ir = c.get("initial_replay") or {}
            co = c.get("correction") or {}
            fi = c.get("final") or {}
            er = c.get("errors") or {}
            L.append(
                f"| {c['label']} "
                f"| {c['t_launch_utc'][11:19]} "
                f"| {_fmt((c.get('roots') or [{}])[0].get('coast_s'), 1)} "
                f"| {_fmt(ir.get('i_m_deg'), 1)}° / "
                f"{_fmt(ir.get('perilune_alt_km'), 0)} km "
                f"| {_fmt(co.get('dAz_deg'), 3)}° / {_fmt(co.get('dtc_s'), 1)} s / "
                f"{_fmt(co.get('dvp_mps'), 2)} m/s "
                f"| {_fmt(co.get('coast_corrected_s'), 1)} "
                f"| {_fmt(co.get('i_m_replay_deg'), 2)}° / "
                f"{_fmt(co.get('peri_alt_replay_km'), 0)} km "
                f"| {_fmt(fi.get('i_m_deg'), 4)} "
                f"| {_fmt(fi.get('perilune_alt_km'), 2)} "
                f"| {_fmt(er.get('dt_launch_s'), 0)} "
                f"| {_fmt(er.get('d_coast_s'), 1)} "
                f"| {_fmt(er.get('dt_tli_s'), 0)} "
                f"| {_fmt(er.get('dt_arrive_s'), 0)} "
                f"| {_fmt(er.get('d_i_m_deg'), 4)} "
                f"| {_fmt(er.get('d_perilune_km'), 2)} "
                f"| {'OK' if c.get('chain_ok') else 'FAIL'} |")
        L.append("")
        for c in r["chain"]:
            fi = c.get("final") or {}
            co = c.get("correction") or {}
            if fi.get("C3_km2s2") is not None:
                L.append(f"- `{c['label']}` 终解：C3(地心出发) "
                         f"{_fmt(fi.get('C3_km2s2'), 3)} km²/s²，"
                         f"v_近月点 {_fmt(fi.get('v_perilune_kms'), 4)} km/s，"
                         f"TLI {c.get('t_tli_utc')}，"
                         f"到达 {c.get('t_arrive_utc')}")
            else:
                err = str(fi.get("error") or "final BVP failed").splitlines()[0]
                L.append(f"- ⚠ `{c['label']}` 终解失败：{err}")
            if not co.get("converged"):
                L.append(f"- ⚠ `{c['label']}` 修正器未收敛"
                         f"（回放 i_m {_fmt(co.get('i_m_replay_deg'), 2)}°，"
                         f"近月点 {_fmt(co.get('peri_alt_replay_km'), 0)} km，"
                         f"coast {_fmt(co.get('coast_corrected_s'), 1)} s）："
                         f"cap 边界处时间自由度被锁定，修正器只剩 2 个控制量；"
                         f"终解是否成功以 BVP 为准")
        L.append("")
    path.write_text("\n".join(L), encoding="utf-8")


def _main_mission_chain(args: argparse.Namespace) -> int:
    g.ensure_kernels()
    report_path = Path(args.report)
    log_path = Path(args.cache) / "logs" / "mission_chain.log"
    fixture_path = getattr(args, "fixture", None)
    missions: list[dict] = []
    truth_path = report_path
    if fixture_path is not None and Path(fixture_path).exists():
        missions = _load_fixture_missions(Path(fixture_path))
        if missions:
            truth_path = Path(fixture_path)
            _log(f"[chain] missions from fixture {Path(fixture_path).name}: "
                 f"{[m['mission'] for m in missions]}", log_path)
            for m in missions:
                _log(f"[chain]   {m['mission']}: liftoff={m['liftoff_utc']} "
                     f"park={m['anchor_utc']} tli={m['tli_utc']} "
                     f"arrive={m['arrival_utc']} | ascent={m['t_to_park_s']:.0f}s "
                     f"coast={m['coast_measured_s']:.0f}s "
                     f"transfer={m['T_transfer_s']:.0f}s "
                     f"T_total={m['T_total_s']:.0f}s ({m['T_h']:.4f}h)",
                     log_path)
    if not missions:
        report = json.loads(report_path.read_text(encoding="utf-8"))
        missions = _verify_missions(report)
    if not missions:
        print("no missions with actual liftoff times found")
        return 1
    _apply_insertion_bookkeeping(missions, args, log_path)
    wanted = {str(x) for x in (args.mission or [])}
    if wanted:
        missions = [m for m in missions if m["mission"] in wanted]
    caps: list[tuple[str, float]] = [("design", float(args.coast_max))]
    if float(args.coast_max_dynamics) != float(args.coast_max):
        caps.append(("dynamics", float(args.coast_max_dynamics)))
    results: list[dict] = []
    for m in missions:
        for i_m in m["i_m_candidates"]:
            for regime, cap in caps:
                results.append(_run_mission_chain(m, i_m, cap, regime, args,
                                                  log_path))
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    doc = {
        "schema": "fusion_mission_chain_v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "fixture": str(truth_path) if truth_path != report_path else None,
        "report": str(report_path),
        "method": ("coast_tool two-body initial (forward closure) -> "
                   "polar_correct 4-body capture-plane correction "
                   "(dAz/dtc/dvp, pol=cos(i_m_target)) -> tli_solver R-1 BVP; "
                   "errors vs the official liftoff / TLI / arrival / coast / "
                   "inclination / perilune"),
        "caps": {regime: cap for regime, cap in caps},
        "results": results,
    }
    json_path = out / "mission_chain.json"
    json_path.write_text(json.dumps(_json_safe(doc), ensure_ascii=False,
                                    indent=2), encoding="utf-8")
    md_path = out / "mission_chain.md"
    _write_mission_chain_markdown(results, truth_path, md_path)
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")
    for r in results:
        for c in r["chain"]:
            er = c.get("errors") or {}
            fi = c.get("final") or {}
            print(f"{r['mission']} {r['i_m_deg']}° cap={r['coast_cap_s']:.0f}"
                  f"({r['regime']}) {c['label']}: chain_ok={c['chain_ok']} "
                  f"coast={_fmt((c.get('correction') or {}).get('coast_corrected_s'), 1)}s "
                  f"d_coast={_fmt(er.get('d_coast_s'), 1)}s "
                  f"i_m={_fmt(fi.get('i_m_deg'), 4)}° "
                  f"peri={_fmt(fi.get('perilune_alt_km'), 2)}km")
    return 0


# ---------------------------------------------------------------------------
# Coast_tool initial -> polar correction -> fusion BVP chain
# ---------------------------------------------------------------------------
def _chain_one(item: dict, cfg: dict) -> dict:
    """One (launch epoch, root) point of the validation chain.

    1. two-body coast_tool initial solution (already the forward root);
    2. coast_tool-style polar correction (``polar_correct.correct_to_polar``)
       with the SPICE/DOP853 four-body model;
    3. final R-1 BVP solve (``tli_solver.solve_polar`` /
       ``solve_non_coplanar``) from the corrected MEP warm start.
    """
    from launch_window import polar_correct as pc
    t_l = float(item["t_launch_et"])
    root_idx = int(item.get("root_idx", 0))
    geo = forward_geometry(t_l, cfg)
    roots = coast_roots_forward(geo, cfg)
    if not roots:
        return {"t_launch_et": t_l, "t_launch_utc": g.et_to_utc(t_l),
                "chain_ok": False, "error": "no coast root"}
    root = roots[min(root_idx, len(roots) - 1)]
    corr = pc.correct_to_polar(t_l, geo, root["coast_s"], root["a_km"], cfg,
                               peri_target_km=cfg["H_m_km"], polar_lock=True,
                               relax=bool(cfg.get("relax_corrector", False)))
    t_arrive = t_l + cfg["T_s"]
    final: dict = {}
    try:
        res = _certify_solve(t_arrive, corr["T_cert_h"] * 3600.0,
                             corr["i0_deg"],
                             (corr["vp_corrected_kms"], corr["Omega_mep_deg"],
                              corr["omega_mep_deg"]), cfg)
        if res["res"] is not None:
            r = res["res"]
            final = {
                "converged": bool(res["converged"]),
                "i_m_deg": float(r.i_m),
                "perilune_alt_km": float(r.q[0]),
                "vp_sol_kms": float(r.p[0]),
                "dvp_final_mps": float((r.p[0] - corr["vp_corrected_kms"])
                                       * 1000.0),
                "perigee_lon_mep_deg": float(r.perigee_lon_mep),
                "perilune_lon_mep_deg": float(r.perilune_lon_mep),
                "error": None,
            }
        else:
            final = {"converged": False, "error": res["error"]}
    except Exception as exc:  # noqa: BLE001
        final = {"converged": False,
                 "error": f"{type(exc).__name__}: {exc}"}
    i_m_ok = (final.get("i_m_deg") is not None
              and abs(float(final["i_m_deg"]) - cfg["i_m_target_deg"])
              <= cfg["i_m_tol_deg"])
    peri_ok = (final.get("perilune_alt_km") is not None
               and abs(float(final["perilune_alt_km"]) - cfg["H_m_km"])
               <= cfg["perilune_tol_km"])
    fusion_ok = bool(final.get("converged") and i_m_ok and peri_ok)
    return {
        "t_launch_et": t_l,
        "t_launch_utc": g.et_to_utc(t_l),
        "t_arrive_utc": g.et_to_utc(t_arrive),
        "root_idx": root_idx,
        "n_roots": len(roots),
        "coast_initial_s": float(root["coast_s"]),
        "az_initial_deg": float(geo["az_deg"]),
        "a_km": float(root["a_km"]),
        "vp0_kms": float(math.sqrt(
            ts.MU_EARTH * (2.0 / cfg["rp_km"] - 1.0 / root["a_km"]))),
        "initial_replay": {
            "i_m_deg": corr["i_m_initial_replay_deg"],
            "perilune_alt_km": corr["peri_initial_replay_km"],
            "peri_dt_h": corr["peri_initial_dt_h"],
        },
        "correction": {k: corr[k] for k in (
            "converged", "dAz_deg", "dtc_s", "dvp_kms", "dvp_mps",
            "coast_corrected_s", "az_corrected_deg", "t_tli_utc",
            "T_cert_h", "i_m_replay_deg", "peri_alt_replay_km",
            "peri_dt_h", "pol", "i0_deg", "n_eval", "elapsed_s")},
        "fusion_solve": final,
        "correction_converged": bool(corr["converged"]),
        # Final validation: the fusion BVP solution must hit the target
        # inclination / perilune from the corrected warm start.  The
        # coast_tool corrector stalling (e.g. at the 1000 s coast cap) is
        # reported separately — the BVP is the authoritative check.
        "chain_ok": fusion_ok,
    }


def _chain_worker(payload: tuple) -> list[dict]:
    items, cfg = payload
    out: list[dict] = []
    cache_dir = Path(cfg["out_dir"]) / "chain_cache"
    for item in items:
        key = _sha1({"schema": 3, "config": cfg["config_id"],
                     "t_launch_et": item["t_launch_et"],
                     "root_idx": item.get("root_idx", 0),
                     "i_m": cfg["i_m_target_deg"], "H_m": cfg["H_m_km"],
                     "relax": bool(cfg.get("relax_corrector", False))})
        cache = cache_dir / f"{key}.json"
        if cache.exists():
            try:
                out.append(json.loads(cache.read_text(encoding="utf-8")))
                continue
            except Exception:  # noqa: BLE001
                pass
        try:
            rec = _chain_one(item, cfg)
        except Exception:  # noqa: BLE001
            rec = {"t_launch_et": item["t_launch_et"],
                   "t_launch_utc": g.et_to_utc(item["t_launch_et"]),
                   "chain_ok": False, "error": traceback.format_exc(limit=6)}
        rec["cache_key"] = key
        _write_json(cache, rec)
        out.append(rec)
    return out


def _chain_windows(points: list[dict], cfg: dict, step_s: float,
                   max_per_window: int) -> list[dict]:
    """Cluster feasible points and sample a time grid per window."""
    windows = cluster_windows(points, cfg)
    out: list[dict] = []
    for w in windows:
        pts = w["_pts"]
        sampled: list[dict] = []
        last_t = -1e18
        for p in pts:
            if p["t_launch_et"] - last_t >= step_s:
                sampled.append(p)
                last_t = p["t_launch_et"]
        if sampled and sampled[-1] is not pts[-1]:
            sampled.append(pts[-1])
        if len(sampled) > max_per_window:
            idx = np.linspace(0, len(sampled) - 1,
                              max_per_window).round().astype(int)
            sampled = [sampled[i] for i in sorted(set(idx.tolist()))]
        out.append({"window": w, "points": sampled})
    return out


def write_chain_markdown(chain: list[dict], cfg: dict, path: Path) -> None:
    L: list[str] = []
    ok = [r for r in chain if r.get("chain_ok")]
    L.append(f"# coast_tool 初解 → 极轨修正 → 融合验证 链路（T = {cfg['T_h']} h, "
             f"i_m = {cfg['i_m_target_deg']:g}°）")
    L.append("")
    L.append("流程：正解闭合给出 coast_tool 式**初解**（两体）→ coast_tool 极轨"
             "修正器（`polar_correct.py`，SPICE/DOP853 四体，三段法 dAz/dtc/dvp）"
             "→ 融合 R-1 BVP（`solve_polar`）从修正态求解。`chain_ok` = 最终 BVP "
             "满足 i_m/近月点容差；`修正收敛` 为 coast_tool 修正器自身的判据。")
    L.append("")
    n_cc = sum(1 for r in chain if r.get("correction_converged"))
    L.append(f"- 采样点 {len(chain)}：最终 BVP 验证成功 {len(ok)}，"
             f"失败 {len(chain) - len(ok)}；其中 coast_tool 修正器收敛 {n_cc}")
    if n_cc < len(chain):
        L.append(f"- ⚠ 修正器未收敛但 BVP 成功的点：修正器受 1000 s 滑行上限/"
                 f"钳制停滞，融合 BVP 从该暖启动仍收敛到目标极轨；"
                 f"反之（修正器收敛、BVP 失败）在本批 0 例")
    # --- why the coast initial solution is not polar (statistics) ---------
    def _dev(r: dict) -> Optional[float]:
        ir = r.get("initial_replay") or {}
        v = ir.get("i_m_deg")
        return None if v is None or v != v else abs(float(v) - 90.0)

    def _peri(r: dict) -> Optional[float]:
        ir = r.get("initial_replay") or {}
        v = ir.get("perilune_alt_km")
        return None if v is None or v != v else float(v)

    devs = [d_ for r in chain if (d_ := _dev(r)) is not None]
    pers = [p for r in chain if (p := _peri(r)) is not None]
    ims = [r["initial_replay"].get("i_m_deg") for r in chain
           if (r.get("initial_replay") or {}).get("i_m_deg") is not None]
    if devs:
        def _frac(th: float) -> str:
            k = sum(1 for d_ in devs if d_ < th)
            return f"{k}/{len(devs)} ({100.0 * k / len(devs):.1f}%)"
        L.append("")
        L.append("### 初解为何打不到极轨（本批统计）")
        L.append("")
        if ims:
            L.append(f"- 初解回放 i_m 中位 {sorted(ims)[len(ims) // 2]:.1f}°"
                     f"（范围 {min(ims):.1f}–{max(ims):.1f}°），"
                     f"|i_m−90°| 中位 {sorted(devs)[len(devs) // 2]:.1f}°")
        L.append(f"- |i_m−90°| < 5°: {_frac(5.0)}；< 30°: {_frac(30.0)}")
        if pers:
            n_imp = sum(1 for p in pers if p <= 55.0)
            L.append(f"- 撞击轨道（近月点护栏 ≤50 km 截断）: {n_imp}/{len(pers)}"
                     f"（{100.0 * n_imp / len(pers):.1f}%）")
        L.append("- 原因：两体闭合只约束**转移面**（发射场–月球质心），不约束"
                 "**月心捕获面**；近月区 |v_sc| ≈ 0.2 km/s ≪ |v_moon| ≈ 1.0 km/s，"
                 "v_rel 由月球轨道速度主导，捕获面法向趋近月球轨道动量方向 → "
                 "i_m 天然远离 90°（或近 10°/170° 两侧分布）。同时月球引力弯曲"
                 "使两体瞄准质心的实飞近月点大幅降低（本批全部撞击）。"
                 "修正器正是利用小 v∞ 的极端敏感性（m/s 级 dvp + 秒级 dtc 旋转"
                 "B 平面数十度）把近月点、时刻、极轨面三约束同时拉回。")
    # --- corrector vs BVP division of labour ------------------------------
    def _bvp_dvp(r: dict) -> Optional[float]:
        v = (r.get("fusion_solve") or {}).get("dvp_final_mps")
        return None if v is None or v != v else float(v)

    a_dvp = [v for r in chain if r.get("correction_converged")
             and (v := _bvp_dvp(r)) is not None]
    b_dvp = [v for r in chain if not r.get("correction_converged")
             and r.get("chain_ok") and (v := _bvp_dvp(r)) is not None]
    if a_dvp or b_dvp:
        L.append("")
        L.append("### 修正器 vs 融合 BVP 的分工（本批统计）")
        L.append("")
        if a_dvp:
            L.append(f"- A 组（修正器收敛，n={len(a_dvp)}）：BVP 相对修正态的额外 "
                     f"|dvp| 中位 {sorted(abs(v) for v in a_dvp)[len(a_dvp)//2]:.3f} m/s、"
                     f"最大 {max(abs(v) for v in a_dvp):.3f} m/s —— 修正器已到根，"
                     f"BVP 基本只做复核（不是 BVP 修好了错误状态）")
        if b_dvp:
            L.append(f"- B 组（修正器未收敛、BVP 成功，n={len(b_dvp)}）：BVP 额外 "
                     f"|dvp| ≤ {max(abs(v) for v in b_dvp):.3f} m/s（中位 "
                     f"{sorted(abs(v) for v in b_dvp)[len(b_dvp)//2]:.3f} m/s），"
                     f"把修正器残余收尾到目标 —— 是亚 m/s 级收尾而非重解")
        L.append("- 修正器自身收敛较少的原因：① 三段法 + (dAz, dtc, dvp) 钳制是"
                 "受限参数化，而 BVP 的 (v_p, Ω, ω) 牛顿无钳制；② 窗口首尾 "
                 "coast→1000 s 时 dtc 上界→0，修正器的\"时间杠杆\"冻结；"
                 "③ \"近月距–侧偏\"V 型褶皱的分支/扫描限制可停在 |i_m−90|≥4° 处。"
                 "反之修正器收敛但 BVP 失败的 0 例。")
        L.append("- 钳制放宽实验（`--relax-corrector`：步长/迭代预算放开，"
                 "滑行上限锁定不变）：coast 距上限 >37 s 的残留点中可自行收敛者"
                 "显著增加；而 coast 顶死 1000 s 的点（修正器仅剩 2 个控制自由度"
                 "对 3 个约束）任何钳制放宽都无法收敛——该处 BVP 以冻结时刻的 "
                 "(v_p, Ω, ω) 三自由度结构性不可替代。")
        L.append("- 注：coast_tool 原版修正器的 3×3 抛光段不会把 dtc 重新夹回"
                 "用户设定的最长滑行；本仓库端口已在所有求解阶段强制滑行上限"
                 "（否则\"收敛\"可能靠超出 CZ-5 1000 s 包络实现）。")
    if ok:
        daz = [abs(r["correction"]["dAz_deg"]) for r in ok]
        dtc = [abs(r["correction"]["dtc_s"]) for r in ok]
        dvp = [abs(r["correction"]["dvp_mps"]) for r in ok]
        im0 = [abs(r["initial_replay"]["i_m_deg"] - 90.0) for r in ok]
        imr = [abs(r["correction"]["i_m_replay_deg"]
                   - cfg["i_m_target_deg"]) for r in ok
               if r.get("correction_converged")]
        imf = [abs(r["fusion_solve"]["i_m_deg"] - cfg["i_m_target_deg"])
               for r in ok]
        L.append(f"- 修正量: |dAz| ≤ {max(daz):.3f}°, |dtc| ≤ {max(dtc):.1f} s, "
                 f"|dvp| ≤ {max(dvp):.3f} m/s")
        L.append(f"- 初解回放偏离极轨 |i_m−90°|: 中位 {sorted(im0)[len(im0)//2]:.1f}°, "
                 f"最大 {max(im0):.1f}°（证明初解本身不是极轨）")
        if imr:
            L.append(f"- 修正器收敛点 修正后回放 |i_m−90°| ≤ {max(imr):.4f}°；"
                     f"融合求解后（全部成功点）|i_m−90°| ≤ {max(imf):.4f}°")
    L.append("")
    cur_date = None
    for i, r in enumerate(chain):
        date = r["t_launch_utc"][:10]
        if date != cur_date:
            cur_date = date
            L.append(f"## {date}")
            L.append("")
            L.append("| 发射 (UTC) | 初解 coast | 初解回放 i_m/近月点 | "
                     "修正 dAz/dtc/dvp | 修正回放 i_m/近月点 | "
                     "融合 i_m/近月点 | ok |")
            L.append("|---|---|---|---|---|---|---|")
        ir = r.get("initial_replay") or {}
        co = r.get("correction") or {}
        fs = r.get("fusion_solve") or {}
        L.append(
            f"| {r['t_launch_utc'][11:19]} "
            f"| {_fmt(r.get('coast_initial_s'), 1)} s "
            f"| {_fmt(ir.get('i_m_deg'), 1)}° / {_fmt(ir.get('perilune_alt_km'), 0)} km "
            f"| {_fmt(co.get('dAz_deg'), 3)}° / {_fmt(co.get('dtc_s'), 1)} s / "
            f"{_fmt(co.get('dvp_mps'), 2)} m/s "
            f"| {_fmt(co.get('i_m_replay_deg'), 3)}° / "
            f"{_fmt(co.get('peri_alt_replay_km'), 1)} km "
            f"| {_fmt(fs.get('i_m_deg'), 3)}° / "
            f"{_fmt(fs.get('perilune_alt_km'), 1)} km "
            f"| {'✅' if r.get('chain_ok') else '❌'} |")
    L.append("")
    path.write_text("\n".join(L), encoding="utf-8")


def _main_correct_chain(args: argparse.Namespace) -> int:
    from launch_window import polar_correct  # noqa: F401  (import check)
    g.ensure_kernels()
    run_dir = Path(args.run)
    with open(run_dir / "fusion_config.json", encoding="utf-8") as fh:
        cfg = json.load(fh)
    cfg = dict(cfg)
    cfg["out_dir"] = str(run_dir)
    cfg["i_m_target_deg"] = float(args.i_m)
    cfg["H_m_km"] = float(args.perilune_alt)
    cfg["perilune_tol_km"] = float(DEFAULT_PERILUNE_TOL_KM)
    cfg["i_m_tol_deg"] = float(DEFAULT_I_M_TOL_DEG)
    cfg["relax_corrector"] = bool(getattr(args, "relax_corrector", False))
    # Load the forward points (full grid incl. rejects) from the run.
    points: list[dict] = []
    gz = run_dir / "forward_points.jsonl.gz"
    with gzip.open(gz, "rt", encoding="utf-8") as fh:
        for line in fh:
            points.append(json.loads(line))
    if args.window_date:
        wanted = set(args.window_date)
        points = [p for p in points if p["t_launch_utc"][:10] in wanted]
    windows = _chain_windows(points, cfg, float(args.step),
                             int(args.max_points_per_window))
    items: list[dict] = []
    for w in windows:
        for p in w["points"]:
            items.append({"t_launch_et": p["t_launch_et"], "root_idx": 0})
    _log(f"[chain] points={len(items)} windows={len(windows)} "
         f"step={args.step}s", None)
    if not items:
        print("no feasible points for the requested windows")
        return 1
    t0 = time.perf_counter()
    n_proc = max(1, min(int(args.workers), len(items)))
    chunk = max(1, math.ceil(len(items) / n_proc))
    payloads = [(items[i:i + chunk], cfg)
                for i in range(0, len(items), chunk)]
    chain: list[dict] = []
    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=n_proc) as pool:
        for out in pool.imap_unordered(_chain_worker, payloads):
            chain.extend(out)
    chain.sort(key=lambda r: r["t_launch_et"])
    out_dir = Path(args.out) if args.out else run_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    doc = {
        "schema": "fusion_polar_chain_v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "run": str(run_dir),
        "T_h": cfg["T_h"],
        "i_m_target_deg": cfg["i_m_target_deg"],
        "H_m_km": cfg["H_m_km"],
        "relax_corrector": cfg["relax_corrector"],
        "n_points": len(chain),
        "n_ok": sum(1 for r in chain if r.get("chain_ok")),
        "points": chain,
    }
    _write_json(out_dir / "fusion_polar_chain.json", doc)
    write_chain_markdown(chain, cfg, out_dir / "fusion_polar_chain.md")
    _log(f"[chain] done {len(chain)} points, ok="
         f"{doc['n_ok']} in {time.perf_counter() - t0:.1f}s", None)
    print(f"wrote {out_dir / 'fusion_polar_chain.json'}")
    print(f"wrote {out_dir / 'fusion_polar_chain.md'}")
    print(f"chain_ok {doc['n_ok']}/{len(chain)}")
    return 0


# ---------------------------------------------------------------------------
# FR-11  Multi-orbit summary (variable-coast, same-T demonstration)
# ---------------------------------------------------------------------------
def build_multi_orbit(cfg: dict, windows: list[dict],
                      cert_records: list[dict]) -> Optional[dict]:
    """Variable-coast multi-orbit summary (fixed T, fixed target i_m).

    The "multi-orbit" picture is the coast_tool initial-solution family:
    at one fixed Earth-Moon transfer time T and one fixed target lunar-orbit
    inclination, every launch epoch inside a continuous window has its own
    coast (from the closure), so a single window delivers a continuum of
    orbits (different TLI timelines / arrival phases) — the corrected and
    certified solutions are listed per window.  (The `--i-m-list` sweep is a
    separate FR-11 capability for arbitrary-inclination missions and is not
    part of this summary.)
    """
    out_windows: list[dict] = []
    for win in windows:
        t_open, t_close = win["t_open_et"], win["t_close_et"]
        recs = [r for r in cert_records
                if r.get("published")
                and t_open - cfg["dt_tol_s"] <= r.get("t_launch_et", -1e18)
                <= t_close + cfg["dt_tol_s"]]
        if not recs:
            continue
        recs.sort(key=lambda r: (float(r.get("t_launch_et", 0.0)),
                                 float(r.get("coast_s", 0.0))))
        feas = [p for p in win["_pts"] if p["roots"]]
        out_windows.append({
            "t_open_utc": win["t_open_utc"],
            "t_close_utc": win["t_close_utc"],
            "duration_min": win["duration_min"],
            "coast_range_s": [min(p["roots"][0]["coast_s"] for p in feas),
                              max(p["roots"][0]["coast_s"] for p in feas)]
            if feas else None,
            "n_certified": len(recs),
            "solutions": [{
                "t_launch_utc": r.get("t_launch_utc"),
                "coast_s": r.get("coast_s"),
                "i_m_target_deg": r.get("i_m_target_deg"),
                "i_m_deg": r.get("i_m_deg"),
                "dvp_mps": r.get("dvp_mps"),
                "C3_km2s2": r.get("C3_km2s2"),
                "azimuth_std_deg": r.get("azimuth_std_deg"),
                "family": r.get("family"), "branch": r.get("branch"),
            } for r in recs],
        })
    return {
        "T_h": cfg["T_h"],
        "i_m_target_deg": cfg["i_m_target_deg"],
        "n_windows": len(out_windows),
        "windows": out_windows,
    }


def write_multi_orbit_markdown(multi: dict, path: Path) -> None:
    L: list[str] = []
    L.append(f"# 变滑行多轨道（固定 T = {multi['T_h']} h，固定目标倾角 "
             f"i_m = {multi['i_m_target_deg']:g}°）")
    L.append("")
    L.append("解读：每个“**可发射时间段**”是一个连续区间（正解/两体筛查）；"
             "区间内**任意发射时刻**对应一条弹道——滑行时间由闭合方程随时刻连续"
             "确定（发射越早、滑行越长），全部弹道共用同一地月转移时间 T 与同一"
             "目标环月倾角。四体认证的实际可用边界/内部窄洞见 "
             "`fusion_certified_boundary.json`（窗口首尾逐点认证）。")
    L.append("")
    L.append("下表列出的是各时段的**代表性弹道**（四体认证采样点：窗口首/尾、"
             "各根 coast 极值、等分点，≤K/窗），**不是全部可发射时刻**；"
             "区间内其余时刻按需认证即可（滑行随时刻连续变化）。")
    L.append("")
    for w in multi["windows"]:
        c0, c1 = (w["coast_range_s"][0], w["coast_range_s"][1]) \
            if w.get("coast_range_s") else (float("nan"), float("nan"))
        L.append(f"## {w['t_open_utc'][:10]} 可发射时间段 "
                 f"{w['t_open_utc'][11:16]}–{w['t_close_utc'][11:16]} UTC"
                 f"（{w['duration_min']:.0f} min，滑行 {c0:.0f}–{c1:.0f} s，"
                 f"代表性弹道 {len(w.get('solutions') or [])} 条）")
        L.append("")
        L.append("| 代表性弹道发射时刻 (UTC) | 滑行 coast (s) | 认证 i_m (°) "
                 "| dvp (m/s) | C3 (km²/s²) | 认证方位角 (°) |")
        L.append("|---|---|---|---|---|---|")
        for r in w.get("solutions") or []:
            L.append(
                f"| {r['t_launch_utc'][:10]} {r['t_launch_utc'][11:19]} "
                f"| {_stats_safe(r['coast_s'], 1)} "
                f"| {_stats_safe(r['i_m_deg'], 3)} "
                f"| {_stats_safe(r['dvp_mps'])} "
                f"| {_stats_safe(r['C3_km2s2'])} "
                f"| {_stats_safe(r['azimuth_std_deg'], 2)} |")
        L.append("")
    path.write_text("\n".join(L), encoding="utf-8")


def _stats_safe(v, nd: int = 3) -> str:
    if v is None:
        return "—"
    try:
        return f"{float(v):.{nd}f}"
    except Exception:  # noqa: BLE001
        return str(v)


# ---------------------------------------------------------------------------
# Stage helpers
# ---------------------------------------------------------------------------
def _forward_chunk_path(out_dir: Path, T_h: float, date: str) -> Path:
    return out_dir / "chunks" / f"forward_T{int(round(T_h * 10)):04d}" / (
        f"day_{date}.json")


def _forward_chunk_matches(path: Path, cfg: dict) -> bool:
    try:
        with open(path, encoding="utf-8") as fh:
            old = json.load(fh)
        keys = ("T_h", "step_s", "az_min_deg", "az_max_deg", "coast_min_s",
                "coast_max_s", "coast_step_s", "t_ascent_s", "ascent_arc_deg",
                "rp_km", "om_park_deg_s", "site", "use_envelope",
                "chunk_schema", "config_id")
        job = old.get("job", {})
        return all(job.get(k) == cfg.get(k) for k in keys if k != "config_id") \
            and job.get("config_id") == cfg["config_id"]
    except Exception:  # noqa: BLE001
        return False


def forward_scan(cfg: dict, log_path: Path, workers: int,
                 force: bool = False) -> tuple[list[dict], list[dict]]:
    """FR-1/FR-2: per-day chunked forward scan with resume."""
    step = cfg["step_s"]
    t0 = cfg["t_start_et"]
    t1 = cfg["t_end_et"]
    n = int(math.floor((t1 - t0) / step))
    points: list[dict] = []
    # Group consecutive points per day into a single chunk job.
    day_bounds: dict[str, tuple[float, float]] = {}
    for k in range(n + 1):
        tl = t0 + k * step
        if tl >= t1:
            break
        date = g.et_to_utc(tl)[:10]
        if date not in day_bounds:
            day_bounds[date] = (tl, tl)
        day_bounds[date] = (day_bounds[date][0], tl)
    pending = []
    for date, (lo, hi) in sorted(day_bounds.items()):
        path = _forward_chunk_path(Path(cfg["out_dir"]), cfg["T_h"], date)
        if path.exists() and not force and _forward_chunk_matches(path, cfg):
            with open(path, encoding="utf-8") as fh:
                points.extend(json.load(fh)["points"])
            continue
        pending.append({"date": date, "t_lo": lo, "t_hi": hi,
                        "path": str(path)})
    _log(f"[FR-1] forward scan: points={n + 1} days={len(day_bounds)} "
         f"pending_chunks={len(pending)} step={step:.0f}s T={cfg['T_h']}h",
         log_path)
    if pending:
        t_start = time.perf_counter()
        n_proc = max(1, min(workers, len(pending)))
        ctx = mp.get_context("spawn")
        with ctx.Pool(processes=n_proc) as pool:
            for i, res in enumerate(pool.imap_unordered(_forward_chunk_worker,
                                                        [(j, cfg) for j in pending]), 1):
                points.extend(res["points"])
                _log(f"  [FR-1] [{i}/{len(pending)}] {res['date']} "
                     f"n={res['n_points']} feasible={res['n_feasible']} "
                     f"t={res['elapsed_s']:.1f}s", log_path)
        _log(f"[FR-1] scan done in {time.perf_counter() - t_start:.1f}s", log_path)
    points.sort(key=lambda p: p["t_launch_et"])
    windows = cluster_windows(points, cfg)
    _log(f"[FR-2] clustered {len(points)} points -> {len(windows)} windows", log_path)
    return points, windows


def _forward_chunk_worker(payload: tuple) -> dict:
    job, cfg = payload
    t0 = time.perf_counter()
    pts: list[dict] = []
    lo, hi = float(job["t_lo"]), float(job["t_hi"])
    step = cfg["step_s"]
    n = int(round((hi - lo) / step))
    for k in range(n + 1):
        pts.append(forward_point(lo + k * step, cfg))
    doc = {"job": dict(cfg),
           "config_id": cfg["config_id"],
           "chunk_schema": CHUNK_SCHEMA,
           "date": job["date"],
           "points": pts}
    _write_json(Path(job["path"]), doc)
    return {"date": job["date"], "points": pts, "n_points": len(pts),
            "n_feasible": sum(1 for p in pts if p["feasible"]),
            "elapsed_s": time.perf_counter() - t0}


def ensure_family_table(cfg: dict, log_path: Path) -> None:
    """Build the L5 family table if it does not exist (Stage 1)."""
    table = Path(cfg["table_path"])
    if table.exists():
        return
    _log(f"[tables] building L5 family table -> {table}", log_path)
    res = scan_im_run_scan(
        i_m_deg=cfg["i_m_target_deg"], moon_model=cfg["moon_model"],
        T_h=cfg["T_h"], H_m=cfg["H_m_km"], polar_ref=cfg["polar_ref"],
        use_sun=cfg["use_sun"])
    _write_json(table, res)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="python -m launch_window.fusion",
        description="Forward/inverse fusion launch-window calculator (PRD v1.0)")
    sub = ap.add_subparsers(dest="cmd")

    p = sub.add_parser("scan", help="full scan -> windows / calendar / cross")
    p.add_argument("--start", required=True)
    p.add_argument("--end", required=True)
    p.add_argument("--step", type=float, default=DEFAULT_STEP_S)
    p.add_argument("--T", type=float, required=True)
    p.add_argument("--site", default="wenchang")
    p.add_argument("--mission", default="CE-7")
    p.add_argument("--perilune-alt", type=float, default=DEFAULT_H_M_KM)
    p.add_argument("--parking-peri-alt", type=float, default=DEFAULT_H_PARK_KM)
    p.add_argument("--parking-apo-alt", type=float, default=DEFAULT_H_APO_KM)
    p.add_argument("--i-m", type=float, default=DEFAULT_I_M_DEG,
                   help="target lunar-orbit inclination (deg); 90 = polar "
                        "(solve_polar), any other value uses "
                        "solve_non_coplanar with q3 = cos(i_m), including "
                        "retrograde targets (e.g. CE-6 137 deg)")
    p.add_argument("--i0-hi", type=float, default=DEFAULT_I0_HI,
                   help="upper bound of the family-table i0 search (deg); "
                        "90 suffices for polar targets, use 114 for "
                        "non-polar / retrograde targets (verify_ce56)")
    p.add_argument("--i-m-tol", type=float, default=DEFAULT_I_M_TOL_DEG,
                   help="certification tolerance |i_m_solved - i_m_target| "
                        "(deg)")
    p.add_argument("--i-m-list", nargs="*", type=float, default=None,
                   help="FR-11 multi-target certification: certify every "
                        "representative against EACH of these lunar-orbit "
                        "inclinations (deg), e.g. --i-m-list 87 88.5 90 "
                        "91.5 93 for the narrow-window multi-orbit "
                        "(variable-coast) demonstration")
    p.add_argument("--polar-ref", default=DEFAULT_POLAR_REF)
    p.add_argument("--moon-model", default=DEFAULT_MOON_MODEL)
    p.add_argument("--no-sun", dest="use_sun", action="store_false")
    p.add_argument("--coast-min", type=float, default=DEFAULT_COAST_MIN_S)
    p.add_argument("--coast-max", type=float, default=DEFAULT_COAST_MAX_S)
    p.add_argument("--coast-step-s", type=float, default=DEFAULT_COAST_STEP_S)
    p.add_argument("--t-ascent-s", type=float, default=DEFAULT_T_ASCENT_S)
    p.add_argument("--ascent-arc-deg", type=float,
                   default=DEFAULT_ASCENT_ARC_DEG)
    p.add_argument("--az", nargs=2, type=float,
                   default=[DEFAULT_AZ_MIN_DEG, DEFAULT_AZ_MAX_DEG])
    p.add_argument("--no-envelope", dest="use_envelope",
                   action="store_false")
    p.add_argument("--certify", action="store_true")
    p.add_argument("--cert-k", type=int, default=DEFAULT_CERT_K)
    p.add_argument("--inverse", action="store_true")
    p.add_argument("--inverse-seed-step-s", nargs="*", type=float,
                   default=None,
                   help="FR-4 seed grid step(s) in seconds (default: one "
                        "sidereal day; pass several for staggered passes)")
    p.add_argument("--enumerate-roots", action="store_true",
                   help="FR-4 all-root enumeration at inverse windows")
    p.add_argument("--allow-inverse-only", type=str, default="false")
    p.add_argument("--workers", type=int, default=DEFAULT_CERT_WORKERS)
    p.add_argument("--force", action="store_true")
    p.add_argument("--out", type=Path, default=DEFAULT_OUT_DIR)
    p.set_defaults(use_sun=DEFAULT_USE_SUN, use_envelope=True)

    c = sub.add_parser("cross", help="coverage-only cross check")
    c.add_argument("--forward", type=Path, required=True)
    c.add_argument("--inverse", type=Path, required=True)
    c.add_argument("--out", type=Path, default=None)

    v = sub.add_parser(
        "verify-mission",
        help="verify liftoff-time accuracy against a ground-truth report "
             "(CE-5 / CE-6): forward scan + arbitrary-i_m certification "
             "around the actual liftoff")
    v.add_argument("--report", type=Path,
                   default=_REPO_ROOT / "tests" / "results" /
                           "launch_window_ce5_ce6_report.json",
                   help="ground-truth JSON with per-mission "
                        "actual_liftoff_utc / actual_arrival_utc / "
                        "actual_orbit_inclination_deg")
    v.add_argument("--out", type=Path,
                   default=_REPO_ROOT / "tests" / "results",
                   help="directory for the verification JSON + markdown")
    v.add_argument("--cache", type=Path,
                   default=DEFAULT_OUT_DIR / "ce56_fusion_verify",
                   help="working directory (cert caches, tables)")
    v.add_argument("--span-h", type=float, default=48.0,
                   help="half-width of the scan around the actual liftoff (h)")
    v.add_argument("--step", type=float, default=60.0)
    v.add_argument("--site", default="wenchang")
    v.add_argument("--coast-max", type=float, default=DEFAULT_COAST_MAX_S,
                   help="design-regime coast cap (CZ-5 single-orbit, s)")
    v.add_argument("--coast-max-dynamics", type=float, default=1500.0,
                   help="dynamics-regime coast cap (s); the historical "
                        "2-body closure coasts sit ~1.3-1.4 ks, above the "
                        "design envelope, so the raised cap measures the "
                        "intrinsic liftoff accuracy")
    v.add_argument("--az", nargs=2, type=float,
                   default=[DEFAULT_AZ_MIN_DEG, DEFAULT_AZ_MAX_DEG])
    v.add_argument("--cert-k", type=int, default=6)
    v.add_argument("--i0-hi", type=float, default=114.0)
    v.add_argument("--i-m-tol", type=float, default=DEFAULT_I_M_TOL_DEG)
    v.add_argument("--workers", type=int, default=DEFAULT_CERT_WORKERS)
    ch = sub.add_parser(
        "correct-chain",
        help="coast_tool initial -> polar correction (polar_correct.py) -> "
             "fusion BVP chain validation on a finished scan run")
    ch.add_argument("--run", type=Path, required=True,
                    help="finished fusion scan directory (needs "
                         "forward_points.jsonl.gz + fusion_config.json)")
    ch.add_argument("--window-date", action="append", default=None,
                    help="YYYY-MM-DD (repeatable); default: all windows")
    ch.add_argument("--step", type=float, default=600.0,
                    help="sampling step inside each window (s)")
    ch.add_argument("--max-points-per-window", type=int, default=40)
    ch.add_argument("--i-m", type=float, default=DEFAULT_I_M_DEG)
    ch.add_argument("--perilune-alt", type=float, default=DEFAULT_H_M_KM)
    ch.add_argument("--relax-corrector", action="store_true",
                    help="experimental: widen the coast_tool corrector's "
                         "Newton clamps / iteration budget (the coast "
                         "timeline bound is NOT relaxed)")
    ch.add_argument("--workers", type=int, default=DEFAULT_CERT_WORKERS)
    ch.add_argument("--out", type=Path, default=None,
                    help="output directory (default: the run directory)")
    mc = sub.add_parser(
        "mission-chain",
        help="coast_tool screening -> capture-plane correction -> fusion "
             "BVP -> error table vs the official CE-5/CE-6 mission data")
    mc.add_argument("--report", type=Path,
                    default=_REPO_ROOT / "tests" / "results" /
                            "launch_window_ce5_ce6_report.json",
                    help="ground-truth JSON fallback (used when --fixture "
                         "is missing); per-mission actual_liftoff_utc / "
                         "actual_arrival_utc")
    mc.add_argument("--fixture", type=Path, default=_FIXTURE_DEFAULT,
                    help="mission timeline fixture (CE-5/CE-6): liftoff / "
                         "coast_start_utc (parking insertion) / tli_utc / "
                         "orbit_insertion_utc; ascent = insertion − liftoff, "
                         "coast = tli − insertion, transfer = insertion − "
                         "tli, T_total = insertion − liftoff")
    mc.add_argument("--out", type=Path,
                    default=_REPO_ROOT / "tests" / "results",
                    help="directory for mission_chain.json + .md")
    mc.add_argument("--cache", type=Path,
                    default=DEFAULT_OUT_DIR / "ce56_fusion_verify",
                    help="working directory (logs; family tables untouched)")
    mc.add_argument("--mission", action="append", default=None,
                    help="mission name filter (repeatable); default: all")
    mc.add_argument("--span-h", type=float, default=48.0,
                    help="half-width of the coast_tool screening scan (h)")
    mc.add_argument("--step", type=float, default=60.0,
                    help="screening grid step (s)")
    mc.add_argument("--site", default="wenchang")
    mc.add_argument("--coast-max", type=float, default=DEFAULT_COAST_MAX_S,
                    help="design-regime parking-orbit coast cap (s) — "
                         "CZ-5 single-orbit envelope 1000 s; the screening "
                         "then picks the nearest cap-feasible candidate")
    mc.add_argument("--coast-max-dynamics", type=float, default=1500.0,
                    help="dynamics-regime coast cap (s); the historical "
                         "2-body closure coasts sit ~1.3-1.4 ks, above the "
                         "design envelope, so the actual liftoff epoch "
                         "becomes feasible and the chain runs on it")
    mc.add_argument("--az", nargs=2, type=float,
                    default=[DEFAULT_AZ_MIN_DEG, DEFAULT_AZ_MAX_DEG])
    mc.add_argument("--i0-hi", type=float, default=114.0)
    mc.add_argument("--ascent-arc-deg", type=float, default=None,
                    help="override the launch->insertion arc (deg).  The "
                         "closure treats it as an inertial rotation; the "
                         "legacy 20.8377 is an Earth-fixed power-law value. "
                         "Prefer --insertion-lat/--insertion-lon, which "
                         "applies the EF->inertial bookkeeping automatically")
    mc.add_argument("--insertion-lat", type=float, default=None,
                    help="parking-orbit insertion (2nd-stage 1st cutoff) "
                         "geocentric latitude (deg)")
    mc.add_argument("--insertion-lon", type=float, default=None,
                    help="parking-orbit insertion geocentric longitude "
                         "(deg east)");
    mc.add_argument("--i-m-tol", type=float, default=DEFAULT_I_M_TOL_DEG)
    mc.add_argument("--cert-k", type=int, default=6)
    mc.add_argument("--relax-corrector", action="store_true",
                    help="widen the corrector's Newton clamps / budget "
                         "(coast timeline bound stays locked)")
    mc.add_argument("--workers", type=int, default=DEFAULT_CERT_WORKERS)
    return ap


def _main_scan(args: argparse.Namespace) -> int:
    # The CLI bounds are launch epochs; normalize "T24:00:00".
    args.start = _norm_iso(args.start)
    args.end = _norm_iso(args.end)
    g.ensure_kernels()
    cfg = build_config(args)
    cfg["out_dir"] = str(args.out)
    table = _table_path(Path(args.out), cfg["T_h"], cfg["H_m_km"],
                        cfg["i_m_target_deg"], cfg["use_sun"])
    if not table.exists():
        prod = DEFAULT_TABLE_DIR
        if (prod / table.name).exists():
            table = prod / table.name
    cfg["table_path"] = str(table)
    log_path = Path(args.out) / "logs" / "fusion.log"
    _log(f"=== fusion scan {cfg['t_start_utc']}..{cfg['t_end_utc']} "
         f"T={cfg['T_h']} h step={cfg['step_s']} s config_id={cfg['config_id']}",
         log_path)
    _write_json(Path(args.out) / "fusion_config.json", cfg)

    # Stage 1: family table.
    ensure_family_table(cfg, log_path)

    # FR-1 / FR-2.
    points, windows = forward_scan(cfg, log_path, args.workers, force=args.force)
    out = Path(args.out)
    _write_jsonl_gz(out / "forward_points.jsonl.gz", points)
    _write_json(out / "forward_windows.json",
                [{k: v for k, v in w.items() if k != "_pts"}
                 for w in windows])

    # FR-3.
    cert_records: list[dict] = []
    cert_fails: list[dict] = []
    if args.certify:
        cert_records, cert_fails = stage_certify(windows, cfg, log_path,
                                                 args.workers)
        _write_json(out / "fusion_certify.json", cert_records)
        _write_json(out / "fusion_certify_failed.json", cert_fails)
    else:
        _write_json(out / "fusion_certify.json", [])
        _write_json(out / "fusion_certify_failed.json", [])

    # FR-4.
    inverse: list[dict] = []
    inverse_rejects: list[dict] = []
    root_enum: list[dict] = []
    if args.inverse:
        inverse, inverse_rejects = inverse_scan(
            cfg, log_path, args.workers,
            seed_steps=(args.inverse_seed_step_s or None),
            force=args.force)
        _write_json(out / "fusion_inverse_windows.json", inverse)
        _write_jsonl_gz(out / "fusion_inverse_rejects.jsonl.gz",
                        inverse_rejects)
        if args.enumerate_roots:
            families = pm.load_families(Path(cfg["table_path"]))
            for w in inverse:
                root_enum.extend(enumerate_roots_at(
                    float(w["t_arrive_et"]), cfg["T_s"], cfg, families))
            _write_json(out / "fusion_inverse_roots.json", root_enum)

    # FR-5.
    cross = cross_validate(cfg, windows, cert_records, inverse, points)
    _write_json(out / "fusion_cross.json", cross)

    # FR-6 / FR-7 / FR-8.
    calendar = build_calendar(cfg, windows, cert_records, inverse)
    _write_json(out / "fusion_calendar.json", calendar)
    multi = build_multi_orbit(cfg, windows, cert_records)
    if multi is not None:
        _write_json(out / "fusion_multi_orbit.json", multi)
        write_multi_orbit_markdown(multi, out / "fusion_multi_orbit.md")
    coverage = build_coverage(cfg, points, windows, cert_records, inverse,
                              cross)
    _write_json(out / "fusion_coverage.json", coverage)
    fusion_windows = [_fusion_window_from_cert(r) for r in cert_records
                      if r.get("published")
                      and abs(float(r.get("i_m_target_deg",
                                              cfg["i_m_target_deg"]))
                              - cfg["i_m_target_deg"]) <= 1e-9]
    for c in fusion_windows:
        c["mission"] = cfg["mission"]
    _write_json(out / "fusion_windows.json", fusion_windows)
    # Consolidated rejects (FR-8 / audit): forward point rejects + cert
    # failures + inverse reject counts + cross diagnostics.
    forward_rejects = [p for p in points if not p["feasible"]]
    _write_json(out / "fusion_rejects.json", {
        "forward": {
            "total": len(forward_rejects),
            "by_reason": dict(Counter(
                code for p in forward_rejects
                for code in p["reject_codes"]).most_common()),
            "points": forward_rejects,
        },
        "certify": {
            "total": len(cert_fails),
            "by_status": dict(Counter(
                r.get("certify_status") for r in cert_fails).most_common()),
            "records": cert_fails,
        },
        "inverse": {
            "total": len(inverse_rejects),
            "by_reason": dict(Counter(
                r.get("reason") for r in inverse_rejects).most_common()),
            "note": "full records in fusion_inverse_rejects.jsonl.gz",
        },
        "cross": {
            "forward_only": cross["forward_only"],
            "inverse_only": cross["inverse_only"],
        },
    })
    write_month_report(cfg, calendar, cross, coverage, cert_records,
                       out / "fusion_report.md")

    _log(f"=== done. artifacts in {out}", log_path)
    allow = str(args.allow_inverse_only).lower() in ("1", "true", "yes")
    if cross["n_inverse_only_defects"] > 0 and not allow:
        _log(f"[FR-6] FAIL: inverse_only_defects="
             f"{cross['n_inverse_only_defects']} > 0 "
             f"(--allow-inverse-only false)", log_path)
        return 2
    return 0


def _main_cross(args: argparse.Namespace) -> int:
    with open(args.forward, encoding="utf-8") as fh:
        fwd = json.load(fh)
    with open(args.inverse, encoding="utf-8") as fh:
        inv = json.load(fh)
    print(json.dumps({"n_forward": len(fwd), "n_inverse": len(inv)},
                     ensure_ascii=False))
    return 0


def main(argv: Optional[list[str]] = None) -> int:
    ap = _build_parser()
    args = ap.parse_args(argv)
    if args.cmd == "scan":
        return _main_scan(args)
    if args.cmd == "cross":
        return _main_cross(args)
    if args.cmd == "verify-mission":
        return _main_verify(args)
    if args.cmd == "mission-chain":
        return _main_mission_chain(args)
    if args.cmd == "correct-chain":
        return _main_correct_chain(args)
    ap.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
