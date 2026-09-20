#!/usr/bin/env python3
"""融合反向求解（fusion FR-4）后台作业 — 单次运行、输出 JSON。

由 ``viewer/server.py`` 以子进程方式调用（SPICE 非线程安全，长作业必须
隔离到独立进程，避免阻塞单线程 Web 服务）。输入参数已由 Web 层严格校验，
本脚本只做数值计算与结果落盘。

用法::

    python scripts/fusion_inverse_job.py \
        --start 2026-08-24 --end 2026-08-25 \
        --T 112.8 --site wenchang --coast-max 1200 \
        --out viewer/jobs/<job_id>

产物（``--out`` 目录下）：
- ``result.json``   作业结果（元信息 + 窗口列表 + 拒绝统计），Web 层直接读取
- ``fusion_config.json`` 完整求解配置（溯源）
- ``job.log``       stdout/stderr 由 Web 层重定向（本脚本同时打印）

支持配置：i_m = 90°（极轨捕获）、polar_ref = moon_equator、H_m = 200 km、
use_sun = True，转移时间 T 必须命中 ``viewer/data/fusion_tables/`` 中
预生成的 L5 家族表（112.0–113.0 h，0.1 h 步长，另含 112.61 h）。
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from launch_window import fusion as fu  # noqa: E402

DEFAULT_TABLE_DIR = Path(
    os.environ.get("FUSION_TABLE_DIR", str(REPO_ROOT / "viewer" / "data" / "fusion_tables"))
)


def _table_path(T_h: float, H_m: float, i_m_deg: float, use_sun: bool) -> Path:
    """按 fusion._table_path 命名规则在预生成表目录中定位家族表。"""
    name = fu._table_path(Path("."), T_h, H_m, i_m_deg, use_sun).name
    return DEFAULT_TABLE_DIR / name


def _build_cfg(args: argparse.Namespace) -> dict:
    ns = argparse.Namespace(
        site=args.site,
        T=args.T,
        parking_peri_alt=200.0,
        parking_apo_alt=380.0,
        az=[args.az_lo, args.az_hi],
        mission="CE-7",
        start=args.start,
        end=args.end,
        step=120.0,
        use_envelope=True,
        coast_min=0.0,
        coast_max=args.coast_max,
        coast_step_s=2.0,
        t_ascent_s=740.0,
        ascent_arc_deg=20.8377,
        i_m=90.0,
        i_m_list=[90.0],
        perilune_alt=200.0,
        polar_ref="moon_equator",
        moon_model="spice",
        use_sun=True,
        i0_hi=90.0,
        cert_k=8,
        allow_inverse_only=False,
    )
    cfg = fu.build_config(ns)
    out_dir = Path(args.out)
    cfg["out_dir"] = str(out_dir)
    cfg["table_path"] = str(_table_path(cfg["T_h"], cfg["H_m_km"], 90.0, True))
    return cfg


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="融合反向求解后台作业（FR-4）")
    ap.add_argument("--start", required=True, help="发射窗口起始日期 UTC（YYYY-MM-DD）")
    ap.add_argument("--end", required=True, help="发射窗口结束日期 UTC（YYYY-MM-DD，右开区间）")
    ap.add_argument("--T", type=float, required=True, help="转移时间（小时）")
    ap.add_argument("--site", default="wenchang")
    ap.add_argument("--coast-max", type=float, default=1200.0, help="停泊滑行上限（秒）")
    ap.add_argument("--az-lo", type=float, default=90.0)
    ap.add_argument("--az-hi", type=float, default=175.0)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)

    out_dir: Path = args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "job.log"

    table = _table_path(args.T, 200.0, 90.0, True)
    if not table.exists():
        print(
            f"[fusion] 不支持的转移时间 T={args.T} h（缺少家族表 {table.name}）",
            file=sys.stderr,
        )
        return 2

    t0 = time.perf_counter()
    cfg = _build_cfg(args)
    print(
        f"[fusion] start={args.start} end={args.end} T={cfg['T_h']} h "
        f"site={args.site} coast_max={args.coast_max}s table={table.name}",
        flush=True,
    )
    with open(out_dir / "fusion_config.json", "w", encoding="utf-8") as f:
        json.dump(fu._json_safe(cfg), f, ensure_ascii=False, indent=1)

    windows, rejects = fu.inverse_scan(cfg, log_path, max_workers=1)
    elapsed = time.perf_counter() - t0

    from collections import Counter

    reason_counts = Counter(r.get("reason") for r in rejects)
    result = {
        "schema": 1,
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "elapsed_s": round(elapsed, 1),
        "config": {
            "start": args.start,
            "end": args.end,
            "T_h": cfg["T_h"],
            "site": args.site,
            "coast_max_s": args.coast_max,
            "i_m_deg": 90.0,
            "polar_ref": cfg["polar_ref"],
            "moon_model": cfg["moon_model"],
            "use_sun": cfg["use_sun"],
            "table": table.name,
        },
        "n_windows": len(windows),
        "windows": windows,
        "rejects_total": len(rejects),
        "rejects_by_reason": dict(reason_counts.most_common()),
    }
    with open(out_dir / "result.json", "w", encoding="utf-8") as f:
        json.dump(fu._json_safe(result), f, ensure_ascii=False)
    print(
        f"[fusion] done: windows={len(windows)} rejects={len(rejects)} "
        f"elapsed={elapsed:.1f}s",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
