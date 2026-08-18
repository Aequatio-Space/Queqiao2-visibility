#!/usr/bin/env python3
"""P0 科学分析一键执行入口。
按依赖顺序依次运行 P0 分析，并生成汇总 ``p0_summary.json``：

1. ``earth_elevation.py``   — P0-2 地球高度角（先生成，供 P0-1 图表对比线使用）
2. ``relay_elevation.py``   — P0-3b 中继卫星三站点高度角/方位角时间序列
   （先生成，供 P0-1 图表叠加逐方位角高度角带）
3. ``horizon_mask.py``      — P0-1 坑沿遮挡角（叠加地球高度角带 + 中继高度角带）
4. ``queqiao_relay.py``     — P0-3 中继窗口估算 + 条件表述
5. 汇总 ``p0_summary.json`` — 三项结论 + 第一幕论点判断

每个脚本以独立子进程运行（``subprocess.run``），各自加载/释放 SPICE 内核，
避免内核状态串扰。任一脚本失败则中止并返回非零退出码。

用法::

    python analysis/elevation_video/run_all.py [--output-dir figures-and-logs]
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime, timezone

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from analysis.elevation_video.utils import (  # noqa: E402
    DEM_PATH,
    SHACKLETON_LAT_DEG,
    SHACKLETON_LON_DEG,
)

DEFAULT_OUTPUT_DIR = os.path.join(_REPO_ROOT, "figures-and-logs")
HORIZON_JSON = "horizon_mask_center.json"
EARTH_JSON = "earth_elevation_timeseries.json"
QUEQIAO_JSON = "queqiao_relay_analysis.json"

# 每步的默认 CLI 参数
STEPS: tuple[dict, ...] = (
    {
        "name": "P0-2 earth_elevation",
        "script": "earth_elevation.py",
        "args": [
            "--lat", str(SHACKLETON_LAT_DEG), "--lon", str(SHACKLETON_LON_DEG),
            "--h", "0.0",
            "--start", "2026-08-01", "--end", "2027-12-31",
            "--interval-h", "1.0",
            "--theoretical-extent-years", "18.6",
        ],
    },
    {
        "name": "P0-3b relay_elevation",
        "script": "relay_elevation.py",
        "args": [
            "--lat", str(SHACKLETON_LAT_DEG), "--lon", str(SHACKLETON_LON_DEG),
            "--h", "0.0",
            "--start", "2026-08-01", "--end", "2027-12-31",
            "--interval-h", "0.5",
            # 地形遮挡判定：高度角须大于该方位角地形地平线（读取 P0-1
            # horizon_mask_{suffix}.json；缺失时回退用 DEM 计算）
            "--terrain",
        ],
    },
    {
        "name": "P0-1 horizon_mask",
        "script": "horizon_mask.py",
        "args": [
            "--dem", DEM_PATH,
            "--lat", str(SHACKLETON_LAT_DEG), "--lon", str(SHACKLETON_LON_DEG),
            "--max-radius-m", "50000",
            "--n-az", "360", "--n-samples", "300",
        ],
    },
    {
        "name": "P0-3 queqiao_relay",
        "script": "queqiao_relay.py",
        "args": [
            "--lat", str(SHACKLETON_LAT_DEG), "--lon", str(SHACKLETON_LON_DEG),
            "--h", "0.0",
            "--start", "2026-08-01", "--end", "2027-12-31",
            "--interval-h", "0.5",
        ],
    },
)


def _run_step(step: dict, output_dir: str, python: str) -> None:
    """以子进程运行单个分析脚本，输出转发到 stdout。"""
    script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), step["script"])
    cmd = [python, script_path, "--output-dir", output_dir, *step["args"]]
    print(f"\n=== [{step['name']}] {step['script']} ===")
    print(f"    cmd: {' '.join(cmd)}")
    proc = subprocess.run(cmd, cwd=_REPO_ROOT)
    if proc.returncode != 0:
        print(f"[run_all] ERROR: step {step['name']} failed "
              f"(exit {proc.returncode})")
        sys.exit(proc.returncode or 1)


def _load_json(output_dir: str, fname: str) -> dict | None:
    path = os.path.join(output_dir, fname)
    if not os.path.exists(path):
        print(f"[run_all] WARNING: {fname} not found; summary entry omitted")
        return None
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def build_summary(output_dir: str) -> dict:
    """组装 p0_summary.json。"""
    hor = _load_json(output_dir, HORIZON_JSON)
    earth = _load_json(output_dir, EARTH_JSON)
    queqiao = _load_json(output_dir, QUEQIAO_JSON)

    p0_1: dict = {}
    if hor is not None:
        p0_1 = {
            "max_occlusion_deg": hor["max_elev_deg"],
            "mean_occlusion_deg": hor["mean_elev_deg"],
            "azimuth_of_max_deg": hor["azimuth_of_max_deg"],
            "evidence_level": hor.get("evidence_level", "A"),
        }

    p0_2: dict = {}
    if earth is not None:
        cur = earth.get("current_period", {})
        theo = earth.get("theoretical_extreme", {})
        p0_2 = {
            "current_max_elev_deg": cur.get("max_elev_deg"),
            "current_min_elev_deg": cur.get("min_elev_deg"),
            "theoretical_max_elev_deg": theo.get("max_elev_deg"),
            "evidence_level": earth.get("evidence_level", "A"),
        }

    p0_3: dict = {}
    if queqiao is not None:
        op = queqiao.get("orbit_parameters", {})
        vs = queqiao.get("visibility_stats", {})
        p0_3 = {
            "data_level": queqiao.get("data_level", "A"),
            "orbit_period_h": op.get("period_h"),
            "longest_relay_arc_h": vs.get("max_daily_visible_h"),
            "mean_daily_visible_h": vs.get("mean_daily_visible_h"),
            "conditional_statement": queqiao.get("conditional_statement", ""),
            "evidence_level": queqiao.get("evidence_level", "A"),
        }

    # P0-3b：中继卫星三站点高度角统计（relay_elevation.py 产物）
    p0_3b: dict = {}
    relay_stats = []
    for suffix in ("center", "north1km", "south1km"):
        rel = _load_json(output_dir, f"relay_elevation_{suffix}.json")
        if rel is None:
            continue
        st = rel.get("site_stats", {})
        site = rel.get("site", {})
        terrain_meta = rel.get("terrain", {})
        relay_stats.append(
            {
                "site_name": site.get("site_name"),
                "lat_deg": site.get("lat_deg"),
                "lon_deg": site.get("lon_deg"),
                "min_elev_deg": st.get("min_elev_deg"),
                "max_elev_deg": st.get("max_elev_deg"),
                "mean_elev_deg": st.get("mean_elev_deg"),
                "visible_fraction": st.get("visible_fraction"),
                "visible_terrain_fraction": st.get("visible_terrain_fraction"),
                "longest_relay_arc_h": st.get("longest_relay_arc_h"),
                "mean_daily_visible_h": st.get("mean_daily_visible_h"),
                "terrain_available": terrain_meta.get("available"),
                "terrain_source": terrain_meta.get("source"),
            }
        )
    if relay_stats:
        p0_3b = {
            "sites": relay_stats,
            "evidence_level": "A",
            "note": (
                "中继卫星在任务时段内的高度角范围（轨道相位依赖 Ω/M 假设值）。"
                "visible_fraction 为纯几何可见比例（水平线之上即算）；"
                "visible_terrain_fraction 为考虑地形遮挡的真实可见比例"
                "（高度角须大于该方位角地形地平线）"
            ),
        }

    summary = {
        "p0_1_horizon": p0_1,
        "p0_2_earth": p0_2,
        "p0_3_queqiao": p0_3,
        "p0_3b_relay_elevation": p0_3b,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }
    return summary


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR,
                   help="输出目录，默认 figures-and-logs")
    p.add_argument("--python", default=sys.executable,
                   help="Python 解释器路径（默认当前解释器）")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    if not os.path.exists(DEM_PATH):
        print(f"[run_all] ERROR: DEM not found: {DEM_PATH}")
        sys.exit(2)

    for step in STEPS:
        _run_step(step, args.output_dir, args.python)

    print(f"\n=== [run_all] generating {os.path.join(args.output_dir, 'p0_summary.json')} ===")
    summary = build_summary(args.output_dir)
    summary_path = os.path.join(args.output_dir, "p0_summary.json")
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print("[run_all] ALL P0 ANALYSES COMPLETE")


if __name__ == "__main__":
    main()
