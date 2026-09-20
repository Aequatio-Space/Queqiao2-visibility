#!/usr/bin/env python3
"""重新生成杨维廉论文 24 条转移轨道的 XYZ JSON 资产。

产物被 ``viewer/server.py`` 的 ``/api/yang_orbits`` 与
``/api/yang_orbits_stream`` 端点直接读取。由于一次完整构建（24 条 R-1
边值问题 + DOP853 采样）约需 20 s 单核 CPU，部署时不在请求内实时计算，
而是预生成静态 JSON 资产（安全性：消除慢端点 DoS 面；稳定性：不受
Serverless 超时限制）。

用法::

    python scripts/build_yang_orbits.py                 # 写入默认路径
    python scripts/build_yang_orbits.py --out /tmp/o.json

任一条轨道求解失败时返回非零退出码且不写文件，保证部署资产完整。
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from launch_window.orbits_xyz import (  # noqa: E402
    DEFAULT_N_SAMPLES,
    H_M_KM,
    H_P_KM,
    T_ARRIVE_UTC,
    build_yang_16_orbits_xyz,
    write_records_json,
)

DEFAULT_OUT = REPO_ROOT / "viewer" / "data" / "yang_orbits_xyz.json"


def main() -> int:
    ap = argparse.ArgumentParser(description="生成 Yang 论文轨道 XYZ JSON 资产")
    ap.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_OUT,
        help=f"输出路径（默认 {DEFAULT_OUT.relative_to(REPO_ROOT)}）",
    )
    args = ap.parse_args()

    records = build_yang_16_orbits_xyz()
    failed = [r for r in records if r.error]
    if failed:
        for r in failed:
            print(f"[FAIL] {r.family} #{r.orbit}: {r.error}", file=sys.stderr)
        print(f"{len(failed)}/{len(records)} 条轨道求解失败，未写出资产", file=sys.stderr)
        return 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_records_json(records, args.out)
    print(
        f"OK: {len(records)} 条轨道 · t_arrive={T_ARRIVE_UTC} · "
        f"H_m={H_M_KM} km · h_p={H_P_KM} km · {DEFAULT_N_SAMPLES} 采样/条"
    )
    print(f"Wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
