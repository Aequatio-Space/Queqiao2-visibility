"""Vercel Serverless 入口 — 将 Flask app 包装为 ASGI/WSGI 处理器。

Vercel Python Runtime 通过 ``vercel.json`` 指向本文件。DEM 使用仓库内置的
压缩子集（``ldem_87s_5mpp_shackleton_50km_25m.tif``，52 MB，随仓库分发），
无需持久盘，所有端点（含 DEM 类：/api/dem_tile、/api/horizon_mask、
/api/terrain_gltf）均可直接工作。

访问统计（/api/visit、/api/click、/api/stats）可工作，但 CSV 写入实例的
临时文件系统，实例回收后数据丢失；如需持久统计，请接入外部存储
（Vercel Blob / KV / Postgres）改写 server.py 的 _append_csv/_read_csv。
"""

import os
import sys

# 仓库根加入 sys.path（本文件位于 api/index.py → 上溯一级即仓库根）
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from viewer.server import app  # noqa: E402

# Vercel Python Runtime 期望的 WSGI 应用名
app = app

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=int(os.environ.get("PORT", "8000")), debug=False)