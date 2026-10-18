#!/usr/bin/env python3
"""交互式地平线掩码可视化工具 — Flask 后端。

端点一览：
- ``GET /``                         → 站点首页（四入口导航）
- ``GET /visibility``               → 地平线掩码可视化前端（原 ``/``）
- ``GET /orbits``                   → 杨维廉论文轨道 3D 可视化前端
- ``GET /coast_tool/``              → 月球发射窗口计算器（停泊轨道滑行方案）
- ``GET /elements_demo``            → 轨道六根数交互演示（二体/摄动/CR3BP）
- ``GET /assets/<path>``            → elements_demo 贴图资产（仓库根 assets/）
- ``GET /api/dem_tile``             → DEM 瓦片 PNG（20 km × 20 km 裁剪，附元数据头）
- ``GET /api/horizon_mask``         → 360° 地形遮挡角数组（lru_cache 缓存）
- ``GET /api/earth_elevation``      → 单时刻地球高度角 + 方位角
- ``GET /api/earth_daily``          → 一天 288 点地球高度角/方位角日序列
- ``GET /api/relay_elevation``      → 单时刻中继卫星高度角 + 方位角
- ``GET /api/relay_daily``          → 一天 288 点中继卫星高度角/方位角日序列
- ``GET /api/stereo_to_lonlat``     → 极射投影 (x, y) → 经纬度
- ``GET /api/yang_orbits``          → Yang 论文 24 条转移轨道 XYZ JSON（预生成资产）
- ``GET /api/yang_orbits_stream``   → 同上，SSE 逐轨道推送（meta/orbit/done）
- ``POST /api/visit``               → 记录一次页面访问（CSV 持久化）
- ``POST /api/click``               → 记录一次按钮点击（CSV 持久化）
- ``GET /api/stats``                → 按日期聚合的访问量 + 按钮点击量
- ``GET /api/health``               → 服务状态（供 smoke 测试）

环境变量：
- ``DEM_PATH``          DEM GeoTIFF 路径（必填，无默认值；未设置时 DEM 端点返回 503）
- ``SPICE_KERNEL_DIR``  SPICE 内核目录（默认仓库根，见 analysis/elevation_video/utils.py）
- ``PORT``/``VIEWER_PORT``  监听端口（默认 5000）
- ``HOST``              监听地址（默认 0.0.0.0；本地调试可设 127.0.0.1）
- ``VIEWER_STATS_DIR``  访问统计 CSV 落盘目录（默认 viewer/stats，部署时建议指向持久盘）

启动：:

    python viewer/server.py
    # 生产 WSGI（Render/Vercel 等）：
    #   gunicorn "viewer.server:app" --bind 0.0.0.0:$PORT --workers 1 --threads 1
    #   waitress-serve --host 0.0.0.0 --port $PORT viewer.server:app

SPICE 内核在首个需要天历的请求时**惰性加载**（``_ensure_spice``），
保证 WSGI 服务器（不经过 ``__main__``）也能正常工作。
"""

from __future__ import annotations

import csv
import io
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time
import uuid
from datetime import datetime, timezone
from functools import lru_cache

import numpy as np

# 允许以脚本方式从仓库任意位置运行：将仓库根加入 sys.path
# （本文件位于 viewer/server.py → 上溯两级即仓库根）
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors  # noqa: E402
import rasterio  # noqa: E402
from flask import (  # noqa: E402
    Flask,
    Response,
    jsonify,
    render_template,
    request,
    send_file,
    send_from_directory,
    stream_with_context,
)
from werkzeug.exceptions import BadRequest, HTTPException  # noqa: E402
from rasterio.enums import Resampling  # noqa: E402
from rasterio.windows import from_bounds  # noqa: E402

from analysis.elevation_video import utils as u  # noqa: E402
from analysis.elevation_video.horizon_mask import compute_horizon_mask  # noqa: E402
from analysis.elevation_video.queqiao_relay import build_orbital_state, propagate_satellite  # noqa: E402

# ---------------------------------------------------------------------------
# 路径与默认值
# ---------------------------------------------------------------------------
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DEM_PATH = u.DEM_PATH

DEFAULT_LAT_DEG = u.SHACKLETON_LAT_DEG  # -89.67
DEFAULT_LON_DEG = u.SHACKLETON_LON_DEG  # 129.78
SHACKLETON_LAT_DEG = u.SHACKLETON_LAT_DEG  # -89.67
SHACKLETON_LON_DEG = u.SHACKLETON_LON_DEG  # 129.78
DEFAULT_DATE = "2026-08-15"

# 多入口站点资产路径
# - coast_tool：月球发射窗口计算器（纯静态，脚本已外置以符合 CSP）
# - yang_orbits：Yang 论文 24 条转移轨道 XYZ JSON（scripts/build_yang_orbits.py
#   预生成；不在请求内实时求解，避免单线程服务被 ~20 s 慢计算拖垮。
#   存放于 viewer/data/ 而非 static/，强制经 API 读取，保持单一入口）
COAST_TOOL_DIR = os.path.join(_REPO_ROOT, "coast_tool")
# elements_demo 贴图（earth_day/earth_night/moon_surface），页内以 ../../assets/ 引用
ASSETS_DIR = os.path.join(_REPO_ROOT, "assets")
YANG_ORBITS_JSON = os.path.join(BASE_DIR, "data", "yang_orbits_xyz.json")
YANG_ORBITS_DEFAULT_T_ARRIVE = "2026-08-24T00:00:00Z"

# 融合反向求解（fusion FR-4）后台作业：
# 单日约 10 s、整月数分钟的 CPU 密集计算。SPICE 非线程安全，若在 Web
# 进程内运行会阻塞单线程服务，故以子进程运行 scripts/fusion_inverse_job.py，
# Web 层只做参数校验、作业调度与结果读取（子进程隔离 + 超时看门狗）。
FUSION_JOB_DIR = os.environ.get("FUSION_JOB_DIR", os.path.join(BASE_DIR, "jobs"))
FUSION_TABLE_DIR = os.environ.get(
    "FUSION_TABLE_DIR", os.path.join(BASE_DIR, "data", "fusion_tables")
)
FUSION_SCRIPT = os.path.join(_REPO_ROOT, "scripts", "fusion_inverse_job.py")
FUSION_MAX_SPAN_DAYS = 31
FUSION_JOB_TIMEOUT_S = int(os.environ.get("FUSION_JOB_TIMEOUT_S", "3600"))
FUSION_EST_SEC_PER_DAY = 15.0
FUSION_SITES = ("wenchang", "xichang", "ksc")
FUSION_I_M_DEG = 90.0

# 极射投影 → 经纬度使用月球赤道半径（与 utils.lonlat_to_polar_stereo_xy 一致）
MOON_R_M = u.MOON_R_EQ_M

# 日序列采样：每 5 分钟一个点 → 24h × 12 = 288 点
DAILY_N_POINTS = 288
DAILY_INTERVAL_MIN = 5

# 轨道状态锚定历元：与 relay_elevation.py / queqiao_relay.py 一致 ——
# 在**任务历元**（默认 2026-08-01 00:00 UTC）用假设平近点角 M=215° 构造
# 轨道初值，然后以月心二体连续传播到任意时刻。
MISSION_EPOCH_UTC = "2026-08-01T00:00:00"

# 单次渲染 DEM 瓦片尺寸（降采样目标，800px 仍远小于源 5m/px 全分辨率，
# 传输量 ~2.4 MB PNG，lru_cache 缓存后仅首帧慢）。
TILE_PX = 800

# 3D 第一人称天空视角地形（/api/terrain_gltf → GLTF → three.js 渲染）
TERRAIN_HALF_M = 15000.0   # 半幅（m）→ 30 km × 30 km 地形块（覆盖 20 km 地图 + 余量）
TERRAIN_GRID = 288         # 网格分辨率（288×288 ≈ 83k 顶点，~4 MB GLTF）
EYE_HEIGHT_M = 2.0         # 第一人称视角眼高（m）

app = Flask(__name__)
app.config["TEMPLATES_AUTO_RELOAD"] = bool(os.environ.get("FLASK_DEBUG"))
app.config["MAX_CONTENT_LENGTH"] = 1 * 1024 * 1024  # 请求体上限 1 MB（本服务仅收小 JSON）
app.json.ensure_ascii = False  # 中文 API 响应不做 \uXXXX 转义

# 隐藏 Werkzeug 开发服务器的 Server 头指纹。
# 开发服务器在 send_response() 时会自动写入 "Server: Werkzeug/x.y Python/z"，
# 与 after_request 设置的值重复。覆写 version_string 使其不暴露版本号；
# after_request 再统一设为 "queqiao2"。
try:
    import werkzeug.serving

    werkzeug.serving.WSGIRequestHandler.server_version = "queqiao2"
    werkzeug.serving.WSGIRequestHandler.sys_version = ""
except Exception:  # pragma: no cover - 非 Werkzeug 部署时无影响
    pass


# ---------------------------------------------------------------------------
# 安全响应头：after_request 钩子统一注入
# ---------------------------------------------------------------------------
@app.after_request
def _security_headers(resp):
    """为所有响应注入安全头并隐藏服务器指纹。

    - 隐藏 ``Server`` 头中的 Werkzeug/Python 版本，防止 CVE 精准匹配；
    - ``X-Frame-Options: DENY`` + CSP ``frame-ancestors 'none'`` → 防点击劫持；
    - ``X-Content-Type-Options: nosniff`` → 防 MIME 嗅探；
    - ``Strict-Transport-Security`` → 强制 HTTPS（防降级中间人）；
    - ``Referrer-Policy`` / ``Permissions-Policy`` → 最小化信息与权限暴露。
    """
    resp.headers["Server"] = "queqiao2"  # 隐藏 Werkzeug/Python 版本
    resp.headers["X-Content-Type-Options"] = "nosniff"
    resp.headers["X-Frame-Options"] = "DENY"
    resp.headers["Referrer-Policy"] = "strict-origin-when-cross-origin"
    resp.headers["Permissions-Policy"] = "geolocation=(), microphone=(), camera=()"
    resp.headers["Strict-Transport-Security"] = (
        "max-age=31536000; includeSubDomains"
    )
    # connect-src 需允许 data: —— 3D 地形 GLTF（/api/terrain_gltf）由 PyVista
    # export_gltf(inline_data=True) 生成，四个几何 buffer 以 data: URI 内联在
    # 响应体中；three.js GLTFLoader 通过 fetch(data:...) 加载 buffer，
    # 若 connect-src 仅 'self' 会拦截 data: URI，导致 3D 地形无法渲染。
    # data: URI 源自同源响应体（非外部资源），放行不扩大攻击面。
    resp.headers["Content-Security-Policy"] = (
        "default-src 'self'; script-src 'self'; "
        "style-src 'self' 'unsafe-inline'; "
        "img-src 'self' data: blob:; connect-src 'self' data:; "
        "object-src 'none'; base-uri 'none'; form-action 'self'; "
        "frame-ancestors 'none'"
    )
    return resp


# ---------------------------------------------------------------------------
# CSRF 防护：POST 请求校验 Origin/Referer 同源
# ---------------------------------------------------------------------------
# 显式允许的 Origin 白名单（部署在反代后、Host 头与公网域名不一致时设置）：
#   CSRF_ALLOWED_ORIGINS=https://queqiao2.onrender.com,https://example.com
CSRF_ALLOWED_ORIGINS = {
    o.strip().rstrip("/")
    for o in os.environ.get("CSRF_ALLOWED_ORIGINS", "").split(",")
    if o.strip()
}


def _csrf_origin_allowed(header_value: str) -> bool:
    """判断 Origin/Referer 头是否允许（同源或显式白名单）。

    比较 host:port（容忍 http/https scheme 差异，适配反向代理场景），
    避免因代理内部 HTTP 连接与公网 HTTPS 的 scheme 不一致而误杀合法请求。
    """
    from urllib.parse import urlparse

    if not header_value:
        return False
    parsed = urlparse(header_value)
    netloc = parsed.netloc
    if not netloc:
        return False
    # 显式白名单（完整 origin 匹配）
    origin = f"{parsed.scheme}://{netloc}"
    if origin in CSRF_ALLOWED_ORIGINS:
        return True
    # 同源匹配：比较 host:port（容忍 scheme 差异）
    return netloc == request.host


@app.before_request
def _limit_body_size():
    """请求体上限（MAX_CONTENT_LENGTH 仅对解析请求体的端点生效，这里统一拦截）。"""
    limit = app.config.get("MAX_CONTENT_LENGTH")
    if limit and request.content_length and request.content_length > limit:
        return jsonify({"error": "请求体过大"}), 413
    return None


@app.before_request
def _csrf_protect():
    """CSRF 防护：对 POST 请求校验 Origin/Referer 同源。

    浏览器在跨站 POST 时必定发送 ``Origin`` 头；同源请求发送 ``Origin``
    或 ``Referer``。若二者存在但与本服务不同源 → 拒绝（403）。
    两者均缺失（如 curl 等非浏览器客户端）→ 放行（非 CSRF 攻击向量）。
    """
    if request.method != "POST":
        return None
    origin = request.headers.get("Origin")
    if origin:
        if not _csrf_origin_allowed(origin):
            return jsonify({"error": "跨站请求被拒绝 (CSRF)"}), 403
        return None
    referer = request.headers.get("Referer")
    if referer and not _csrf_origin_allowed(referer):
        return jsonify({"error": "跨站请求被拒绝 (CSRF)"}), 403
    # Origin 与 Referer 均缺失：非浏览器客户端，放行
    return None

# 全局 SPICE 半径（首次需要天历时惰性加载）
R_EQ: float = u.MOON_R_EQ_M
R_POL: float = u.MOON_R_EQ_M
_SPICE_LOCK = threading.Lock()
_SPICE_READY = False


# ---------------------------------------------------------------------------
# SPICE 初始化（幂等、线程安全、惰性）
# ---------------------------------------------------------------------------
def init_spice() -> None:
    """加载 SPICE 内核并读取月球半径（幂等，可重复调用）。

    生产 WSGI 服务器（gunicorn/waitress）不会执行 ``__main__``，
    因此内核改由 ``_ensure_spice`` 在首个天历请求时惰性加载。
    """
    global R_EQ, R_POL, _SPICE_READY
    if _SPICE_READY:
        return
    with _SPICE_LOCK:
        if _SPICE_READY:
            return
        u.load_spice_kernels()
        R_EQ, R_POL = u.moon_radii_m()
        _SPICE_READY = True


def _ensure_spice() -> None:
    """惰性 SPICE 初始化：所有需要天历/月固系坐标的端点调用。"""
    init_spice()


# ---------------------------------------------------------------------------
# 访问统计（按 IP + 日期粒度，CSV 持久化）
# ---------------------------------------------------------------------------
# 数据落盘目录：viewer/stats/（可通过环境变量 VIEWER_STATS_DIR 重定向，
# 便于测试与部署时集中管理，如指向 Render 持久盘）。CSV 每次请求追加一行，
# 天然支持按日期/按钮聚合。
STATS_DIR = os.environ.get("VIEWER_STATS_DIR", os.path.join(BASE_DIR, "stats"))
VISIT_CSV = os.path.join(STATS_DIR, "visits.csv")
CLICK_CSV = os.path.join(STATS_DIR, "clicks.csv")

VISIT_HEADER = ["date", "time_utc", "ip", "path"]
CLICK_HEADER = ["date", "time_utc", "ip", "button"]

# CSV 追加写需要互斥，避免并发请求交错（Flask 开发服务器虽单线程，
# 但 gunicorn 等生产部署多 worker 时仍需保护）。
_STATS_LOCK = threading.Lock()

# 可信反向代理网段：仅当请求直连来源（request.remote_addr）属于此集合时，
# 才采信 X-Forwarded-For / X-Real-IP 头，防止客户端直连时伪造 IP。
# 部署时设置环境变量 TRUSTED_PROXIES=10.0.0.0/8,127.0.0.1（Render/Cloudflare 内网）。
# 留空 = 不信任任何代理头，始终使用套接字地址（最安全，适合直连场景）。
TRUSTED_PROXIES = {
    p.strip() for p in os.environ.get("TRUSTED_PROXIES", "").split(",") if p.strip()
}


def _ip_in_trusted_proxies(ip: str) -> bool:
    """判断 IP 是否属于可信代理集合（支持精确 IP 与 CIDR 网段）。"""
    if not ip or ip == "unknown":
        return False
    if ip in TRUSTED_PROXIES:
        return True
    # 支持 CIDR 网段匹配（如 10.0.0.0/8）
    import ipaddress

    try:
        addr = ipaddress.ip_address(ip)
    except ValueError:
        return False
    for cidr in TRUSTED_PROXIES:
        try:
            net = ipaddress.ip_network(cidr, strict=False)
        except ValueError:
            continue
        if addr in net:
            return True
    return False


def _client_ip() -> str:
    """获取客户端 IP（仅可信代理来源才采信转发头）。

    直连场景下 ``X-Forwarded-For`` / ``X-Real-IP`` 可被客户端任意伪造，
    故仅当 ``request.remote_addr`` 属于 ``TRUSTED_PROXIES`` 时才读取转发头；
    否则回退到套接字地址，杜绝 IP 伪造。
    """
    remote = request.remote_addr or "unknown"
    if _ip_in_trusted_proxies(remote):
        forwarded = request.headers.get("X-Forwarded-For")
        if forwarded:
            return forwarded.split(",")[0].strip()
        real_ip = request.headers.get("X-Real-IP")
        if real_ip:
            return real_ip.strip()
    return remote


# ---------------------------------------------------------------------------
# 速率限制：内存滑动窗口限流器，无需额外依赖
# ---------------------------------------------------------------------------
class _RateLimiter:
    """按 key（通常为客户端 IP）的滑动窗口内存限流器。

    线程安全；单进程内存实现，适合单 worker 部署（本项目因 SPICE
    非线程安全而强制单线程/单 worker）。多实例部署需替换为 Redis 后端。
    """

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._buckets: dict[str, list[float]] = {}

    def allow(self, key: str, limit: int, window_sec: float) -> bool:
        """返回 True 表示允许通过，False 表示超限（应返回 429）。"""
        import time

        now = time.monotonic()
        cutoff = now - window_sec
        with self._lock:
            hits = self._buckets.get(key, [])
            # 丢弃窗口外的旧记录
            hits = [t for t in hits if t > cutoff]
            if len(hits) >= limit:
                self._buckets[key] = hits
                return False
            hits.append(now)
            self._buckets[key] = hits
            return True

    def cleanup(self, max_keys: int = 10000) -> None:
        """定期清理过期 key，防止内存无限增长（防限流器自身被 DoS）。"""
        import time

        now = time.monotonic()
        with self._lock:
            if len(self._buckets) < max_keys:
                return
            # 保留最近有活动的 key
            self._buckets = {
                k: v for k, v in self._buckets.items() if v and v[-1] > now - 3600
            }


_rate_limiter = _RateLimiter()


def _rate_limit(limit: int, window_sec: float = 60.0, scope: str | None = None):
    """装饰器：按客户端 IP 限流，超限返回 429 JSON。

    :param limit: 窗口内允许的最大请求数。
    :param window_sec: 滑动窗口长度（秒，默认 60）。
    :param scope: 限流桶作用域（默认取被装饰函数名）。不同端点使用独立
        计数桶，避免高频端点（如 /api/stats 60/min）挤占其他端点的配额。
    """

    def decorator(func):
        from functools import wraps

        bucket = scope or func.__name__

        @wraps(func)
        def wrapper(*args, **kwargs):
            ip = _client_ip()
            if not _rate_limiter.allow(f"{bucket}:{ip}", limit, window_sec):
                resp = jsonify({"error": "请求过于频繁，请稍后再试"})
                resp.status_code = 429
                resp.headers["Retry-After"] = str(int(window_sec))
                return resp
            _rate_limiter.cleanup()
            return func(*args, **kwargs)

        return wrapper

    return decorator


# ---------------------------------------------------------------------------
# 慢端点并发守卫：限制同时处理的慢请求数，超限立即返回 503
# ---------------------------------------------------------------------------
# SPICE 非线程安全 → 慢端点（horizon_mask/terrain_gltf/dem_tile）须串行。
# 单线程部署下本守卫为前瞻性防护；多线程/多 worker 部署时防止慢请求
# 并发拖垮服务。非阻塞获取：若已有慢请求在处理，新请求立即返回 503，
# 避免队列堆积导致全局阻塞。
_SLOW_MAX_CONCURRENT = int(os.environ.get("SLOW_MAX_CONCURRENT", "1"))
_slow_semaphore = threading.BoundedSemaphore(_SLOW_MAX_CONCURRENT)


def _slow_endpoint_guard():
    """装饰器：慢端点并发限制，超限返回 503。

    确保同一时刻最多 ``_SLOW_MAX_CONCURRENT`` 个慢请求在处理（默认 1，
    因 SPICE 非线程安全）。超限请求立即返回 503 + Retry-After，不排队，
    防止慢请求堆积阻塞全部其他请求。
    """

    def decorator(func):
        from functools import wraps

        @wraps(func)
        def wrapper(*args, **kwargs):
            if not _slow_semaphore.acquire(blocking=False):
                resp = jsonify({"error": "服务繁忙，请稍后再试"})
                resp.status_code = 503
                resp.headers["Retry-After"] = "5"
                return resp
            try:
                return func(*args, **kwargs)
            finally:
                _slow_semaphore.release()

        return wrapper

    return decorator


def _stats_utc_now() -> tuple[str, str]:
    """当前 UTC 日期与时间（统计按 UTC 日界分组）。"""
    now = datetime.now(timezone.utc)
    return now.strftime("%Y-%m-%d"), now.strftime("%H:%M:%S")


def _sanitize_csv_field(value: str) -> str:
    """中和 CSV 公式注入。

    以 ``=``、``+``、``-``、``@``、Tab、CR 开头的字段会被 Excel/WPS/
    LibreOffice 解释为公式，可触发 DDE 命令执行、恶意超链接或数据外泄。
    在此类字段前前置单引号（``'``）即可令电子表格将其视为纯文本。

    参见 OWASP CSV Injection 指南。
    """
    s = str(value)
    if s and s[0] in ("=", "+", "-", "@", "\t", "\r"):
        return "'" + s
    return s


def _append_csv(path: str, header: list[str], row: list[str]) -> None:
    """线程安全地追加一行 CSV；文件不存在时先写表头。

    所有字段在写入前经 ``_sanitize_csv_field`` 中和公式注入字符，
    确保无论数据来源（用户参数 / IP 头）均不可触发电子表格公式执行。
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    safe_row = [_sanitize_csv_field(f) for f in row]
    with _STATS_LOCK:
        is_new = not os.path.exists(path)
        with open(path, "a", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            if is_new:
                writer.writerow(header)
            writer.writerow(safe_row)


def _read_csv(path: str) -> list[dict]:
    """读取 CSV 为 dict 列表；文件缺失或格式损坏时返回空列表。"""
    if not os.path.exists(path):
        return []
    try:
        with open(path, "r", newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            return [dict(r) for r in reader if any(r.values())]
    except (OSError, csv.Error):
        return []


def record_visit(ip: str, path: str) -> None:
    """记录一次页面访问（IP + 时间 + 路径）。"""
    date, t = _stats_utc_now()
    _append_csv(VISIT_CSV, VISIT_HEADER, [date, t, ip, path])


def record_click(ip: str, button: str) -> None:
    """记录一次按钮点击（IP + 时间 + 按钮名）。"""
    date, t = _stats_utc_now()
    _append_csv(CLICK_CSV, CLICK_HEADER, [date, t, ip, button])


def stats_summary() -> dict:
    """按日期聚合统计：访问量（总数 + 唯一 IP）+ 各按钮点击量。

    :return::
        {
          "days": [{"date", "visits", "unique_ips", "clicks": {btn: n}}],
          "total_visits": int, "total_unique_ips": int, "total_clicks": int,
        }
    """
    visits = _read_csv(VISIT_CSV)
    clicks = _read_csv(CLICK_CSV)

    daily_visits: dict[str, int] = {}
    daily_unique_ips: dict[str, set[str]] = {}
    for r in visits:
        d = r.get("date", "") or "unknown"
        daily_visits[d] = daily_visits.get(d, 0) + 1
        daily_unique_ips.setdefault(d, set()).add(r.get("ip", ""))

    daily_clicks: dict[str, dict[str, int]] = {}
    for r in clicks:
        d = r.get("date", "") or "unknown"
        btn = (r.get("button", "") or "unknown").strip() or "unknown"
        day = daily_clicks.setdefault(d, {})
        day[btn] = day.get(btn, 0) + 1

    dates = sorted(set(list(daily_visits) + list(daily_unique_ips) + list(daily_clicks)))
    rows = []
    for d in dates:
        rows.append(
            {
                "date": d,
                "visits": daily_visits.get(d, 0),
                "unique_ips": len(daily_unique_ips.get(d, set())),
                "clicks": daily_clicks.get(d, {}),
            }
        )

    return {
        "days": rows,
        "total_visits": len(visits),
        "total_unique_ips": len({r.get("ip") for r in visits if r.get("ip")}),
        "total_clicks": len(clicks),
    }


# ---------------------------------------------------------------------------
# DEM 可用性（缺失时优雅降级，健康检查暴露状态）
# ---------------------------------------------------------------------------
def dem_available() -> bool:
    return bool(DEM_PATH) and os.path.isfile(DEM_PATH)


# ---------------------------------------------------------------------------
# 工具函数
# ---------------------------------------------------------------------------
def _lunar_cmap():
    """月球灰度配色（LDEM 月面纹理）：低洼深灰 → 高原中灰（纯灰阶，避免黑白反差突兀）。"""
    return mcolors.LinearSegmentedColormap.from_list(
        "lunar_gray",
        [
            (0.00, (0.20, 0.20, 0.20)),   # 最深阴影（坑底）深灰
            (0.50, (0.38, 0.38, 0.38)),   # 中灰
            (1.00, (0.55, 0.55, 0.55)),   # 高原/坑沿 中浅灰
        ],
    )


def _mjd_to_et(t_mjd: float) -> float:
    """MJD → SPICE ET（秒）。

    所有 SPICE 调用都经过本函数（地球/太阳/中继位置），因此在这里惰性
    初始化内核，保证 WSGI 服务器（不执行 ``__main__``）首次请求即可用。
    """
    _ensure_spice()
    import spiceypy

    return spiceypy.unitim(float(t_mjd) + 2400000.5, "JDTDB", "ET")


def _parse_float_arg(name: str, default: float) -> float:
    """安全解析查询参数为有限浮点数；非法/非有限输入抛 BadRequest。

    拒绝非数字字符串（``abc``）、溢出值（``1e999`` → inf）以及显式
    ``inf``/``NaN``，避免 ``math domain error`` 500 与 NaN 静默传播。
    """
    raw = request.args.get(name)
    if raw is None:
        return float(default)
    try:
        val = float(raw)
    except (TypeError, ValueError):
        raise BadRequest(f"参数 {name} 必须为数字")
    if not math.isfinite(val):
        raise BadRequest(f"参数 {name} 必须为有限数（不接受 inf/NaN）")
    return val


def _parse_location() -> tuple[float, float, float]:
    """解析经纬度查询参数，返回 (lat, lon, h)。非法输入抛 BadRequest。"""
    lat = _parse_float_arg("lat", SHACKLETON_LAT_DEG)
    lon = _parse_float_arg("lon", SHACKLETON_LON_DEG)
    h = _parse_float_arg("h", 0.0)
    if not (-90.0 <= lat <= 90.0):
        raise BadRequest("lat 必须在 -90..90 范围内")
    if not (-180.0 <= lon <= 180.0):
        raise BadRequest("lon 必须在 -180..180 范围内")
    return lat, lon, h


def _parse_time() -> tuple[str, float, str]:
    """解析日期 + 小时查询参数，返回 (t_utc, t_mjd, date)。非法输入抛 BadRequest。"""
    date = request.args.get("date", DEFAULT_DATE)
    try:
        datetime.strptime(date, "%Y-%m-%d")
    except (TypeError, ValueError):
        raise BadRequest("date 必须为 YYYY-MM-DD 格式的有效日期")
    hour = _parse_float_arg("hour", 12.0)
    hour = max(0.0, min(23.999, hour))
    hh = int(hour)
    mm = int(round((hour - hh) * 60.0))
    if mm >= 60:
        mm = 0
        hh += 1
    t_utc = f"{date}T{hh:02d}:{mm:02d}:00"
    t_mjd = u.utc_iso_to_mjd(t_utc)
    return t_utc, t_mjd, date


def _earth_pos_iau_moon_m(t_mjd: float) -> np.ndarray:
    """地球在 IAU_MOON 月固系中的位置（米）。"""
    import spiceypy

    et = _mjd_to_et(t_mjd)
    state, _lt = spiceypy.spkezr("EARTH", et, "IAU_MOON", "LT+S", "MOON")
    return np.asarray(state[:3], dtype=float) * 1000.0


def _sun_pos_iau_moon_m(t_mjd: float) -> np.ndarray:
    """太阳在 IAU_MOON 月固系中的位置（米）。

    SPICE 优先；内核缺失/调用失败时回退到低精度解析模型（与
    ``illumination._analytical_sun_vector`` 同口径：子太阳点纬度 0、
    经度按朔望月相位推进），保证任意环境都能渲染光照瓦片。
    """
    try:
        import spiceypy

        et = _mjd_to_et(t_mjd)
        state, _lt = spiceypy.spkezr("SUN", et, "IAU_MOON", "LT+S", "MOON")
        return np.asarray(state[:3], dtype=float) * 1000.0
    except Exception:
        synodic_days = 29.53059
        phase = 2.0 * math.pi * (t_mjd - 62502.0) / synodic_days
        sub_solar_lat = 0.0
        sub_solar_lon = phase
        sun_r = 1.495978707e11
        return np.array(
            [
                sun_r * math.cos(sub_solar_lat) * math.cos(sub_solar_lon),
                sun_r * math.cos(sub_solar_lat) * math.sin(sub_solar_lon),
                sun_r * math.sin(sub_solar_lat),
            ],
            dtype=float,
        )


def _sun_elev_az(lat_deg: float, lon_deg: float, h_m: float, t_mjd: float) -> tuple[float, float]:
    """太阳高度角 + 方位角（度；方位角 0°=北、90°=东，与地平线图约定一致）。"""
    sun_pos_m = _sun_pos_iau_moon_m(t_mjd)
    return u.enu_elevation_azimuth(sun_pos_m, lon_deg, lat_deg, h_m, R_EQ, R_POL)


def _relay_pos_iau_moon_m(et_list: np.ndarray, pos_j2000_km: np.ndarray) -> np.ndarray:
    """J2000 位置（km）→ IAU_MOON 位置（米）。"""
    import spiceypy

    out = np.empty_like(pos_j2000_km)
    for i, et in enumerate(et_list):
        rot = spiceypy.pxform("J2000", "IAU_MOON", et)
        out[i] = spiceypy.mxv(rot, pos_j2000_km[i]) * 1000.0
    return out


def _relay_mission_state() -> tuple[float, np.ndarray, float]:
    """在**固定任务历元**构造中继轨道初值（与 relay_elevation.py 口径一致）。

    关键：用 ``MISSION_EPOCH_UTC``（默认 2026-08-01 00:00，假设平近点角
    M=215°）构造轨道状态，然后由月心二体连续传播到任意查询时刻。

    :return: (et_mission, state0 [km, km/s], mu_moon)
    """
    et_mission = float(_mjd_to_et(u.utc_iso_to_mjd(MISSION_EPOCH_UTC)))
    state0, mu = build_orbital_state(et_mission)
    return et_mission, state0, mu


def _relay_daily_positions(date: str) -> tuple[np.ndarray, np.ndarray]:
    """预计算一天 288 点的中继卫星 IAU_MOON 位置（米）。

    轨道状态在**固定任务历元**构造后连续传播到当日各采样点（与
    ``relay_elevation.py`` 一致），避免每日重排轨道相位。

    :return: (et_list, pos_iau_m) — pos_iau_m 形状 (n, 3)，单位米
    """
    t_mjd_start = u.utc_iso_to_mjd(f"{date}T00:00:00")
    t_mjd_arr = t_mjd_start + np.arange(DAILY_N_POINTS) * (
        DAILY_INTERVAL_MIN / 60.0 / 24.0
    )
    et_list = np.array([_mjd_to_et(float(t)) for t in t_mjd_arr], dtype=float)
    et_mission, state0, mu = _relay_mission_state()
    pos_j2000 = propagate_satellite(state0, mu, et_mission, et_list)
    pos_iau = _relay_pos_iau_moon_m(et_list, pos_j2000)
    return et_list, pos_iau


def _grid_center_elevation(data: np.ndarray) -> float:
    """网格中心 (0, 0) 的双线性插值高程（与 meshgrid 局部坐标一致）。

    偶数网格（如 288×288）的整数中心索引 (GRID//2, GRID//2) 对应局部坐标
    (+52 m, +52 m) 而非原点；在陡坡处该偏移会造成数十米高程误差，导致
    GLTF 相机眼高低于地表（3D 视角被灰色地形填满）。此函数对连续网格
    中心索引 (N-1)/2 进行双线性插值，返回观测点正下方地表的精确高程。
    """
    n = data.shape[0]
    c = (n - 1) / 2.0
    i = int(np.floor(c))
    f = c - i
    return float(
        data[i, i] * (1 - f) * (1 - f)
        + data[i + 1, i] * f * (1 - f)
        + data[i, i + 1] * (1 - f) * f
        + data[i + 1, i + 1] * f * f
    )


def build_terrain_gltf(lat_deg: float, lon_deg: float) -> tuple[bytes, dict]:
    """构建观测点周围 3D 地形网格（月球灰度纹理）并导出 GLTF bytes。

    坐标系：以观测点为原点，X=东、Y=北、Z=高程（VTK Z-up）。PyVista 的
    ``export_gltf`` 自动做 Z-up → glTF Y-up 转换（X→Z, Y→X, Z→Y），
    由浏览器 three.js 加载渲染（第一人称天空视角）。

    :return: (gltf_bytes, meta) — meta 含观测点高程等，用于前端相机定位
    """
    import pyvista as pv

    x0, y0 = u.lonlat_to_polar_stereo_xy(lon_deg, lat_deg)
    half = TERRAIN_HALF_M

    # 读取 DEM 窗口；越界/无数据时回退为"平坦参考地形"（高程 0），
    # 保证任意观测点都能渲染（其余要素——星空/地球/鹊桥/地平线——不受影响）。
    fallback = False
    try:
        with rasterio.open(DEM_PATH) as ds:
            win = from_bounds(
                x0 - half, y0 - half, x0 + half, y0 + half, ds.transform
            )
            data = ds.read(
                1, window=win, out_shape=(TERRAIN_GRID, TERRAIN_GRID),
                resampling=Resampling.bilinear, masked=True,
            )
            data = np.where(data.mask, np.nan, data.filled(np.nan))
    except Exception as exc:  # RasterioIOError / WindowError 等
        print(f"[terrain_gltf] DEM 读取失败，回退平坦地形: {type(exc).__name__}: {exc}")
        fallback = True

    if not fallback:
        valid = np.isfinite(data)
        if not valid.any():
            print("[terrain_gltf] 观测点周围无有效 DEM 数据，回退平坦地形")
            fallback = True
        else:
            # 无效像素用有效中值填充，避免网格空洞
            data = np.where(valid, data, float(np.nanmedian(data[valid])))

    if fallback:
        data = np.zeros((TERRAIN_GRID, TERRAIN_GRID), dtype=np.float32)

    # 回退后重新计算有效掩码（zeros 全有效），并保证掩码形状正确
    valid = np.isfinite(data)

    # DEM 行序：第 0 行 = 最北（y 最大）→ 翻转为 y 递增，与网格一致
    data = data[::-1, :]

    # 局部坐标网格（单位米）
    xs = np.linspace(-half, half, TERRAIN_GRID)
    ys = np.linspace(-half, half, TERRAIN_GRID)
    xx, yy = np.meshgrid(xs, ys)

    grid = pv.StructuredGrid(xx, yy, data.astype(np.float32))

    # 月球灰度顶点色（按高程归一化，与 2D 瓦片同色表）
    zmin, zmax = np.percentile(data[valid], [2, 98])
    norm = np.clip((data - zmin) / (zmax - zmin + 1e-9), 0.0, 1.0)
    rgba = _lunar_cmap()(norm)
    rgb = (rgba[:, :, :3] * 255.0).astype(np.uint8)
    grid["RGB"] = rgb.reshape(-1, 3)

    surf = grid.extract_surface(algorithm="dataset_surface")

    plotter = pv.Plotter(off_screen=True)
    try:
        plotter.add_mesh(surf, rgb=True, smooth_shading=True)
        # 使用 NamedTemporaryFile 原子创建不可预测的临时文件名，
        # 防止符号链接劫持覆盖任意文件。
        with tempfile.NamedTemporaryFile(suffix=".gltf", delete=False) as tmp:
            tmp_path = tmp.name
        try:
            plotter.export_gltf(tmp_path)
            with open(tmp_path, "rb") as f:
                gltf_bytes = f.read()
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
    finally:
        plotter.close()

    # 观测点高程（网格中心），用于前端相机眼高定位。
    # 对网格中心 (0, 0) 双线性插值取精确高程，保证相机始终位于地表上方。
    obs_z = _grid_center_elevation(data)
    meta = {
        "terrain_half_m": float(half),
        "grid_size": TERRAIN_GRID,
        "obs_z_m": obs_z,
        "lat_deg": float(lat_deg),
        "lon_deg": float(lon_deg),
    }
    return gltf_bytes, meta


# ---------------------------------------------------------------------------
# 缓存
# ---------------------------------------------------------------------------
@lru_cache(maxsize=128)
def cached_horizon_mask(lat_deg: float, lon_deg: float) -> dict:
    """缓存地平线掩码计算结果（按 0.01° 精度的经纬度键）。"""
    x0, y0 = u.lonlat_to_polar_stereo_xy(lon_deg, lat_deg)
    return compute_horizon_mask(
        dem_path=DEM_PATH,
        x0=x0, y0=y0,
        max_radius_m=50000.0,
        n_az=360,
        n_samples=300,
    )


def _render_illuminated_tile(
    data: np.ndarray, sun_elev_deg: float, sun_az_deg: float, pixel_size_m: float
) -> tuple[np.ndarray, float, float]:
    """把 DEM 高程数组渲染为 RGBA 光照图（余弦山体阴影 + 月球灰阶反照率）。

    - 表面法向由 DEM 梯度计算（行序：第 0 行 = 最北，行增 = 南向）；
    - 亮度 = 环境光下限 + 直射余弦项 max(0, n·s)。太阳略低于几何地平线
      时，朝向太阳的低倾角坡面仍可受光，呈现月南极"坑缘低角度日照、
      坑底深阴影"的月球质感，无需额外的地平线判断；
    - 无数据像素保持透明（RGB 清零，仅 alpha 区分）。
    """
    valid = np.isfinite(data)
    fill = float(np.nanmedian(data[valid])) if valid.any() else 0.0
    z = np.where(valid, data, fill).astype(np.float64)

    # 梯度：axis=1 = 列（东向），axis=0 = 行（南向）→ 北向梯度取负
    dz_ds, dz_de = np.gradient(z, pixel_size_m, pixel_size_m, edge_order=2)
    dz_dn = -dz_ds

    # 表面向上法向（ENU 局部坐标）
    n = np.empty(z.shape + (3,), dtype=np.float64)
    n[..., 0] = -dz_de
    n[..., 1] = -dz_dn
    n[..., 2] = 1.0
    norm = np.linalg.norm(n, axis=-1, keepdims=True)
    n /= np.maximum(norm, 1e-12)

    # 太阳方向（ENU：E=cos(elev)·sin(az), N=cos(elev)·cos(az), U=sin(elev)）
    elev_rad = math.radians(sun_elev_deg)
    az_rad = math.radians(sun_az_deg)
    s = np.array(
        [
            math.cos(elev_rad) * math.sin(az_rad),
            math.cos(elev_rad) * math.cos(az_rad),
            math.sin(elev_rad),
        ]
    )

    cos_i = np.clip(n @ s, 0.0, 1.0)

    # 基础反照率：月球灰阶（高程 2%–98% 归一化）
    if valid.any():
        vmin, vmax = np.percentile(data[valid], [2, 98])
    else:
        vmin, vmax = -3000.0, 1000.0
    base = np.clip((data - vmin) / (vmax - vmin + 1e-9), 0.0, 1.0)
    base_rgb = np.asarray(_lunar_cmap()(base))[..., :3]

    # 光照合成：环境光下限（避免纯黑，代表地照/天光）+ 直射余弦项
    ambient = 0.08
    shade = ambient + (1.0 - ambient) * cos_i
    rgb = np.where(valid[..., None], base_rgb * shade[..., None], 0.0)
    alpha = np.where(valid, 1.0, 0.0)
    img = (np.dstack([rgb, alpha]) * 255.0).astype(np.uint8)
    return img, float(vmin), float(vmax)


@lru_cache(maxsize=128)
def cached_dem_tile_png(
    center_lat: float,
    center_lon: float,
    width_m: float,
    height_m: float,
    t_mjd: float,
) -> tuple[bytes, dict]:
    """渲染 20 km × 20 km DEM 瓦片 PNG：月球灰阶反照率 × 当前太阳方向山体阴影。

    太阳高度/方位按瓦片中心经纬度 + 请求时刻（MJD）计算（SPICE 优先、
    解析回退）；时间进入缓存键 → 拖动时间滑块可实时改变光照方向。
    """
    x0, y0 = u.lonlat_to_polar_stereo_xy(center_lon, center_lat)
    with rasterio.open(DEM_PATH) as ds:
        win = from_bounds(
            x0 - width_m / 2.0, y0 - height_m / 2.0,
            x0 + width_m / 2.0, y0 + height_m / 2.0,
            ds.transform,
        )
        data = ds.read(
            1, window=win, out_shape=(TILE_PX, TILE_PX),
            resampling=Resampling.bilinear, masked=True,
        )
        data = np.where(data.mask, np.nan, data.filled(np.nan))

    sun_elev_deg, sun_az_deg = _sun_elev_az(center_lat, center_lon, 0.0, t_mjd)
    img, vmin, vmax = _render_illuminated_tile(
        data, sun_elev_deg, sun_az_deg, width_m / TILE_PX
    )
    from PIL import Image

    pil = Image.fromarray(img, mode="RGBA")
    buf = io.BytesIO()
    pil.save(buf, format="PNG")
    buf.seek(0)

    pixel_size = width_m / TILE_PX
    meta = {
        "center_x_stereo_m": float(x0),
        "center_y_stereo_m": float(y0),
        "center_lat_deg": center_lat,
        "center_lon_deg": center_lon,
        "width_m": float(width_m),
        "height_m": float(height_m),
        "pixel_size_m": pixel_size,
        "img_width_px": TILE_PX,
        "img_height_px": TILE_PX,
        "vmin_m": vmin,
        "vmax_m": vmax,
        "sun_elev_deg": float(sun_elev_deg),
        "sun_az_deg": float(sun_az_deg),
    }
    return buf.getvalue(), meta


@lru_cache(maxsize=128)
def cached_terrain_gltf(lat_deg: float, lon_deg: float) -> tuple[bytes, dict]:
    """缓存观测点周围 3D 地形 GLTF（按 0.01° 精度经纬度键）。"""
    return build_terrain_gltf(lat_deg, lon_deg)


# ---------------------------------------------------------------------------
# Yang 论文轨道资产（预生成 JSON，只读；不暴露文件系统路径）
# ---------------------------------------------------------------------------
def _json_safe(value):
    """递归把非有限浮点（NaN/±Inf）转为 None，保证输出为严格 JSON。

    ``OrbitRecord.to_jsonable()`` 在求解失败（``error`` 记录）时会写入
    ``NaN``；Python ``json`` 默认允许输出裸 ``NaN``（非标准 JSON），
    严格解析器（如浏览器 ``JSON.parse``）会报错，故在服务端统一净化。
    """
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {k: _json_safe(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_json_safe(v) for v in value]
    return value


@lru_cache(maxsize=1)
def _yang_orbits_data() -> dict | None:
    """读取并缓存轨道 JSON 资产；缺失/损坏时返回 None（端点返回 503）。"""
    if not os.path.isfile(YANG_ORBITS_JSON):
        return None
    try:
        with open(YANG_ORBITS_JSON, "r", encoding="utf-8") as f:
            payload = json.load(f)
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(payload, dict) or not isinstance(payload.get("orbits"), list):
        return None
    return _json_safe(payload)


@lru_cache(maxsize=1)
def _yang_orbits_json_bytes() -> bytes | None:
    """缓存序列化后的轨道 JSON（严格模式，拒绝 NaN/Inf）。"""
    data = _yang_orbits_data()
    if data is None:
        return None
    return json.dumps(data, ensure_ascii=False, allow_nan=False).encode("utf-8")


def _parse_t_arrive() -> str:
    """解析并严格校验 ``t_arrive`` 查询参数（UTC ISO-8601）。

    仅接受 ``YYYY-MM-DD`` 或 ``YYYY-MM-DDTHH:MM:SSZ``，且年份须落在
    SPICE 内核（de442s.bsp）覆盖的 1900–2050；非法输入抛 BadRequest（400），
    不把用户字符串直接传给 SPICE。
    """
    raw = request.args.get("t_arrive")
    if raw is None:
        return YANG_ORBITS_DEFAULT_T_ARRIVE
    s = raw.strip()
    if not s or len(s) > 32:
        raise BadRequest("t_arrive 必须为 UTC ISO-8601 时间字符串")
    try:
        if "T" in s:
            dt = datetime.strptime(s, "%Y-%m-%dT%H:%M:%SZ")
        else:
            dt = datetime.strptime(s, "%Y-%m-%d")
    except (TypeError, ValueError):
        raise BadRequest("t_arrive 必须为 YYYY-MM-DD 或 YYYY-MM-DDTHH:MM:SSZ（UTC）")
    if not (1900 <= dt.year <= 2050):
        raise BadRequest("t_arrive 必须在 1900–2050 年（SPICE 内核覆盖范围）")
    return dt.strftime("%Y-%m-%dT%H:%M:%SZ")


def _yang_orbits_meta(data: dict) -> dict:
    """轨道资产元信息（供 SSE meta 事件与前端状态栏使用）。"""
    return {
        "t_arrive_utc": data.get("t_arrive_utc"),
        "h_p_km": data.get("h_p_km"),
        "H_m_km": data.get("H_m_km"),
        "n_samples_per_orbit": data.get("n_samples_per_orbit"),
        "moon_model": data.get("moon_model"),
        "n_orbits_expected": len(data["orbits"]),
    }


# ---------------------------------------------------------------------------
# 融合反向求解（fusion FR-4）作业调度
# ---------------------------------------------------------------------------
@lru_cache(maxsize=1)
def _fusion_supported_T() -> tuple[float, ...]:
    """从预生成家族表目录解析支持的转移时间 T（小时）。

    家族表命名由 ``fusion._table_path`` 决定（i_m=90000、Hm200、meq、sun），
    对应极轨捕获（i_m=90°）、月球真极轴参考、近月点 200 km、含太阳摄动。
    """
    values: set[float] = set()
    if os.path.isdir(FUSION_TABLE_DIR):
        for name in os.listdir(FUSION_TABLE_DIR):
            m = re.fullmatch(r"r1s_scan_im90000_T(\d{4})_Hm200_meq_sun\.json", name)
            if m:
                values.add(int(m.group(1)) / 10.0)
    return tuple(sorted(values))


def _fusion_snap_T(t_h: float) -> float | None:
    """把请求的 T 吸附到最近的预生成表格；偏差 > 0.35 h 视为不支持。"""
    supported = _fusion_supported_T()
    if not supported:
        return None
    nearest = min(supported, key=lambda s: abs(s - t_h))
    return nearest if abs(nearest - t_h) <= 0.35 else None


def _parse_date_field(value, name: str) -> datetime:
    """严格解析作业日期字段（YYYY-MM-DD，1900–2050）。"""
    if not isinstance(value, str):
        raise BadRequest(f"{name} 必须为 YYYY-MM-DD 格式")
    try:
        dt = datetime.strptime(value.strip(), "%Y-%m-%d")
    except (TypeError, ValueError):
        raise BadRequest(f"{name} 必须为 YYYY-MM-DD 格式的有效日期")
    if not (1900 <= dt.year <= 2050):
        raise BadRequest(f"{name} 必须在 1900–2050 年（SPICE 内核覆盖范围）")
    return dt


def _fusion_job_dir(job_id: str) -> str:
    """作业目录（job_id 已通过 16 位十六进制白名单校验，无路径注入面）。"""
    return os.path.join(FUSION_JOB_DIR, job_id)


def _fusion_job_meta(job_id: str) -> dict | None:
    path = os.path.join(_fusion_job_dir(job_id), "meta.json")
    if not os.path.isfile(path):
        return None
    try:
        with open(path, "r", encoding="utf-8") as f:
            meta = json.load(f)
    except (OSError, json.JSONDecodeError):
        return None
    return meta if isinstance(meta, dict) else None


def _write_fusion_marker(job_id: str, status: str, error: str | None = None) -> None:
    try:
        with open(
            os.path.join(_fusion_job_dir(job_id), "state.json"), "w", encoding="utf-8"
        ) as f:
            json.dump({"status": status, "error": error}, f, ensure_ascii=False)
    except OSError:
        pass


def _read_fusion_marker(job_id: str) -> dict | None:
    path = os.path.join(_fusion_job_dir(job_id), "state.json")
    try:
        with open(path, "r", encoding="utf-8") as f:
            marker = json.load(f)
    except (OSError, json.JSONDecodeError):
        return None
    return marker if isinstance(marker, dict) else None


def _terminate_fusion_proc(proc: subprocess.Popen) -> None:
    """终止作业子进程组（start_new_session=True，可整组回收）。"""
    try:
        os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
    except (ProcessLookupError, PermissionError, OSError):
        try:
            proc.terminate()
        except OSError:
            pass
    try:
        proc.wait(timeout=10)
    except (subprocess.TimeoutExpired, OSError):
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except (ProcessLookupError, PermissionError, OSError):
            pass


def _fusion_cleanup_old_jobs(max_age_s: float = 48 * 3600.0) -> None:
    """清理超过 48 h 的旧作业目录（防止磁盘无限增长）。"""
    if not os.path.isdir(FUSION_JOB_DIR):
        return
    now = time.time()
    try:
        entries = os.listdir(FUSION_JOB_DIR)
    except OSError:
        return
    for name in entries:
        path = os.path.join(FUSION_JOB_DIR, name)
        try:
            if os.path.isdir(path) and now - os.path.getmtime(path) > max_age_s:
                shutil.rmtree(path, ignore_errors=True)
        except OSError:
            continue


_fusion_lock = threading.Lock()
_fusion_active: dict | None = None


@app.route("/coast_tool/api/fusion/config")
@_rate_limit(120, 60.0, scope="fusion_config")
def fusion_config():
    """融合反解可用配置（支持的 T 列表、站点、规模与耗时估计）。"""
    return jsonify(
        {
            "supported_T_h": list(_fusion_supported_T()),
            "sites": list(FUSION_SITES),
            "i_m_deg": FUSION_I_M_DEG,
            "polar_ref": "moon_equator",
            "max_span_days": FUSION_MAX_SPAN_DAYS,
            "est_sec_per_day": FUSION_EST_SEC_PER_DAY,
            "ready": bool(_fusion_supported_T()) and os.path.isfile(FUSION_SCRIPT),
        }
    )


@app.route("/coast_tool/api/fusion/jobs", methods=["POST"])
@_rate_limit(10, 600.0, scope="fusion_start")
def fusion_start_job():
    """启动一次融合反向求解作业（长耗时，须 confirm=true）。

    请求体（JSON）::

        {"start": "2026-08-24", "end": "2026-08-25", "T_h": 112.8,
         "site": "wenchang", "coast_max_s": 1200, "confirm": true}

    仅支持预生成家族表覆盖的 T（自动吸附 ±0.35 h）；同一时刻只允许一个
    作业（CPU 密集），重复提交返回 409。
    """
    global _fusion_active
    body = request.get_json(silent=True)
    if not isinstance(body, dict):
        raise BadRequest("请求体必须为 JSON 对象")
    if str(body.get("confirm", "")).lower() not in ("1", "true", "yes", "on"):
        raise BadRequest("长耗时计算须显式确认（confirm=true）")

    start_dt = _parse_date_field(body.get("start"), "start")
    end_dt = _parse_date_field(body.get("end"), "end")
    span_days = (end_dt - start_dt).days
    if span_days < 1:
        raise BadRequest("end 必须晚于 start（右开区间，至少 1 天）")
    if span_days > FUSION_MAX_SPAN_DAYS:
        raise BadRequest(f"扫描范围最长 {FUSION_MAX_SPAN_DAYS} 天")

    try:
        t_requested = float(body.get("T_h"))
    except (TypeError, ValueError):
        raise BadRequest("T_h 必须为数字（小时）")
    if not math.isfinite(t_requested):
        raise BadRequest("T_h 必须为有限数")
    t_h = _fusion_snap_T(t_requested)
    if t_h is None:
        supported = ", ".join(f"{s:g}" for s in _fusion_supported_T()) or "无"
        raise BadRequest(f"当前仅支持 T ∈ {{{supported}}} h（自动吸附 ±0.35 h）")

    site = str(body.get("site", "wenchang")).strip().lower()
    if site not in FUSION_SITES:
        raise BadRequest(f"site 必须为 {', '.join(FUSION_SITES)} 之一")

    try:
        coast_max_s = float(body.get("coast_max_s", 1200.0))
    except (TypeError, ValueError):
        raise BadRequest("coast_max_s 必须为数字")
    if not math.isfinite(coast_max_s):
        raise BadRequest("coast_max_s 必须为有限数")
    coast_max_s = min(max(coast_max_s, 0.0), 3600.0)

    if not os.path.isfile(FUSION_SCRIPT):
        return jsonify({"error": "融合求解脚本缺失"}), 503
    if not _fusion_supported_T():
        return jsonify({"error": "融合家族表资产缺失"}), 503

    with _fusion_lock:
        if _fusion_active is not None and _fusion_active["proc"].poll() is None:
            return jsonify({"error": "已有融合反解作业在运行，请稍后再试"}), 409
        _fusion_cleanup_old_jobs()
        job_id = uuid.uuid4().hex[:16]
        job_dir = _fusion_job_dir(job_id)
        os.makedirs(job_dir, exist_ok=False)
        config = {
            "start": start_dt.strftime("%Y-%m-%d"),
            "end": end_dt.strftime("%Y-%m-%d"),
            "T_h": t_h,
            "T_requested_h": t_requested,
            "site": site,
            "coast_max_s": coast_max_s,
            "i_m_deg": FUSION_I_M_DEG,
            "polar_ref": "moon_equator",
            "span_days": span_days,
        }
        with open(
            os.path.join(job_dir, "meta.json"), "w", encoding="utf-8"
        ) as f:
            json.dump(
                {
                    "id": job_id,
                    "created_utc": datetime.now(timezone.utc).strftime(
                        "%Y-%m-%dT%H:%M:%SZ"
                    ),
                    "config": config,
                },
                f,
                ensure_ascii=False,
            )
        # 最小化子进程环境（不继承服务端可能存在的密钥类变量）
        env = {
            "PATH": os.environ.get("PATH", "/usr/bin:/bin"),
            "HOME": os.environ.get("HOME", "/tmp"),
            "PYTHONUNBUFFERED": "1",
            "FUSION_TABLE_DIR": FUSION_TABLE_DIR,
        }
        if os.environ.get("SPICE_KERNEL_DIR"):
            env["SPICE_KERNEL_DIR"] = os.environ["SPICE_KERNEL_DIR"]
        stdout_path = os.path.join(job_dir, "stdout.log")
        with open(stdout_path, "ab") as out:
            proc = subprocess.Popen(  # noqa: S603 - 参数均来自白名单校验
                [
                    sys.executable,
                    FUSION_SCRIPT,
                    "--start",
                    config["start"],
                    "--end",
                    config["end"],
                    "--T",
                    f"{t_h:g}",
                    "--site",
                    site,
                    "--coast-max",
                    f"{coast_max_s:g}",
                    "--out",
                    job_dir,
                ],
                cwd=_REPO_ROOT,
                stdout=out,
                stderr=subprocess.STDOUT,
                env=env,
                start_new_session=True,
            )
        _fusion_active = {
            "id": job_id,
            "proc": proc,
            "started": time.monotonic(),
        }

    return (
        jsonify(
            {
                "job_id": job_id,
                "status": "running",
                "config": config,
                "estimated_s": int(span_days * FUSION_EST_SEC_PER_DAY),
            }
        ),
        202,
    )


@app.route("/coast_tool/api/fusion/jobs/<job_id>")
@_rate_limit(240, 60.0, scope="fusion_status")
def fusion_job_status(job_id: str):
    """查询作业状态；完成后返回 ``result``（窗口列表 + 拒绝统计）。"""
    global _fusion_active
    if not re.fullmatch(r"[0-9a-f]{16}", job_id or ""):
        raise BadRequest("job_id 非法")
    meta = _fusion_job_meta(job_id)
    if meta is None:
        return jsonify({"error": "作业不存在或已过期"}), 404

    with _fusion_lock:
        active = (
            _fusion_active
            if (_fusion_active is not None and _fusion_active["id"] == job_id)
            else None
        )
        proc = active["proc"] if active else None
        if proc is not None and proc.poll() is not None:
            _fusion_active = None
            active, proc = None, None
        elapsed = (time.monotonic() - active["started"]) if active else None
        if (
            proc is not None
            and elapsed is not None
            and elapsed > FUSION_JOB_TIMEOUT_S
        ):
            _terminate_fusion_proc(proc)
            _fusion_active = None
            _write_fusion_marker(job_id, "failed", "作业超时被终止")
            proc, elapsed = None, None

    payload: dict = {
        "job_id": job_id,
        "status": "unknown",
        "config": meta.get("config"),
    }
    if proc is not None:
        payload.update(status="running", elapsed_s=round(elapsed or 0.0, 1))
        return jsonify(payload)

    result_path = os.path.join(_fusion_job_dir(job_id), "result.json")
    marker = _read_fusion_marker(job_id)
    if os.path.isfile(result_path):
        try:
            with open(result_path, "r", encoding="utf-8") as f:
                result = json.load(f)
            payload.update(status="done", result=result)
        except (OSError, json.JSONDecodeError):
            payload.update(status="failed", error="结果文件损坏")
    elif marker and marker.get("status") == "cancelled":
        payload.update(status="cancelled")
    elif marker and marker.get("status") == "failed":
        payload.update(status="failed", error=marker.get("error") or "作业失败")
    else:
        payload.update(
            status="unknown", error="作业进程句柄已丢失（服务可能重启过），可稍后重试查询"
        )
    return jsonify(payload)


@app.route("/coast_tool/api/fusion/jobs/<job_id>", methods=["DELETE"])
@_rate_limit(10, 60.0, scope="fusion_cancel")
def fusion_cancel_job(job_id: str):
    """取消运行中的作业（终止子进程组）。"""
    global _fusion_active
    if not re.fullmatch(r"[0-9a-f]{16}", job_id or ""):
        raise BadRequest("job_id 非法")
    if _fusion_job_meta(job_id) is None:
        return jsonify({"error": "作业不存在或已过期"}), 404
    with _fusion_lock:
        if _fusion_active is not None and _fusion_active["id"] == job_id:
            proc = _fusion_active["proc"]
            if proc.poll() is None:
                _terminate_fusion_proc(proc)
            _fusion_active = None
            _write_fusion_marker(job_id, "cancelled")
    return jsonify({"job_id": job_id, "status": "cancelled"})


# ---------------------------------------------------------------------------
# 页面
# ---------------------------------------------------------------------------
@app.route("/")
def landing() -> str:
    """站点首页：四个功能区入口（可见性 / 轨道可视化 / 发射窗口 / 轨道演示）。"""
    return render_template("landing.html")


@app.route("/visibility")
def index() -> str:
    return render_template("index.html")


@app.route("/orbits")
def orbits_page() -> str:
    return render_template("orbits.html")


@app.route("/elements_demo")
def elements_demo_page():
    """轨道六根数交互演示（开普勒二体 + 数值摄动 + CR3BP，纯前端计算）。

    内联脚本已外置到 static/elements_demo.js 以符合 CSP（script-src 'self'）；
    页内相对路径 static/... 落到 Flask 全局 static，贴图经 /assets/<path> 提供。
    """
    return send_file(os.path.join(BASE_DIR, "elements_demo.html"))


@app.route("/assets/<path:filename>")
def elements_demo_asset(filename: str):
    """elements_demo 贴图资产（send_from_directory 自带路径穿越防护）。"""
    return send_from_directory(ASSETS_DIR, filename)


# ---------------------------------------------------------------------------
# coast_tool 静态子应用（send_from_directory 自带路径穿越防护）
# ---------------------------------------------------------------------------
@app.route("/coast_tool/")
def coast_tool_index():
    if not os.path.isfile(os.path.join(COAST_TOOL_DIR, "index.html")):
        return jsonify({"error": "coast_tool 资产缺失"}), 503
    return send_from_directory(COAST_TOOL_DIR, "index.html")


@app.route("/coast_tool/<path:filename>")
def coast_tool_asset(filename: str):
    """coast_tool 静态资源（vendor/assets/app.js）；非法路径返回 404。"""
    return send_from_directory(COAST_TOOL_DIR, filename)


# ---------------------------------------------------------------------------
# API: Yang 论文轨道（预生成资产，只读）
# ---------------------------------------------------------------------------
@app.route("/api/yang_orbits")
@_rate_limit(30, 60.0)  # 每 IP 每分钟最多 30 次
def yang_orbits():
    """Yang 论文 24 条转移轨道 XYZ JSON（``t_arrive`` 仅支持论文默认历元）。"""
    t_arrive = _parse_t_arrive()
    data = _yang_orbits_data()
    if data is None:
        return jsonify({"error": "轨道数据资产不可用"}), 503
    if t_arrive != data.get("t_arrive_utc"):
        return (
            jsonify(
                {
                    "error": "当前部署仅提供论文默认到达历元 "
                    f"{data.get('t_arrive_utc')}（可用 scripts/build_yang_orbits.py 重新生成）"
                }
            ),
            400,
        )
    body = _yang_orbits_json_bytes()
    if body is None:  # pragma: no cover - 竞态兜底
        return jsonify({"error": "轨道数据资产不可用"}), 503
    return send_file(io.BytesIO(body), mimetype="application/json", max_age=3600)


@app.route("/api/yang_orbits_stream")
@_rate_limit(20, 60.0)  # SSE 连接较贵，限 20 次/分钟
def yang_orbits_stream():
    """SSE 逐轨道推送：``meta`` → 24 × ``orbit`` → ``done``。

    与 ``/api/yang_orbits`` 同源同数据，仅按事件分片，便于前端边收边渲染
    并展示进度。计算已离线完成，本端点不做重计算。
    """
    t_arrive = _parse_t_arrive()
    data = _yang_orbits_data()
    if data is None:
        return jsonify({"error": "轨道数据资产不可用"}), 503
    if t_arrive != data.get("t_arrive_utc"):
        return (
            jsonify(
                {
                    "error": "当前部署仅提供论文默认到达历元 "
                    f"{data.get('t_arrive_utc')}"
                }
            ),
            400,
        )

    meta = _yang_orbits_meta(data)
    orbits = data["orbits"]

    def generate():
        yield "event: meta\ndata: " + json.dumps(
            meta, ensure_ascii=False, allow_nan=False
        ) + "\n\n"
        for rec in orbits:
            yield "event: orbit\ndata: " + json.dumps(
                rec, ensure_ascii=False, allow_nan=False
            ) + "\n\n"
        yield "event: done\ndata: {}\n\n"

    resp = Response(stream_with_context(generate()), mimetype="text/event-stream")
    resp.headers["Cache-Control"] = "no-cache"
    # 反向代理（Render/nginx）禁用缓冲，保证事件即时到达前端
    resp.headers["X-Accel-Buffering"] = "no"
    return resp


# ---------------------------------------------------------------------------
# API: DEM 瓦片
# ---------------------------------------------------------------------------
@app.route("/api/dem_tile")
@_slow_endpoint_guard()  # 慢端点并发限制，超限返回 503
def dem_tile():
    if not dem_available():
        return jsonify({"error": f"DEM 文件不可用: {DEM_PATH}"}), 503
    center_lat = _parse_float_arg("center_lat", DEFAULT_LAT_DEG)
    center_lon = _parse_float_arg("center_lon", DEFAULT_LON_DEG)
    if not (-90.0 <= center_lat <= 90.0):
        raise BadRequest("center_lat 必须在 -90..90 范围内")
    if not (-180.0 <= center_lon <= 180.0):
        raise BadRequest("center_lon 必须在 -180..180 范围内")
    width_m = _parse_float_arg("width_m", 20000.0)
    height_m = _parse_float_arg("height_m", 20000.0)
    # 限制瓦片尺寸在合理范围（100 m – 200 km），防止超大窗口读取/OOM
    width_m = min(max(width_m, 100.0), 200000.0)
    height_m = min(max(height_m, 100.0), 200000.0)
    # 光照瓦片：可选时间参数（date=YYYY-MM-DD, hour=0-23.999）。缺省用默认日期 12:00。
    _t_utc, t_mjd, _date = _parse_time()

    try:
        png_bytes, meta = cached_dem_tile_png(
            center_lat, center_lon, width_m, height_m, t_mjd
        )
    except Exception as exc:  # RasterioIOError / WindowError 等
        return jsonify({"error": f"DEM 瓦片渲染失败: {type(exc).__name__}: {exc}"}), 500
    resp = send_file(
        io.BytesIO(png_bytes),
        mimetype="image/png",
        max_age=3600,
    )
    for k, v in meta.items():
        resp.headers[k] = v
    return resp


# ---------------------------------------------------------------------------
# API: 3D 地形 GLTF（第一人称天空视角）
# ---------------------------------------------------------------------------
@app.route("/api/terrain_gltf")
@_slow_endpoint_guard()  # 慢端点并发限制，超限返回 503
def terrain_gltf():
    if not dem_available():
        return jsonify({"error": f"DEM 文件不可用: {DEM_PATH}"}), 503
    lat = round(_parse_float_arg("lat", SHACKLETON_LAT_DEG), 4)
    lon = round(_parse_float_arg("lon", SHACKLETON_LON_DEG), 4)
    if not (-90.0 <= lat <= 90.0):
        raise BadRequest("lat 必须在 -90..90 范围内")
    if not (-180.0 <= lon <= 180.0):
        raise BadRequest("lon 必须在 -180..180 范围内")
    try:
        gltf_bytes, meta = cached_terrain_gltf(lat, lon)
    except (RuntimeError, ValueError) as exc:
        return jsonify({"error": str(exc)}), 400
    resp = send_file(
        io.BytesIO(gltf_bytes),
        mimetype="model/gltf+json",
        max_age=3600,
    )
    for k, v in meta.items():
        resp.headers[k] = v
    return resp


# ---------------------------------------------------------------------------
# API: 地平线掩码
# ---------------------------------------------------------------------------
@app.route("/api/horizon_mask")
@_slow_endpoint_guard()  # 慢端点并发限制，超限返回 503
def horizon_mask():
    if not dem_available():
        return jsonify({"error": f"DEM 文件不可用: {DEM_PATH}"}), 503
    lat = round(_parse_float_arg("lat", SHACKLETON_LAT_DEG), 4)
    lon = round(_parse_float_arg("lon", SHACKLETON_LON_DEG), 4)
    if not (-90.0 <= lat <= 90.0):
        raise BadRequest("lat 必须在 -90..90 范围内")
    if not (-180.0 <= lon <= 180.0):
        raise BadRequest("lon 必须在 -180..180 范围内")
    try:
        result = cached_horizon_mask(lat, lon)
    except (RuntimeError, ValueError) as exc:
        return jsonify({"error": str(exc)}), 400
    return jsonify(
        {
            "lat_deg": lat,
            "lon_deg": lon,
            "azimuths_deg": result["azimuths_deg"],
            "horizon_elev_deg": result["horizon_elev_deg"],
            "max_elev_deg": result["max_elev_deg"],
            "mean_elev_deg": result["mean_elev_deg"],
            "elev_m": result["elev_m"],
        }
    )


# ---------------------------------------------------------------------------
# API: 地球高度角（单时刻）
# ---------------------------------------------------------------------------
@app.route("/api/earth_elevation")
def earth_elevation():
    lat, lon, h = _parse_location()
    t_utc, t_mjd, _date = _parse_time()

    earth_pos_m = _earth_pos_iau_moon_m(t_mjd)
    elev, az = u.enu_elevation_azimuth(earth_pos_m, lon, lat, h, R_EQ, R_POL)
    return jsonify(
        {
            "elev_deg": round(float(elev), 4),
            "az_deg": round(float(az), 4),
            "visible": bool(elev > 0.0),
            "t_utc": t_utc,
            "lat_deg": lat,
            "lon_deg": lon,
        }
    )


# ---------------------------------------------------------------------------
# API: 地球高度角（日序列 288 点）
# ---------------------------------------------------------------------------
@app.route("/api/earth_daily")
def earth_daily():
    lat, lon, h = _parse_location()
    _t_utc, _t_mjd, date = _parse_time()

    t_mjd_start = u.utc_iso_to_mjd(f"{date}T00:00:00")
    t_mjd_arr = t_mjd_start + np.arange(DAILY_N_POINTS) * (
        DAILY_INTERVAL_MIN / 60.0 / 24.0
    )

    points = []
    for t in t_mjd_arr:
        earth_pos_m = _earth_pos_iau_moon_m(float(t))
        elev, az = u.enu_elevation_azimuth(earth_pos_m, lon, lat, h, R_EQ, R_POL)
        points.append(
            {
                "hour": float((t - t_mjd_start) * 24.0),
                "elev_deg": round(float(elev), 4),
                "az_deg": round(float(az), 4),
                "visible": bool(elev > 0.0),
            }
        )
    return jsonify({"date": date, "lat_deg": lat, "lon_deg": lon, "points": points})


# ---------------------------------------------------------------------------
# API: 中继卫星高度角 + 方位角（单时刻）
# ---------------------------------------------------------------------------
@app.route("/api/relay_elevation")
def relay_elevation():
    import spiceypy

    lat, lon, h = _parse_location()
    t_utc, t_mjd, _date = _parse_time()
    et = _mjd_to_et(t_mjd)

    # 与 relay_elevation.py 一致：轨道初值在固定任务历元构造，连续传播到查询时刻
    et_mission, state0, mu = _relay_mission_state()
    pos_j2000_km = np.asarray(
        propagate_satellite(state0, mu, et_mission, np.array([et]))[0], dtype=float
    )
    rot = spiceypy.pxform("J2000", "IAU_MOON", et)
    pos_iau_m = spiceypy.mxv(rot, pos_j2000_km) * 1000.0

    elev, az = u.enu_elevation_azimuth(pos_iau_m, lon, lat, h, R_EQ, R_POL)
    return jsonify(
        {
            "elev_deg": round(float(elev), 3),
            "az_deg": round(float(az), 3),
            "visible": bool(elev > 0.0),
            "t_utc": t_utc,
            "lat_deg": lat,
            "lon_deg": lon,
        }
    )


# ---------------------------------------------------------------------------
# API: 中继卫星（日序列 288 点）
# ---------------------------------------------------------------------------
@app.route("/api/relay_daily")
def relay_daily():
    lat, lon, h = _parse_location()
    _t_utc, _t_mjd, date = _parse_time()

    et_list, pos_iau_m = _relay_daily_positions(date)
    t_mjd_start = u.utc_iso_to_mjd(f"{date}T00:00:00")

    points = []
    for i, p in enumerate(pos_iau_m):
        elev, az = u.enu_elevation_azimuth(p, lon, lat, h, R_EQ, R_POL)
        points.append(
            {
                "hour": i * DAILY_INTERVAL_MIN / 60.0,
                "elev_deg": round(float(elev), 3),
                "az_deg": round(float(az), 3),
                "visible": bool(elev > 0.0),
            }
        )
    return jsonify({"date": date, "lat_deg": lat, "lon_deg": lon, "points": points})


# ---------------------------------------------------------------------------
# API: 极射投影 → 经纬度
# ---------------------------------------------------------------------------
@app.route("/api/stereo_to_lonlat")
def stereo_to_lonlat():
    x = _parse_float_arg("x", 0.0)
    y = _parse_float_arg("y", 0.0)
    lon, lat = u.polar_stereo_xy_to_lonlat(x, y)
    return jsonify({"lon": round(float(lon), 6), "lat": round(float(lat), 6)})


# ---------------------------------------------------------------------------
# API: 访问统计（按 IP + 日期粒度，CSV 持久化）
# ---------------------------------------------------------------------------
@app.route("/api/visit", methods=["POST"])
@_rate_limit(30, 60.0)  # 每 IP 每分钟最多 30 次访问记录
def api_visit():
    """记录一次页面访问。前端页面加载后调用；返回 204。"""
    ip = _client_ip()
    path = request.args.get("path", "/")
    record_visit(ip, path[:256])
    return "", 204


@app.route("/api/click", methods=["POST"])
@_rate_limit(30, 60.0)  # 每 IP 每分钟最多 30 次点击记录
def api_click():
    """记录一次按钮点击。按钮名经白名单过滤（防乱码/超长）。"""
    ip = _client_ip()
    button = request.args.get("button", "")
    button = (button or "").strip()[:64] or "unknown"
    record_click(ip, button)
    return "", 204


@app.route("/api/stats")
@_rate_limit(60, 60.0)  # 每 IP 每分钟最多 60 次统计查询
def api_stats():
    """按日期聚合的访问量 + 各按钮点击量（供前端统计面板展示）。"""
    return jsonify(stats_summary())


# ---------------------------------------------------------------------------
# API: 健康检查（smoke 测试用）
# ---------------------------------------------------------------------------
@app.route("/api/health")
def health():
    return jsonify(
        {
            "status": "ok",
            "dem_exists": dem_available(),  # 仅返回布尔值，不暴露 dem_path
            "spice_ready": _SPICE_READY,
            "orbit_data": os.path.isfile(YANG_ORBITS_JSON),
            "coast_tool": os.path.isfile(os.path.join(COAST_TOOL_DIR, "index.html")),
            "fusion_ready": bool(_fusion_supported_T())
            and os.path.isfile(FUSION_SCRIPT),
            "r_eq_m": round(float(R_EQ), 1),
            "r_pol_m": round(float(R_POL), 1),
        }
    )


# ---------------------------------------------------------------------------
# 全局错误处理（非法输入 → 400 JSON 而非 500；统一 API 错误格式）
# ---------------------------------------------------------------------------
@app.errorhandler(HTTPException)
def _handle_http_exception(e: HTTPException):
    """所有 HTTP 客户端错误统一返回 JSON，避免 HTML 错误页与信息泄露。"""
    return jsonify({"error": e.description}), e.code


@app.errorhandler(500)
def _handle_server_error(e):
    """未捕获异常统一返回 JSON 500，不泄露堆栈（debug=False 已保证）。"""
    return jsonify({"error": "内部错误"}), 500


if __name__ == "__main__":
    init_spice()
    # 端口可配置：环境变量 PORT / VIEWER_PORT（默认 5000；macOS 上 5000 可能被
    # AirPlay 接收器占用，可用 PORT=5050 等避开）
    port = int(os.environ.get("PORT", os.environ.get("VIEWER_PORT", "5050")))
    # 监听地址可配置：HOST 环境变量（默认 0.0.0.0 便于容器/云平台访问）
    host = os.environ.get("HOST", "0.0.0.0")
    print(
        f"SPICE kernels loaded. Moon radii: R_eq={R_EQ:.0f} m, "
        f"R_pol={R_POL:.0f} m"
    )
    print(f"DEM: {DEM_PATH} (exists={os.path.exists(DEM_PATH)})")
    print(f"Serving http://{host}:{port}")
    # 单线程：SPICE 非线程安全，Flask 开发服务器默认单线程即可。
    # 生产环境请改用 gunicorn/waitress（见模块 docstring）。
    app.run(host=host, port=port, debug=False, threaded=False)
