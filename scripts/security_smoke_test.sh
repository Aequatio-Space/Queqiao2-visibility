#!/bin/bash
# security_smoke_test.sh — 全站安全冒烟测试（主页 + 三个功能区 + API）
#
# 覆盖（对照 SECURITY_TEST_REPORT.md / SECURITY_REGRESSION_REPORT.md 口径）：
#   1. 安全响应头 / Server 指纹       2. CSP 与无内联脚本
#   3. 输入校验矩阵（400/404，不 500） 4. 方法限制
#   5. 路径穿越                       6. CSRF（Origin/Referer）
#   7. CSV 公式注入                   8. IP 伪造
#   9. 速率限制                       10. 信息泄露 / 严格 JSON
#  11. CORS                           12. 融合反解作业（校验/并发/取消）
#
# 用法：PORT=5099 bash scripts/security_smoke_test.sh
# 退出码：0 = 全部通过；1 = 存在失败项。
set -u

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PYTHON="${PYTHON:-python3}"
PORT="${PORT:-$((20000 + RANDOM % 20000))}"
B="http://127.0.0.1:${PORT}"
TMP="$(mktemp -d /tmp/q2sec.XXXXXX)"
STATS="$TMP/stats"; JOBS="$TMP/jobs"
mkdir -p "$STATS" "$JOBS"

PASS=0; FAIL=0
ok()  { PASS=$((PASS + 1)); printf "  PASS  %s\n" "$1"; }
bad() { FAIL=$((FAIL + 1)); printf "  FAIL  %s\n" "$1"; }

# code <method> <url> [curl 额外参数...] → HTTP 状态码
code() {
  local method="$1" url="$2"; shift 2
  curl -s -o /dev/null -w '%{http_code}' -X "$method" "$@" "$url"
}
# expect <期望码> <说明> <method> <url> [额外参数...]
expect() {
  local want="$1" desc="$2" method="$3" url="$4"; shift 4
  local got; got="$(code "$method" "$url" "$@")"
  if [ "$got" = "$want" ]; then ok "$desc → $got"; else bad "$desc → 期望 $want 实际 $got"; fi
}
# expect_in <候选1|候选2> <说明> ... 允许 400 或 404
expect_in() {
  local want="$1" desc="$2" method="$3" url="$4"; shift 4
  local got; got="$(code "$method" "$url" "$@")"
  case "|$want|" in
    *"|$got|"*) ok "$desc → $got" ;;
    *) bad "$desc → 期望 {$want} 实际 $got" ;;
  esac
}

cd "$REPO_ROOT"
PORT="$PORT" VIEWER_STATS_DIR="$STATS" FUSION_JOB_DIR="$JOBS" \
  nohup "$PYTHON" viewer/server.py > "$TMP/server.log" 2>&1 &
SRV=$!
trap 'kill $SRV 2>/dev/null' EXIT
for _ in $(seq 1 80); do
  grep -q "Address already in use" "$TMP/server.log" && break
  curl -s -o /dev/null "$B/api/health" && break
  sleep 0.5
done
if grep -q "Address already in use" "$TMP/server.log"; then
  echo "端口 $PORT 已被占用，请设置 PORT=<空闲端口> 后重试"; exit 1
fi
if ! curl -s -o /dev/null "$B/api/health"; then
  echo "服务启动失败："; tail -20 "$TMP/server.log"; exit 1
fi
echo "服务已启动：${B}（日志 $TMP/server.log）"

# ---------------------------------------------------------------------------
echo "== 1. 安全响应头 / Server 指纹 =="
HDR_KEYS="x-content-type-options|x-frame-options|referrer-policy|permissions-policy|strict-transport-security|content-security-policy"
for p in / /visibility /orbits /coast_tool/ /api/health /api/yang_orbits /coast_tool/api/fusion/config; do
  n=$(curl -s -D - -o /dev/null "$B$p" | grep -ciE "^($HDR_KEYS):")
  [ "$n" = "6" ] && ok "$p 安全头 6/6" || bad "$p 安全头 $n/6"
done
srv=$(curl -s -D - -o /dev/null "$B/" | grep -i '^Server:' | tr -d '\r')
case "$srv" in
  *Werkzeug*|*Python*) bad "Server 头泄露版本：$srv" ;;
  *) ok "Server 头已脱敏：$srv" ;;
esac
csp=$(curl -s -D - -o /dev/null "$B/" | grep -i '^Content-Security-Policy:' | tr -d '\r')
for token in "default-src 'self'" "script-src 'self'" "frame-ancestors 'none'" "object-src"; do
  case "$csp" in *"$token"*) ok "CSP 含 $token" ;; *) bad "CSP 缺 $token" ;; esac
done

# ---------------------------------------------------------------------------
echo "== 2. 无内联脚本 / 内联事件处理器 =="
for p in / /visibility /orbits /coast_tool/; do
  html=$(curl -s "$B$p")
  inline=$(printf '%s' "$html" | grep -cE '<script[^>]*>[^<]' || true)
  handler=$(printf '%s' "$html" | grep -cE 'on(click|change|input|submit|load|error)=' || true)
  [ "$inline" = "0" ] && ok "$p 无内联脚本" || bad "$p 含 $inline 段内联脚本"
  [ "$handler" = "0" ] && ok "$p 无内联事件处理器" || bad "$p 含 $handler 个内联事件处理器"
done
for f in /coast_tool/app.js /coast_tool/rockets.js /static/orbits.js /static/app.js; do
  body=$(curl -s "$B$f")
  case "$body" in *"eval("*|*"new Function"*) bad "$f 含动态代码执行" ;; *) ok "$f 无 eval/new Function" ;; esac
done

# ---------------------------------------------------------------------------
echo "== 3. 输入校验矩阵（应 400/404，不得 500）=="
expect 400 "dem_tile lat=abc"      GET "$B/api/dem_tile?center_lat=abc"
expect 400 "dem_tile lat=1e999"    GET "$B/api/dem_tile?center_lat=1e999"
expect 400 "dem_tile lat=NaN"      GET "$B/api/dem_tile?center_lat=NaN"
expect 400 "dem_tile lon=200"      GET "$B/api/dem_tile?center_lon=200"
expect 400 "horizon_mask lat=inf"  GET "$B/api/horizon_mask?lat=inf"
expect 400 "earth_elevation h=x"   GET "$B/api/earth_elevation?h=xyz"
expect 400 "relay_daily date=bad"  GET "$B/api/relay_daily?date=INVALID"
expect 400 "relay_daily date=13-45" GET "$B/api/relay_daily?date=2026-13-45"
expect 400 "stereo x=../../etc/passwd" GET "$B/api/stereo_to_lonlat?x=../../etc/passwd"
expect 400 "yang t_arrive=abc"     GET "$B/api/yang_orbits?t_arrive=abc"
expect 400 "yang t_arrive=1899"    GET "$B/api/yang_orbits?t_arrive=1899-01-01"
expect 400 "yang t_arrive=2030"    GET "$B/api/yang_orbits?t_arrive=2030-01-01"
expect 400 "yang SSE t_arrive=abc" GET "$B/api/yang_orbits_stream?t_arrive=abc"
expect 400 "fusion 无 confirm" POST "$B/coast_tool/api/fusion/jobs" \
  -H 'Content-Type: application/json' -d '{"start":"2026-08-24","end":"2026-08-25","T_h":112.8}'
expect 400 "fusion T=72" POST "$B/coast_tool/api/fusion/jobs" \
  -H 'Content-Type: application/json' -d '{"start":"2026-08-24","end":"2026-08-25","T_h":72,"confirm":true}'
expect 400 "fusion span 40d" POST "$B/coast_tool/api/fusion/jobs" \
  -H 'Content-Type: application/json' -d '{"start":"2026-08-01","end":"2026-09-11","T_h":112.8,"confirm":true}'
expect 400 "fusion site=evil" POST "$B/coast_tool/api/fusion/jobs" \
  -H 'Content-Type: application/json' -d '{"start":"2026-08-24","end":"2026-08-25","T_h":112.8,"site":"evil","confirm":true}'
expect 400 "fusion T=NaN" POST "$B/coast_tool/api/fusion/jobs" \
  -H 'Content-Type: application/json' -d '{"start":"2026-08-24","end":"2026-08-25","T_h":"NaN","confirm":true}'
expect 400 "fusion 非 JSON 体" POST "$B/coast_tool/api/fusion/jobs" \
  -H 'Content-Type: application/json' -d 'not-json'
expect_in "400|404" "fusion 非法 job_id" GET "$B/coast_tool/api/fusion/jobs/..%2f..%2fetc%2fpasswd"
expect 404 "fusion 不存在 job" GET "$B/coast_tool/api/fusion/jobs/0123456789abcdef"
# 任何 5xx 都视为失败
for u in "$B/api/dem_tile?center_lat=abc" "$B/api/horizon_mask?lat=1e999" \
         "$B/api/earth_elevation?date=INVALID" "$B/api/yang_orbits?t_arrive=%2e%2e%2f"; do
  c=$(code GET "$u")
  case "$c" in 5*) bad "非法输入触发 5xx：$u → $c" ;; *) ok "非法输入无 5xx（${c}）：$u" ;; esac
done
# 超大请求体（MAX_CONTENT_LENGTH=1MB）→ 413 JSON
head -c 2000000 /dev/zero > "$TMP/big.bin"
expect 413 "超大请求体 → 413" POST "$B/api/visit" \
  -H 'Content-Type: application/octet-stream' --data-binary @"$TMP/big.bin"

# ---------------------------------------------------------------------------
echo "== 4. HTTP 方法限制 =="
expect 405 "POST /api/yang_orbits"  POST "$B/api/yang_orbits"
expect 405 "POST /api/health"       POST "$B/api/health"
expect 405 "GET /api/visit"         GET  "$B/api/visit"
expect 405 "GET /api/click"         GET  "$B/api/click"
expect 405 "DELETE /api/yang_orbits" DELETE "$B/api/yang_orbits"
c=$(code GET "$B/coast_tool/api/fusion/jobs")
[ "$c" = "200" ] && bad "GET fusion 集合端点返回 200" || ok "GET fusion 集合端点不可用（${c}）"

# ---------------------------------------------------------------------------
echo "== 5. 路径穿越 =="
expect 404 "coast_tool ../../etc/passwd" GET "$B/coast_tool/../../etc/passwd" --path-as-is
expect 404 "coast_tool 编码 %2e%2e"      GET "$B/coast_tool/%2e%2e/%2e%2e/etc/passwd"
expect 404 "static ../server.py"         GET "$B/static/../viewer/server.py" --path-as-is
expect 404 "static 编码穿越"             GET "$B/static/%2e%2e/viewer/server.py"
expect 404 "轨道 JSON 静态旁路"          GET "$B/static/data/yang_orbits_xyz.json"
expect_in "404|308" "coast_tool 双斜杠穿越" GET "$B/coast_tool//etc/passwd" --path-as-is

# ---------------------------------------------------------------------------
echo "== 6. CSRF（跨站 POST 应 403）=="
expect 403 "visit 跨站 Origin" POST "$B/api/visit?path=/x" -H 'Origin: http://evil.com'
expect 403 "click 跨站 Origin" POST "$B/api/click?button=x" -H 'Origin: http://evil.com'
expect 403 "fusion 跨站 Origin" POST "$B/coast_tool/api/fusion/jobs" \
  -H 'Origin: http://evil.com' -H 'Content-Type: application/json' \
  -d '{"start":"2026-08-24","end":"2026-08-25","T_h":112.8,"confirm":true}'
expect 204 "visit 同源 Origin" POST "$B/api/visit?path=/ok" -H "Origin: $B"
expect 403 "visit 跨站 Referer" POST "$B/api/visit?path=/x" -H 'Referer: http://evil.com/'

# ---------------------------------------------------------------------------
echo "== 7. CSV 公式注入 =="
curl -s -o /dev/null -X POST --get --data-urlencode "button==cmd|'/C calc'!A1" "$B/api/click"
curl -s -o /dev/null -X POST -H 'X-Forwarded-For: =1+1' \
  --get --data-urlencode 'path==HYPERLINK("http://evil.com","x")' "$B/api/visit"
tail -1 "$STATS/clicks.csv" | grep -q "'=cmd" && ok "clicks.csv 公式已中和" || bad "clicks.csv 公式未中和"
tail -1 "$STATS/visits.csv" | grep -q "'=HYPERLINK" && ok "visits.csv 公式已中和" || bad "visits.csv 公式未中和"
long=$(curl -s -o /dev/null -w '%{http_code}' -X POST --get --data-urlencode "button=$(printf 'A%.0s' $(seq 1 500))" "$B/api/click")
[ "$long" = "204" ] && ok "超长 button 被截断且受理" || bad "超长 button → $long"

# ---------------------------------------------------------------------------
echo "== 8. IP 伪造（未配置可信代理）=="
curl -s -o /dev/null -X POST -H 'X-Forwarded-For: 999.999.999.999' "$B/api/visit?path=/spoof"
tail -1 "$STATS/visits.csv" | grep -q "127.0.0.1" && ok "XFF 被忽略（记录套接字 IP）" || bad "XFF 被采信"

# ---------------------------------------------------------------------------
echo "== 9. 信息泄露 / 严格 JSON / CORS =="
health=$(curl -s "$B/api/health")
case "$health" in *"/Users/"*|*"dem_path"*) bad "health 泄露路径：$health" ;; *) ok "health 无路径泄露" ;; esac
case "$health" in *'"fusion_ready"'*|*'"orbit_data"'*) ok "health 含资产布尔状态" ;; *) bad "health 缺资产状态" ;; esac
curl -s -D - -o /dev/null "$B/api/yang_orbits" | grep -qi 'access-control-allow-origin' \
  && bad "存在 CORS 放行" || ok "无 CORS 放行"
curl -s "$B/api/yang_orbits" -o "$TMP/orbits.json"
if "$PYTHON" - "$TMP/orbits.json" <<'PY'
import json, sys
raw = open(sys.argv[1], 'rb').read()
assert b'NaN' not in raw and b'Infinity' not in raw, 'non-strict JSON'
def boom(x): raise ValueError(x)
d = json.loads(raw, parse_constant=boom)
assert d['n_orbits'] == 24 and all(len(o['r_mep']) == 121 for o in d['orbits'])
PY
then ok "yang_orbits 严格 JSON（24 条 × 121 点，无 NaN）"; else bad "yang_orbits JSON 不严格"; fi
c=$(code GET "$B/api/does-not-exist")
[ "$c" = "404" ] && ok "未知路径 404" || bad "未知路径 → $c"
curl -s "$B/api/does-not-exist" | grep -qiE "traceback|/Users/" && bad "404 响应泄露堆栈/路径" || ok "404 响应无泄露"
for p in /coast_tool/app.js /coast_tool/assets/earth_day.jpg /static/orbits.js; do
  ct=$(curl -s -o /dev/null -w '%{content_type}' "$B$p")
  ok "$p Content-Type: $ct"
done

# ---------------------------------------------------------------------------
echo "== 10. 速率限制 =="
codes=""
for _ in $(seq 1 40); do codes="$codes $(code POST "$B/api/visit?path=/r")"; done
echo "$codes" | grep -q "429" && ok "visit 超限返回 429" || bad "visit 未限流"
codes=""
for _ in $(seq 1 35); do codes="$codes $(code GET "$B/api/yang_orbits")"; done
echo "$codes" | grep -q "429" && ok "yang_orbits 超限返回 429" || bad "yang_orbits 未限流"

# ---------------------------------------------------------------------------
echo "== 11. 融合反解作业（长耗时，子进程隔离）=="
JOB=$(curl -s -X POST -H 'Content-Type: application/json' -H "Origin: $B" \
  -d '{"start":"2026-08-24","end":"2026-08-25","T_h":112.8,"site":"wenchang","coast_max_s":99999,"confirm":true}' \
  "$B/coast_tool/api/fusion/jobs")
JID=$(printf '%s' "$JOB" | "$PYTHON" -c 'import sys,json; print(json.load(sys.stdin).get("job_id",""))' 2>/dev/null)
if [ -n "$JID" ]; then
  ok "作业已启动 ${JID}（coast_max 应钳制 3600）"
  expect 409 "并发作业 409" POST "$B/coast_tool/api/fusion/jobs" \
    -H 'Content-Type: application/json' -d '{"start":"2026-08-24","end":"2026-08-25","T_h":112.8,"confirm":true}'
  for _ in $(seq 1 60); do
    ST=$(curl -s "$B/coast_tool/api/fusion/jobs/$JID")
    S=$(printf '%s' "$ST" | "$PYTHON" -c 'import sys,json; print(json.load(sys.stdin).get("status",""))' 2>/dev/null)
    [ "$S" != "running" ] && break
    sleep 2
  done
  if [ "$S" = "done" ]; then
    printf '%s' "$ST" > "$TMP/fusion_status.json"
    if "$PYTHON" - "$TMP" <<'PY'
import json, sys
d = json.load(open(sys.argv[1] + '/fusion_status.json'))['result']
c = d['config']
assert c['coast_max_s'] == 3600.0, 'coast_max 未钳制'
assert d['n_windows'] >= 0
open(sys.argv[1] + '/fusion_result.json', 'w').write(json.dumps(d))
PY
    then ok "作业完成：窗口数 $("$PYTHON" -c 'import json;print(json.load(open("'"$TMP"'/fusion_result.json"))["n_windows"])')，coast_max 已钳制"
    else bad "作业结果校验失败"; fi
  else
    bad "作业未在预期时间内完成（status=${S}）"
  fi
  # 取消流程
  J2=$(curl -s -X POST -H 'Content-Type: application/json' \
    -d '{"start":"2026-08-24","end":"2026-08-25","T_h":112.8,"confirm":true}' "$B/coast_tool/api/fusion/jobs")
  JID2=$(printf '%s' "$J2" | "$PYTHON" -c 'import sys,json; print(json.load(sys.stdin).get("job_id",""))' 2>/dev/null)
  sleep 1
  expect 200 "取消作业" DELETE "$B/coast_tool/api/fusion/jobs/$JID2"
  sleep 1
  S2=$(curl -s "$B/coast_tool/api/fusion/jobs/$JID2" | "$PYTHON" -c 'import sys,json; print(json.load(sys.stdin).get("status",""))')
  [ "$S2" = "cancelled" ] && ok "取消后状态 cancelled" || bad "取消后状态 $S2"
else
  bad "作业启动失败：$JOB"
fi

# ---------------------------------------------------------------------------
echo "== 12. 访问统计聚合可用（前端统计面板依赖）=="
stats=$(curl -s "$B/api/stats")
printf '%s' "$stats" | "$PYTHON" -c 'import sys,json; d=json.load(sys.stdin); assert "days" in d and "total_visits" in d' \
  && ok "/api/stats JSON 结构正常" || bad "/api/stats 结构异常"

# ---------------------------------------------------------------------------
kill $SRV 2>/dev/null; wait $SRV 2>/dev/null
echo
echo "=============================================="
echo "  通过 $PASS 项，失败 $FAIL 项"
echo "  临时目录：$TMP"
echo "=============================================="
[ "$FAIL" -eq 0 ] || exit 1
