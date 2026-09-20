# 多入口站点新增功能区 — 安全评审与回归报告

> **评审对象**：新增首页（`/`）、轨道可见性迁移（`/visibility`）、轨道 3D 可视化
> （`/orbits` + `/api/yang_orbits`、`/api/yang_orbits_stream`）、发射窗口静态子应用
> （`/coast_tool/`）
> **日期**：2026-09-21
> **环境**：本地 `127.0.0.1:5079`，conda `spaceport`，DEM + SPICE 内核就绪
> **方法**：延续 `SECURITY_TEST_REPORT.md` / `SECURITY_REGRESSION_REPORT.md`
> 的黑盒探测 + 白盒审计口径，共 12 组用例

---

## 一、结论

**新增功能区沿用原有安全基线，未引入新的高危/中危问题。**

| 维度 | 状态 |
| --- | --- |
| 安全响应头（6 项，含 CSP） | ✅ 所有新页面/接口 6/6 |
| 输入校验（`t_arrive`） | ✅ 非法/越界/注入串全部 400 |
| 路径穿越（`/coast_tool/` 静态服务） | ✅ 全部 404 |
| 方法限制 | ✅ 405 |
| CSRF（跨站 POST） | ✅ 403 |
| CSV 公式注入 | ✅ `'` 前缀中和 |
| IP 伪造 | ✅ 仅记录套接字地址 |
| 速率限制 | ✅ 第 30 次后 429（按端点独立计数） |
| CORS | ✅ 未放行 |
| 信息泄露（health） | ✅ 仅布尔值，无文件系统路径 |
| 严格 JSON | ✅ 无 `NaN`/`Infinity` |
| 供应链 | ✅ 无 CDN/外部请求，three.js 与贴图本地 vendored |

---

## 二、新增攻击面与处置

### 2.1 `/api/yang_orbits`、`/api/yang_orbits_stream`（轨道数据）

| 风险 | 处置 |
| --- | --- |
| 慢端点 DoS（一次完整求解约 20 s 单核） | **改为离线预生成资产**（`scripts/build_yang_orbits.py` → `viewer/data/yang_orbits_xyz.json`），请求内零重计算，从根上消除该 DoS 面 |
| 用户字符串直入 SPICE（`str2et` 注入/异常） | `_parse_t_arrive()` 仅接受 `YYYY-MM-DD` / `YYYY-MM-DDTHH:MM:SSZ`，且年份限 1900–2050；其余 400 |
| 任意历元导致缓存耗尽 | 仅服务论文默认历元（`2026-08-24T00:00:00Z`），非默认 400 |
| 非严格 JSON（`NaN`） | `_json_safe()` 递归净化 + `allow_nan=False`；`/api/yang_orbits` 实测无 `NaN`/`Infinity` |
| 高频拉取 | 速率限制：JSON 30/min、SSE 20/min（按端点独立计数桶） |
| 资产被静态目录旁路 | JSON 存于 `viewer/data/`（非 `static/`），`/static/data/...` 实测 404，仅 API 单一入口 |

### 2.2 `/coast_tool/`（静态子应用）

| 风险 | 处置 |
| --- | --- |
| 内联脚本违反 CSP | 原 4 段内联脚本按原顺序外置为 `app.js`（内容未改动），全站保持 `script-src 'self'` |
| 路径穿越 | Flask `send_from_directory`（`safe_join`）；`../`、编码 `%2e%2e` 实测 404 |
| 外部依赖/供应链 | three.js r147、OrbitControls、NASA 贴图全部仓库内本地文件；无 `fetch`/XHR/WebSocket/CDN |
| DOM XSS | 页面无内联事件处理器；`innerHTML` 仅用于内部数值格式化（无用户可控字符串） |
| 代码执行面 | `app.js` 无 `eval` / `new Function` / `document.write` |

### 2.3 页面与统计

- 首页与轨道页均无内联脚本/内联事件处理器，CSP 生效；
- 访问/点击统计复用既有 `/api/visit`、`/api/click`（CSV 中和、64 字符截断、可信代理校验、限流不变）；
- `_rate_limit` 增加 `scope` 参数：各端点独立计数桶，避免高频端点挤占其他端点配额（行为兼容：限流阈值与响应不变）。

---

## 三、回归证据（摘要）

```text
T1 页面可用性          /、/visibility、/orbits、/coast_tool/ 及 app.js/vendor/assets 全部 200
T2 安全响应头          各页面与 API 均 6/6（XFO/nosniff/Referrer/Permissions/HSTS/CSP）
T3 JSON 严格性         n_orbits=24，121 点/条，无 NaN/Infinity
T4 SSE 完整性          meta=1，orbit=24，done=1
T5 t_arrive 校验       abc / 2026-13-45 / 2030-01-01 / 1899-12-31 / 1e999 /
                       ../../etc/passwd / +08:00 / 空格 全部 400
T6 路径穿越            ../../etc/passwd → 404；%2e%2e → 404；静态旁路 → 404
T7 方法限制            POST /api/yang_orbits → 405；GET /api/visit → 405
T8 CSRF                Origin: evil.com → 403；同源 → 204
T9 CSV 注入            button='=cmd|... 落盘为 "'=cmd|...'；超长截断 64 字符
T10 IP 伪造            X-Forwarded-For: 999.999.999.999 被忽略，记录 127.0.0.1
T11 限流               35 连发（叠加前置用例消耗）→ 20×200 + 15×429，
                       配额（30/min）耗尽后全部 429
T12 CORS / health      ACAO 0 次；health 仅布尔状态，无 /Users/ 路径
```

复现命令见 `SECURITY_TEST_REPORT.md` 风格；测试脚本要点：以
`VIEWER_STATS_DIR=<tmp>` 隔离统计落盘，`--path-as-is` 绕过 curl 路径规范化，
`json.loads(parse_constant=raise)` 校验严格 JSON。

---

## 四、融合反向求解作业（coast_tool 可选工具，新增）

### 4.1 设计

- 反解为 CPU 密集长作业（单日 ~10 s、整月数分钟）。为不阻塞单线程 Web
  服务（SPICE 非线程安全），Web 层**不**在请求内计算，而是以子进程运行
  `scripts/fusion_inverse_job.py`（`start_new_session=True`，可整组回收）；
- Web 层仅做参数校验、作业调度、状态查询与结果读取；
- 同一时刻仅一个作业（全局锁 + 409），作业可取消，超时（默认 1 h）由
  状态端点看门狗终止，产物 48 h 后自动清理。

### 4.2 攻击面处置

| 风险 | 处置 |
| --- | --- |
| 命令注入 | 子进程参数由**校验后的数值/白名单**拼装（无 shell、`shell=False`）；站点白名单、日期 `strptime`、T 吸附到预生成表格、`coast_max_s` 钳制 0–3600 |
| 路径注入 | `job_id` 由服务端生成（16 位十六进制），查询/取消路径先 `fullmatch` 白名单校验；作业目录由 id 拼接，无用户路径输入 |
| 未确认即运行（资源滥用） | 请求体必须 `confirm=true`，否则 400；前端另有确认对话框 |
| 并发滥用 | 单作业槽位（409）+ 启动限流 10/10 min（独立计数桶） |
| 密钥泄漏给子进程 | 子进程环境**最小化**（仅 PATH/HOME/PYTHONUNBUFFERED/FUSION_TABLE_DIR/SPICE_KERNEL_DIR），不继承服务端环境变量 |
| 作业目录无限增长 | 启动时清理 >48 h 旧目录；单作业产物在作业目录内 |
| 状态接口信息泄露 | 仅返回状态/耗时/窗口数值与配置（无文件系统路径）；不返回原始日志 |
| CSRF | 复用全局 `before_request` Origin/Referer 校验（POST/DELETE 均覆盖） |
| 结果注入前端 | 前端以 `textContent` 渲染全部字段，无 `innerHTML` |

### 4.3 回归结果（2026-09-21，本地 5087/5088 端口）

```text
校验矩阵         no confirm / bad date / start>=end / span 40d / T=NaN / T=72 /
                 site=evil → 全部 400；coast_max_s=-5 钳制为 0 后正常受理
路径安全         job_id ../../etc/passwd（--path-as-is）→ 404；
                 %2e%2e%2f 编码穿越 → 404；不存在的合法 id → 404
CSRF             Origin: http://evil.com → 403
并发             作业运行中再次提交 → 409
限流             启动端点 10/10min 独立计数（测试中超额后 429）
吸附/钳制        T=112.75 → 112.7 h（命中预生成表）；coast_max_s=99999 → 3600
作业生命周期     2 天范围：6 窗口 / 16.6 s；取消后状态 cancelled；
                 result.json 含窗口列表 + 拒绝统计；子进程最小环境验证通过
```

### 4.4 残留风险（非阻塞）

1. **内存**：融合子进程加载 scipy + spiceypy（~150–250 MB）；与 Web 进程
   同机时需确保实例内存充足（Render free 512 MB 较紧张，建议更高配）。
2. **CPU 争用**：作业运行期间单 vCPU 实例的页面响应会变慢；确认对话框已
   明示，建议生产使用 ≥2 vCPU 或独立 worker。
3. **Serverless 不适用**：Vercel 实例在响应后冻结/回收，无法承载分钟级
   后台子进程；该功能需常驻服务（Render / Docker）。

---

## 五、残留风险（纵深防御建议，非阻塞）

1. **SSE 长连接**：单 worker 下每个 SSE 连接占用一个请求线程；已限流 20/min，
   建议反代层（Render/Cloudflare）设置单 IP 并发连接上限。
2. **静态子应用文档**：`/coast_tool/README.md` 可被公开读取（内容为公开算法
   说明），如不希望暴露可在路由中白名单扩展名。
3. **Vercel Serverless**：`/api/yang_orbits` 为纯静态读取，无超时风险；SSE 在
   部分 Serverless 网关可能被缓冲，建议以 Render 常驻服务为主入口。

---

*本报告基于 2026-09-21 本地模拟测试；新增功能未修改既有端点行为。*
