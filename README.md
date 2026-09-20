# Queqiao2-visibility — 月球探测可视化站

基于 Flask 的交互式月球探测可视化站点，首页提供三个功能区入口：

| 功能区 | 路径 | 内容 |
| --- | --- | --- |
| **轨道可见性分析** | `/visibility` | 月球南极极射投影 DEM 光照地形图（实时太阳光照渲染）、极坐标地平线图 / 第一人称 3D 天空视角，底部时间滑块联动展示地球与鹊桥二号（Queqiao-2）中继卫星的全天可见性窗口 |
| **地月转移轨道可视化** | `/orbits` | 杨维廉论文 24 条地月转移轨道（表 1–4）在冻结 MEP 坐标系中的 three.js 三维展示：三栏工作台、搜索过滤与族谱分组、**任务时间轴动画**（航天器沿采样轨迹推进、月球沿真实轨迹移动）、**C3 柱状图**（点击选择）、聚焦选中、相机聚焦与视角预设；地球采用 NASA Blue Marble 昼夜贴图（复用 `coast_tool/assets/` 素材）（`launch_window` 引擎离线求解，SSE 推送） |
| **月球发射窗口计算器** | `/coast_tool/` | 停泊轨道滑行方案（固定射向走廊 → ~200 km 停泊 → 滑行相位闭合 → TLI 入 LTO），四体引力修正支持**任意目标环月倾角 i_M ∈ [0°,180°]**；内置 **CZ-5/CZ-8 火箭能力核对**与**嫦娥五号/鹊桥二号/嫦娥六号三个验证案例**（展示滑行偏长等已知误差）；右侧另附**可选的融合反向求解**（高精度 L6 家族表反解，须确认后运行） |

数据由服务端计算（SPICE 天历 + 月球 DEM），前端全部基于 Canvas/WebGL
客户端渲染。DEM 使用仓库内置的**压缩子集** `ldem_87s_5mpp_shackleton_50km_25m.tif`
（沙克尔顿坑 50 km 半径、25 m/px、52 MB，随仓库分发，**无需持久盘**）；
如需更高精度可自行用完整 DEM（3.2 GB）覆盖 `DEM_PATH`。

站点完全自包含：SPICE 内核、three.js、地球贴图与轨道 JSON 资产均随仓库分发，
无任何外部 CDN / 外部服务依赖。

## 目录结构

```
viewer/                     # Flask 应用
  server.py                 # 后端（生产 WSGI 入口：viewer.server:app）
  templates/                # 页面：landing.html（首页）/ index.html（可见性）
                            #       orbits.html（轨道 3D）/ coast_tool 在仓库根
  static/app.js             # 可见性前端逻辑
  static/orbits.js          # 轨道 3D 前端逻辑
  static/landing.js         # 首页访问统计
  static/vendor/            # three.js / GLTFLoader（本地 vendored）
  data/yang_orbits_xyz.json # Yang 论文 24 条轨道 XYZ 资产（预生成，随仓库分发）
  data/fusion_tables/       # 融合反解 L5 家族表（T=112.0–113.0 h，预生成）
  jobs/                     # 融合反解作业产物（运行时生成，git-ignored）
  stats/                    # 访问统计 CSV（运行时生成，git-ignored）
launch_window/              # 移植的轨道/窗口引擎（离线轨道资产 + 融合反解）
coast_tool/                 # 月球发射窗口计算器（静态子应用，three.js vendored）
scripts/build_yang_orbits.py# 重新生成 viewer/data/yang_orbits_xyz.json
scripts/fusion_inverse_job.py# 融合反向求解后台作业（由 Web 层以子进程调用）
analysis/elevation_video/   # 科学计算模块（地平线掩码 / 轨道传播 / ENU 投影 / 中文字体）
core/environment/data/      # SPICE 内核（随仓库分发）
api/index.py                # Vercel Serverless 入口
render.yaml                 # Render Blueprint
vercel.json                 # Vercel 配置（全路径转发至 Flask）
Procfile                    # Render / 通用平台进程定义
requirements.txt            # Python 依赖
```

功能亮点：
- 首页三入口导航（可见性 / 轨道可视化 / 发射窗口），站点内互相跳转；
- 左侧月球光照地形图（实时太阳光照渲染）、右侧极坐标地平线图 / 第一人称 3D 天空视角；
- 底部时间滑块联动地球与鹊桥二号全天可见性窗口；
- 杨维廉论文 24 条转移轨道 3D 可视化（共面 / 非共面 / 同向 / 极月四族）；
- **访问统计**：页面访问 + 按钮点击按 IP/日期落盘 CSV，前端 📊 统计面板查看。

## 本地运行

要求 Python 3.10+（推荐 conda/mamba 环境）。DEM 使用仓库内置的压缩子集
`ldem_87s_5mpp_shackleton_50km_25m.tif`（随仓库分发），**无需额外下载或持久盘**。

```bash
# 1. 安装依赖
pip install -r requirements.txt

# 2. 启动（开发服务器，默认即用内置 DEM）
python viewer/server.py
# 或指定 DEM 路径 / 端口
DEM_PATH=/path/to/ldem_87s_5mpp.tif PORT=5050 python viewer/server.py

# 3. 浏览器打开
#    http://127.0.0.1:5050            （首页：三入口导航）
#    http://127.0.0.1:5050/visibility （轨道可见性）
#    http://127.0.0.1:5050/orbits     （轨道 3D 可视化）
#    http://127.0.0.1:5050/coast_tool/（发射窗口计算器）
```

轨道 3D 页面读取仓库内置资产 `viewer/data/yang_orbits_xyz.json`（随仓库分发，
无需运行时求解）。页面功能：左侧按族谱分组/搜索/勾选显示，点击轨道查看
p/q 解、C3、近地点/近月点经纬度等详情；底部时间轴可播放全弧段动画
（1×–8×、逐采样步进、空格播放/暂停），航天器与月球按 121 个采样点运动；
右栏 C3 柱状图可直接点选轨道。如需在**其他到达历元**重建 24 条论文轨道：

```bash
# 约 20 s 单核 CPU；需要 scipy + spiceypy（已在 requirements.txt）
python scripts/build_yang_orbits.py            # 写入 viewer/data/yang_orbits_xyz.json
```

### 融合反向求解（coast_tool 可选工具）

`/coast_tool/` 右侧提供可选的高精度**融合反向求解**（`launch_window/fusion.py`
的 FR-4 L6 反解链路）。该计算为 CPU 密集长作业（单日约 10–40 s、整月约
5–15 min），因此：

- 前端须在确认对话框中确认后才提交；
- 后端以**独立子进程**运行 `scripts/fusion_inverse_job.py`（SPICE 非线程安全，
  且避免阻塞单线程 Web 服务），前端轮询作业状态；
- 同一时刻只允许一个作业（重复提交返回 409），作业可取消，超时（默认 1 h）
  自动终止，产物 48 h 后清理。

支持范围（由预生成家族表决定，其余参数取反解默认值）：

| 参数 | 支持值 |
| --- | --- |
| 转移时间 `T` | 112.0–113.0 h（0.1 h 步长）与 112.61 h，自动吸附 ±0.35 h |
| 目标环月倾角 `i_m` | 90°（极轨捕获，月球真极轴参考 `moon_equator`） |
| 近月点高度 `H_m` | 200 km（含太阳摄动） |
| 扫描范围 | 1–31 天（右开区间） |
| 发射场 | wenchang / xichang / ksc |

> 反解链路依赖 `viewer/data/fusion_tables/` 的预生成 L5 家族表（随仓库分发，
> 约 2.6 MB）；如需其他 `T`/`i_m`，可用上游 `scan_im.run_scan` 离线重建后放入
> 该目录（不支持的参数会在前端直接提示）。

生产 WSGI（gunicorn / waitress）：

```bash
gunicorn viewer.server:app --bind 0.0.0.0:$PORT --workers 1 --threads 1 --timeout 300
# 或
waitress-serve --host 0.0.0.0 --port $PORT viewer.server:app
```

> 单线程（`--workers 1 --threads 1`）：SPICE（spiceypy）非线程安全。

## 安全测试

`scripts/security_smoke_test.sh` 对主页与三个功能区及全部 API 做 12 组 88 项
安全检查（安全头/CSP/无内联脚本、输入校验矩阵、方法限制、路径穿越、CSRF、
CSV 注入、IP 伪造、限流、信息泄露、严格 JSON、CORS、融合反解作业）：

```bash
PYTHON=python3 bash scripts/security_smoke_test.sh   # 默认随机空闲端口
# 或指定端口 / 解释器
PORT=5099 PYTHON=/path/to/python bash scripts/security_smoke_test.sh
```

退出码 0 表示全部通过。测试结果与历次报告见
`SECURITY_TEST_REPORT.md`、`SECURITY_REGRESSION_REPORT.md`、
`SECURITY_REVIEW_MULTI_ENTRY.md`。

## 环境变量

| 变量 | 默认值 | 说明                                          |
| --- | --- |-----------------------------------------------|
| `DEM_PATH` | 仓库内置 `ldem_87s_5mpp_shackleton_50km_25m.tif` | DEM GeoTIFF 路径（默认随仓库分发，无需设置；可覆盖为其他 DEM） |
| `SPICE_KERNEL_DIR` | 仓库根 | SPICE 内核目录（含 `core/environment/data/`） |
| `PORT` / `VIEWER_PORT` | `5050` | 监听端口                                      |
| `HOST` | `0.0.0.0` | 监听地址                                      |
| `VIEWER_STATS_DIR` | `viewer/stats` | 访问统计 CSV 落盘目录（部署时指向持久盘）     |
| `FUSION_JOB_DIR` | `viewer/jobs` | 融合反解作业产物目录（部署时建议指向持久盘）  |
| `FUSION_TABLE_DIR` | `viewer/data/fusion_tables` | 融合反解 L5 家族表目录            |
| `FUSION_JOB_TIMEOUT_S` | `3600` | 融合反解作业超时（秒），超时终止子进程      |
| `FLASK_DEBUG` | 空 | 设为非空时开启模板热重载                      |

## 部署到 Render

DEM 使用仓库内置的压缩子集（52 MB，随仓库分发），**无需持久盘**，可直接部署。

方式一：使用 Blueprint（推荐）

1. 把本仓库推送到 GitHub；
2. Render Dashboard → New → Blueprint → 选择本仓库；
3. 自动按 `render.yaml` 创建 Web Service（内置 DEM 开箱即用）。

方式二：手动创建 Web Service

- Root Directory：`/`（仓库根）
- Environment：`Python 3`
- Build Command：`pip install -r requirements.txt`
- Start Command：`gunicorn viewer.server:app --bind 0.0.0.0:$PORT --workers 1 --threads 1 --timeout 300`

> 如需更高精度 DEM，可挂载 Disk 上传完整 DEM（3.2 GB）并设置
> `DEM_PATH=/opt/data/ldem_87s_5mpp.tif`；访问统计如需持久化，
> 设置 `VIEWER_STATS_DIR=/opt/data/stats` 并挂载 Disk。

健康检查：`GET /api/health` 返回 `{"status": "ok", "dem_exists": true, ...}`。

## 部署到 Vercel

DEM 使用仓库内置的压缩子集（52 MB，随仓库分发），**无需持久盘**，
所有端点（含 DEM 类：地形瓦片、地平线掩码、3D 地形）均可直接工作。
`vercel.json` 将全部路径转发至 Flask（`api/index.py`），三个功能区与静态
子应用均可用；轨道 JSON 为预生成资产，Serverless 内不做重计算。

```bash
vercel --prod
```

> **访问统计 CSV 无法持久化**（Serverless 实例文件系统为临时性），统计端点
> 会工作但数据会随实例回收丢失；如需持久统计，可接入外部存储
> （Vercel Blob / KV / Postgres）。
>
> **融合反向求解不可用**：Serverless 实例在响应后会被冻结/回收，无法承载
> 分钟级后台子进程；该功能需部署在常驻服务（Render / Docker / 自有服务器）。

## API 一览

| 端点 | 说明 |
| --- | --- |
| `GET /` | 站点首页（三入口导航） |
| `GET /visibility` | 可见性可视化前端页面 |
| `GET /orbits` | 轨道 3D 可视化前端页面 |
| `GET /coast_tool/` | 发射窗口计算器前端页面（静态子应用） |
| `GET /coast_tool/api/fusion/config` | 融合反解可用配置（支持的 T / 站点 / 规模） |
| `POST /coast_tool/api/fusion/jobs` | 启动作业（须 `confirm=true`；长耗时，返回 202 + job_id） |
| `GET /coast_tool/api/fusion/jobs/<id>` | 作业状态 / 结果（窗口列表 + 拒绝统计） |
| `DELETE /coast_tool/api/fusion/jobs/<id>` | 取消运行中的作业 |
| `GET /api/dem_tile?center_lat=&center_lon=&width_m=&height_m=&date=&hour=` | 光照 DEM 瓦片 PNG（附元数据响应头） |
| `GET /api/horizon_mask?lat=&lon=` | 360° 地形遮挡角 |
| `GET /api/earth_daily?lat=&lon=&date=` | 地球高度角/方位角日序列（288 点） |
| `GET /api/relay_daily?lat=&lon=&date=` | 鹊桥二号日序列（288 点） |
| `GET /api/earth_elevation` | 单时刻地球高度角 |
| `GET /api/relay_elevation` | 单时刻中继卫星高度角 |
| `GET /api/terrain_gltf?lat=&lon=` | 3D 地形 GLTF |
| `GET /api/stereo_to_lonlat?x=&y=` | 极射投影 → 经纬度 |
| `GET /api/yang_orbits?t_arrive=` | Yang 论文 24 条转移轨道 XYZ JSON（仅论文默认历元） |
| `GET /api/yang_orbits_stream?t_arrive=` | 同上，SSE 逐轨道推送（meta/orbit/done） |
| `POST /api/visit?path=` | 记录一次页面访问（CSV 持久化） |
| `POST /api/click?button=` | 记录一次按钮点击（CSV 持久化） |
| `GET /api/stats` | 按日期聚合的访问量 + 按钮点击量 |
| `GET /api/health` | 服务状态 |

## 科学口径

- 观测点默认 Shackleton 坑底（`-89.67°, 129.78°`，IAU/Gazetteer）。
- 鹊桥二号为月球大椭圆太阳同步测控回归冻结轨道（论文参数：a≈10004 km，
  e=0.782，i=119.2°），Ω=115°、M=215° 为假设值（证据等级 B），可见时长
  量级 10–20 h/天不受影响。
- 地平线掩码含月面曲率修正（30–50 km 量级不可忽略）。
- 地月转移轨道：杨维廉论文表 1–4 共 24 条（共面 12 / 非共面 4 / 同向 4 /
  极月 4），冻结 MEP 坐标系（到达时刻月地连线为 X 轴）；动力学为地心
  受限四体（地球 + 月球直接/间接项，月球取论文圆轨道模型），
  DOP853（rtol=1e-12），与论文数值回归一致。
- 发射窗口：停泊轨道滑行方案（LTO 面 = span{发射场@发射, 月球@到达}），
  经嫦娥五号/六号、鹊桥二号实测数据反演验证；四体引力修正支持任意目标
  环月倾角 i_M（含逆行）；算法、火箭能力口径与三个验证案例（含滑行偏长
  等已知误差）详见 `coast_tool/README.md`。

## License 说明与许可分界

本仓库采用**双许可分界**：代码与内容（文字/图表）适用不同的许可协议，
请按材料类型分别遵守。

| 材料 | 许可 | 文件 |
| --- | --- | --- |
| **代码**（Python / JS / HTML / CSS 等程序文件） | [MIT License](LICENSE) | `LICENSE` |
| **文字与图表**（文档、科学口径叙述、图表、渲染图、前端文案等） | [CC BY 4.0](figures-and-logs/LICENSE-CONTENT) | `LICENSE-CONTENT` |
| **第三方材料**（NASA 数据 / SPICE 内核 / 论文引用 / 开源库） | 归原始权利人所有，按各自许可使用 | [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) |

要点：

- **代码**：版权人 方辰，MIT 许可，可自由使用、修改、再分发，需保留版权声明。
- **内容**：版权人 方辰，CC BY 4.0 许可，可复制、修改、再分发（含商用），
  但必须署名并注明许可协议与改动情况。中文摘要见
  https://creativecommons.org/licenses/by/4.0/deed.zh-hans 。
- **第三方材料不在上述许可范围内**：
  - NASA 数据（月球 DEM、SPICE 内核）为美国公共领域作品，使用需致谢
    （详见 THIRD_PARTY_NOTICES.md §1）；
  - 科学口径引用的论文版权归原作者与期刊所有，本项目仅引用事实与数值；
  - 第三方开源库（Flask / three.js / spiceypy 等）按各自许可（MIT / BSD 等）
    使用，完整清单见 THIRD_PARTY_NOTICES.md §3。

移植组件来源（均已在站内自包含，运行期不依赖原仓库）：

- `launch_window/`：移植自作者自有项目 AgentMoons 的 `launch_window` 包
  （v0.6.0）：轨道链路（orbits_xyz / tli_solver / mep_frame / geometry）与
  融合窗口引擎（fusion / launch_epoch / perigee_match / scan_im / r1s_scan /
  launcher / lambert_window / lambert / month_accounting），内核路径与硬编码
  本地路径改指本仓库（`core/environment/data/`、`viewer/data/fusion_tables/`、
  `viewer/jobs/`），科学计算逻辑未改动；
- `coast_tool/`：移植自 lunarlaunchwindow 项目的 `coast_tool/`
  （MIT License, Copyright (c) 2026 Zhuoxiao；原始许可文本见
  `coast_tool/LICENSE`），内联脚本已外置为 `app.js` 以符合 CSP；其第一代
  （已证伪）工具与验证脚本未移植，验证结论保留在 `coast_tool/README.md`；
- 前端 vendored 的 three.js / GLTFLoader / OrbitControls（MIT,
  Copyright 2010-2022 Three.js Authors）许可文本见
  `viewer/static/vendor/LICENSE` 与 `coast_tool/vendor/LICENSE`；
- `viewer/data/yang_orbits_xyz.json`：由 `scripts/build_yang_orbits.py`
  调用 `launch_window` 离线生成，随仓库分发；
- `viewer/data/fusion_tables/`：由上游 `scan_im.run_scan` 离线生成的 L5
  家族表（T=112.0–113.0 h 与 112.61 h，i_m=90°、H_m=200 km、含太阳摄动），
  随仓库分发。

> 若你引用本项目科学参数，请同时引用原始论文（见 `analysis/elevation_video/queqiao_relay.py` 的 `REFERENCES`）。
