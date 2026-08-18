# Queqiao2-visibility — 月球南极地平线掩码交互可视化

基于 Flask + Canvas 2D / three.js 的交互式可视化工具：左侧为南极极射投影
DEM 光照地形图（实时太阳光照渲染），右侧为极坐标地平线图 / 第一人称 3D
天空视角，底部时间滑块联动展示地球与鹊桥二号（Queqiao-2）中继卫星的全天
可见性窗口。

数据由服务端实时计算（SPICE 天历 + 月球 DEM），前端全部基于 Canvas/WebGL
客户端渲染。DEM 使用仓库内置的**压缩子集** `ldem_87s_5mpp_shackleton_50km_25m.tif`
（沙克尔顿坑 50 km 半径、25 m/px、52 MB，随仓库分发，**无需持久盘**）；
如需更高精度可自行用完整 DEM（3.2 GB）覆盖 `DEM_PATH`。

## 目录结构

```
viewer/                     # Flask 应用
  server.py                 # 后端（生产 WSGI 入口：viewer.server:app）
  templates/index.html      # 前端页面
  static/app.js             # 前端逻辑
  static/vendor/            # three.js / GLTFLoader（本地 vendored）
  stats/                    # 访问统计 CSV（运行时生成，git-ignored）
analysis/elevation_video/   # 科学计算模块（地平线掩码 / 轨道传播 / ENU 投影 / 中文字体）
core/environment/data/      # SPICE 内核（随仓库分发）
api/index.py                # Vercel Serverless 入口
render.yaml                 # Render Blueprint
vercel.json                 # Vercel 配置
Procfile                    # Render / 通用平台进程定义
requirements.txt            # Python 依赖
```

功能亮点：
- 左侧月球光照地形图（实时太阳光照渲染）、右侧极坐标地平线图 / 第一人称 3D 天空视角；
- 底部时间滑块联动地球与鹊桥二号全天可见性窗口；
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
#    http://127.0.0.1:5050
```

生产 WSGI（gunicorn / waitress）：

```bash
gunicorn viewer.server:app --bind 0.0.0.0:$PORT --workers 1 --threads 1 --timeout 300
# 或
waitress-serve --host 0.0.0.0 --port $PORT viewer.server:app
```

> 单线程（`--workers 1 --threads 1`）：SPICE（spiceypy）非线程安全。

## 环境变量

| 变量 | 默认值 | 说明                                          |
| --- | --- |-----------------------------------------------|
| `DEM_PATH` | 仓库内置 `ldem_87s_5mpp_shackleton_50km_25m.tif` | DEM GeoTIFF 路径（默认随仓库分发，无需设置；可覆盖为其他 DEM） |
| `SPICE_KERNEL_DIR` | 仓库根 | SPICE 内核目录（含 `core/environment/data/`） |
| `PORT` / `VIEWER_PORT` | `5000` | 监听端口                                      |
| `HOST` | `0.0.0.0` | 监听地址                                      |
| `VIEWER_STATS_DIR` | `viewer/stats` | 访问统计 CSV 落盘目录（部署时指向持久盘）     |
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

```bash
vercel --prod
```

> **访问统计 CSV 无法持久化**（Serverless 实例文件系统为临时性），统计端点
> 会工作但数据会随实例回收丢失；如需持久统计，可接入外部存储
> （Vercel Blob / KV / Postgres）。

## API 一览

| 端点 | 说明 |
| --- | --- |
| `GET /` | 前端页面 |
| `GET /api/dem_tile?center_lat=&center_lon=&width_m=&height_m=&date=&hour=` | 光照 DEM 瓦片 PNG（附元数据响应头） |
| `GET /api/horizon_mask?lat=&lon=` | 360° 地形遮挡角 |
| `GET /api/earth_daily?lat=&lon=&date=` | 地球高度角/方位角日序列（288 点） |
| `GET /api/relay_daily?lat=&lon=&date=` | 鹊桥二号日序列（288 点） |
| `GET /api/earth_elevation` | 单时刻地球高度角 |
| `GET /api/relay_elevation` | 单时刻中继卫星高度角 |
| `GET /api/terrain_gltf?lat=&lon=` | 3D 地形 GLTF |
| `GET /api/stereo_to_lonlat?x=&y=` | 极射投影 → 经纬度 |
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

> 若你引用本项目科学参数，请同时引用原始论文（见 `analysis/elevation_video/queqiao_relay.py` 的 `REFERENCES`）。
