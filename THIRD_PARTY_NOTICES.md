# 第三方来源与许可声明（THIRD_PARTY_NOTICES.md）

本仓库（Queqiao2-visibility）在代码、数据与内容中使用了第三方来源的材料。
本文件逐一说明其来源与许可状态，作为 [LICENSE](LICENSE)（MIT，覆盖代码）与
[LICENSE-CONTENT](figures-and-logs/LICENSE-CONTENT)（CC BY 4.0，覆盖文字与图表）的补充。
凡本文件列出的第三方材料，其权利仍归原始权利人所有；本仓库仅按各自许可
条款使用并致谢，**不代表**原始权利人背书本项目。

---

## 1. NASA 数据（公共领域，需致谢）

美国政府的作品（含 NASA 及其任务产出）在美国不受版权保护，属于公共领域，
可自由使用；NASA 仅要求在使用时予以致谢。本项目使用并分发（或引用）以下
NASA 数据：

### 1.1 月球 DEM —— LDEM_87S_5MPP

| 项目 | 说明 |
| --- | --- |
| 数据 | 月球南极 87°S 以南 5 m/像素 DEM（GeoTIFF，`ldem_87s_5mpp.tif`） |
| 来源 | NASA GSFC / LROC（Lunar Reconnaissance Orbiter Camera）团队 |
| 任务 | Lunar Reconnaissance Orbiter（LRO） |
| 许可 | 公共领域（NASA 数据），需致谢 LROC 团队 |
| 致谢要求 | 使用时请注明：*"Lunar Reconnaissance Orbiter Camera (LROC) data, NASA/GSFC/Arizona State University"* |
| 备注 | 该 DEM 体积约 3.2 GB，**不随仓库分发**；通过环境变量 `DEM_PATH` 指向外部存储（见 README） |

### 1.2 SPICE 内核（随仓库分发于 `core/environment/data/`）

| 文件 | 内容 | 来源 |
| --- | --- | --- |
| `de442s.bsp` | 行星/卫星星历（DE442 精简版） | NASA JPL NAIF |
| `latest_leapseconds.tls` | 闰秒表 | NASA JPL NAIF |
| `moon_pa_de421_1900-2050.bpc` | 月极轴指向（DE421） | NASA JPL NAIF |
| `pck00011_n0066.tpc` | 行星常数内核 | NASA JPL NAIF |

- 来源：NASA JPL **NAIF**（Navigation and Ancillary Information Facility），
  通过其 SPICE 内核分发渠道（https://naif.jpl.nasa.gov/naif/）获取。
- 许可：NASA 公共领域作品，可自由分发与使用。
- 致谢要求：NAIF 建议引用 *"NAIF SPICE kernels, NASA JPL"*。
- 备注：本项目使用 [spiceypy](https://pypi.org/project/spiceypy/)（Python
  封装，BSD-3-Clause）加载上述内核。

### 1.3 官方图片（月球/地球/星空等素材）

如本项目前端页面或 `figures-and-logs/` 中使用了 NASA 官方发布的图片（如月球表面
影像、地球照片等），该类图片同为 NASA 公共领域作品；使用时保留 NASA 原始
图片说明与致谢信息（通常随图附注 *"NASA"* / *"NASA GSFC"* / *"NASA/JPL"*
等字样）。

本项目明确分发的 NASA 图片素材：

| 文件 | 内容 | 来源 |
| --- | --- | --- |
| `coast_tool/assets/earth_day.jpg` | 地球日面贴图（Blue Marble） | NASA Earth Observatory（Blue Marble 计划，MODIS 数据） |
| `coast_tool/assets/earth_night.jpg` | 地球夜面贴图（城市灯光） | NASA Earth Observatory / NOAA NGDC（Black Marble 系列） |

上述贴图仅用于 coast_tool 3D 场景的地球昼夜渲染，致谢要求：*"NASA Earth
Observatory"*；其他行星/月球纹理为程序生成（Canvas 绘制），不涉及第三方素材。

---

## 2. 学术论文引用（版权归原作者与期刊所有）

本项目在科学口径（轨道参数、可见弧段等）中引用了以下公开论文的**科学事实
与数值**。引用论文不属于本仓库创作内容，版权归原作者与出版期刊所有；
本项目**未复制**论文中的图表、版式或大段原文，仅引用事实与数值并注明出处，
符合学术引用惯例。如您引用本项目的科学参数，请同时引用原始论文。

1. **孟占峰.** 月球大椭圆太阳同步测控回归冻结轨道设计[J]. 航空学报,
   2024, 45(18): 210-225.
   - 引用内容：鹊桥二号轨道设计（a≈10004 km、e≈0.782、i≈119.2°、冻结
     ω≈±90°、设计近月点高度 600 km 等）。
2. **Zhang Lihua et al.** 2024. *System design and validation of
   Queqiao-2.* doi:10.16708/j.cnki.1000-758X.2024.0071.
   - 引用内容：鹊桥二号系统设计、测控回归周期 T_R≈24.84 h。
3. **Zhou Wenyan et al.** 2024. *Matching design of trajectory with
   Chang'e-6 mission.* doi:10.16708/j.cnki.1000-758X.2024.0087.
   - 引用内容：鹊桥二号与嫦娥六号任务的轨道匹配设计、最长中继弧段
     ≈20.7 h/天（论文表 4–5）。
4. **杨维廉.** 发射极月卫星的转移轨道研究（北京空间飞行器总体设计部，
   收稿 1997-10-21）.
   - 引用内容：MEP 月地平面坐标系定义、R-1 月球侧边值问题（q/p 变量与
     状态转移矩阵迭代）、表 1–4 的转移轨道数值（近地点 200 km、近月点
     300 km，T=72/96/120 h 等）。本项目移植的 `launch_window` 引擎即
     实现该论文的数值求解框架，`viewer/data/yang_orbits_xyz.json` 为
     按论文数值回归后的轨道资产。
5. **Bill Gray, Find_Orb（projectpluto.com）.** 嫦娥五号/六号、鹊桥二号
   光学定轨密切根数（TLE）.
   - 引用内容：`coast_tool/README.md` 中三次历史任务实测倾角的比对来源；
     原始 TLE 文件未随本仓库分发（仅引用数值事实）。

> 备注：完整引用块见
> `analysis/elevation_video/queqiao_relay.py` 中的 `REFERENCES`。

---

## 3. 开源软件（第三方库 / 前端库）

| 组件 | 版本/位置 | 许可 |
| --- | --- | --- |
| three.js | `viewer/static/vendor/three.min.js` | MIT |
| GLTFLoader | `viewer/static/vendor/GLTFLoader.js`（three.js 示例模块） | MIT |
| three.js（coast_tool 副本） | `coast_tool/vendor/three.min.js`（r147） | MIT |
| OrbitControls | `coast_tool/vendor/OrbitControls.js`（three.js 示例模块） | MIT |
| Flask | PyPI | BSD-3-Clause |
| NumPy | PyPI | BSD-3-Clause |
| rasterio | PyPI | BSD-3-Clause |
| matplotlib | PyPI | PSF-based license |
| Pillow | PyPI | HPND（Historical Permission Notice and Disclaimer） |
| SciPy | PyPI | BSD-3-Clause |
| spiceypy | PyPI | BSD-3-Clause |
| pyvista | PyPI | MIT |

完整依赖清单见 `requirements.txt`。上述库的完整许可文本以各库随附
LICENSE 文件为准。前端 vendored 的 three.js / GLTFLoader / OrbitControls
版权归 Three.js 贡献者所有（MIT, Copyright 2010-2022 Three.js Authors），
版权声明保留于 `three.min.js` 文件头，完整许可文本随文件分发：

- `viewer/static/vendor/LICENSE`（three.min.js / GLTFLoader.js）
- `coast_tool/vendor/LICENSE`（three.min.js / OrbitControls.js）

上游仓库：https://github.com/mrdoob/three.js

---

## 4. 移植组件来源（本项目作者自有项目）

以下组件由本项目作者从自有项目移植并已完全自包含（运行期不引用原项目
路径、数据或服务），在本仓库内继续适用 MIT 许可（`LICENSE`）：

| 组件 | 来源 | 说明 |
| --- | --- | --- |
| `launch_window/` | AgentMoons 项目 `launch_window` 包（v0.6.0） | 仅保留 Yang 论文轨道求解子集（orbits_xyz / tli_solver / mep_frame / geometry）；SPICE 内核路径改指本仓库 `core/environment/data/` |
| `coast_tool/` | lunarlaunchwindow 项目 `coast_tool/`（MIT, Copyright (c) 2026 Zhuoxiao） | 内联脚本已外置为 `app.js` 以符合 CSP；原始 MIT 许可文本见 `coast_tool/LICENSE`；`启动.command` / `node_modules` 等未移植 |
| `viewer/data/yang_orbits_xyz.json` | 由本仓库 `scripts/build_yang_orbits.py` 调用 `launch_window` 离线生成 | 论文默认到达历元（2026-08-24T00:00:00Z）的 24 条轨道 |

---

## 5. 汇总：本项目材料与第三方材料的边界

| 材料 | 归属 | 许可 |
| --- | --- | --- |
| 源代码（Python/JS/HTML/CSS，含 vendored 前端库中的库本体除外） | 本项目作者 方辰 | MIT（`LICENSE`） |
| 移植代码（`launch_window/`、`coast_tool/`，见 §4） | 本项目作者（自有项目） | MIT（`LICENSE`） |
| 文档、科学口径文字、图表、渲染图 | 本项目作者 方辰 | CC BY 4.0（`LICENSE-CONTENT`） |
| SPICE 内核（`core/environment/data/`） | NASA JPL NAIF | 公共领域（致谢 NAIF） |
| 月球 DEM（`ldem_87s_5mpp.tif`，外部存储） | NASA GSFC / LROC | 公共领域（致谢 LROC） |
| 论文中的事实与数值（本项目仅引用） | 原作者与出版期刊 | 版权归原权利人，学术引用 |
| 第三方 Python 库 / three.js | 各自作者 | 见上表及各自 LICENSE |

---

*本文件最后更新：2026.9.21。若发现遗漏或错误，欢迎提交 Issue / PR 修正。*
