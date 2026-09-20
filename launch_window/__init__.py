"""launch_window — 杨维廉论文轨道 R-1 求解与融合发射窗口引擎（移植子集）。

本包移植自作者自有的 AgentMoons 项目 ``launch_window`` 包（v0.6.0），
包含两部分：

1. **Yang 论文轨道可视化链路**（``/orbits`` 页面离线资产）::

       orbits_xyz.py  — 论文轨道（表 1–4，共 24 条）构建为 XYZ 记录 / JSON
       tli_solver.py  — R-1 月球侧边值问题求解（DOP853 + Newton）与轨迹采样
       mep_frame.py   — MEP（月地近焦点）冻结坐标系
       geometry.py    — SPICE 几何与内核加载

2. **融合发射窗口引擎**（``/coast_tool`` 页面的可选反向求解）::

       fusion.py           — 前向/反向融合框架（FR-1..FR-8）
       launch_epoch.py     — L6 反解链路（家族表 + 口袋搜索）
       perigee_match.py    — 近地点匹配家族表
       scan_im.py          — 任意 i_m 家族表构建（预生成资产见
                             viewer/data/fusion_tables/）
       r1s_scan.py         — R-1S 家族连续扫描
       launcher.py         — C3/能力曲线与窗口宽度
       lambert_window.py   — 发射场容器与 Lambert 窗口几何
       lambert.py          — Lambert 初值/经典根数（poliastro 惰性导入，
                             反解链路不需要）
       month_accounting.py — 月份账目工具

依赖仓库内置的 SPICE 内核（``core/environment/data/``），无外部数据依赖。
SPICE 非线程安全：Web 部署须单线程/单 worker（与 viewer/server.py 一致）；
融合反解为 CPU 密集长作业，由 Web 层以**子进程**方式调度，避免阻塞服务。

移植适配：仅内核路径、硬编码本地绝对路径（改为环境变量可覆盖的仓库内
默认路径）与 CLI 默认输出路径有改动，科学计算逻辑未改动。
"""

__version__ = "0.6.0"
