"""launch_window — 杨维廉论文轨道 R-1 求解与 XYZ 轨迹采样引擎（移植子集）。

本包移植自作者自有的 AgentMoons 项目 ``launch_window`` 包（v0.6.0），
仅保留 Yang 论文轨道可视化所需的求解器链路::

    orbits_xyz.py  — 论文轨道（表 1–4，共 24 条）构建为 XYZ 记录 / JSON
    tli_solver.py  — R-1 月球侧边值问题求解（DOP853 + Newton）与轨迹采样
    mep_frame.py   — MEP（月地近焦点）冻结坐标系
    geometry.py    — SPICE 几何与内核加载

依赖仓库内置的 SPICE 内核（``core/environment/data/``），无外部数据依赖。
SPICE 非线程安全；部署须单线程/单 worker（与 viewer/server.py 一致）。
"""

__version__ = "0.6.0"
