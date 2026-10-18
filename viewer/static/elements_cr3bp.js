/*!
 * elements_cr3bp.js — 月球轨道与地月圆形限制性三体（CR3BP）模块。
 *
 * 与 elements_core.js / elements_perturb.js 同一 UMD 双端约定
 * （浏览器挂 window.ElementsCr3bp，node 走 module.exports），
 * 加载时自动向 ElementsCore 注册三个传播器：
 *   kepler-moon    开普勒二体（月球）
 *   perturbed-moon 月心 Cowell 数值摄动（月球 J2–J4 + 地球第三体 + 太阳光压）
 *   cr3bp          地月会合系三体传播器（Halo / NRHO / DRO 族周期轨道，
 *                  及根数初值模式 elements：由开普勒根数直接传播，如 ELFO）
 *
 * ── 力学与坐标系（详细推导与出处见 tools/interactive_viewer/月球三体-算法与验证.md）──
 *
 * 月球选项卡显示系（本模块所有 km 级坐标的参考系）：
 *   - 原点 = 月心；ẑ = 月球北极轴（MOON_PA，SPICE 历元 2026-10-13 锁定）；
 *   - x̂ = 月赤道面与地月轨道面的交线（极轴×轨道法向）；
 *   - 该系在 90 d 窗口内视为惯性系（ Cassini 平面 18.6 a 进动、极轴自由章动忽略）。
 *   - 月球带谐引力场 (J2–J4) 关于 ẑ 精确对称 —— 冻结轨道条件所依赖的
 *     倾角即相对月赤道之倾角，与 B10/S13 文献口径一致。
 *
 * 地球圆化轨道（月心第三体模型）：
 *   - 半径 a_EM = 384400 km（CR3BP 特征长度），角速度 ω_EM = √(G(m⊕+m☾)/a³)；
 *   - 轨道面 = 地月轨道面，对月赤道倾斜 Cassini 角 6.687°（SPICE 实测）；
 *   - 历元方位由 SPICE de442s 锁定（spice_probe.py 输出）。
 *   - 第三体摄动 = 直接项 + 间接项（受摄二体相对运动方程标准形式）。
 *
 * CR3BP（会合坐标系，无量纲）：
 *   - 特征长度 l = a_EM、特征质量 m⊕+m☾、特征时间 1/ω_EM；
 *   - μ = m☾/(m⊕+m☾)；P1(地球) 位于 (−μ,0,0)，P2(月球) 位于 (1−μ,0,0)；
 *   - 运动方程 ẍ−2ẏ = ∂Ω/∂x 等，Ω = (1−μ)/r₁ + μ/r₂ + (x²+y²)/2；
 *   - 雅可比积分 C = 2Ω − v²（旋转系速度）；
 *   - 共线平动点 L1/L2/L3 由 ∂Ω/∂x = 0 的五次方程二分求解。
 *
 * 周期轨道初值（Halo / NRHO / DRO）：
 *   - 微分修正（单打靶）：Halo 族固定 z₀（族参数 Az），自由变元 (x₀, v_y0, T/2)，
 *     约束半周期后回到 xz 平面垂直穿越（y=vx=vz=0）；
 *   - DRO 平面族固定周期 T（2:1 共振 = 恒星月/2），自由变元 (x₀, v_y0)，
 *     约束 y(T/2)=0, vx(T/2)=0；
 *   - 预设初值由 tools/interactive_viewer/tune_cr3bp_presets.cjs 离线解出后固化，
 *     页面内 Az 滑条调整时从当前解热启动续修正（1–3 次牛顿迭代收敛）。
 *
 * 单位约定：km / s / km³·s⁻²；CR3BP 内部无量纲，输出前统一换算。
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) {
    module.exports = factory(require('./elements_core.js'));
  } else {
    root.ElementsCr3bp = factory(root.ElementsCore);
  }
}(typeof globalThis !== 'undefined' ? globalThis : this, function (C) {
  'use strict';

  if (!C) { throw new Error('ElementsCr3bp 需要 ElementsCore 先加载'); }

  /* ══ 常数（全部注明出处；换算过程见算法文档） ══════════════════════ */

  var MU_MOON = 4902.800066;        // km³/s²，DE440 月心引力常数
  var R_MOON = 1737.4;              // km，月球平均半径（高度/撞月判断）
  var R_MOON_GRAV = 1738.0;         // km，GL0660B 参考半径（带谐项专用）

  /* 月球带谐系数（GL0660B，JPL GRAIL 主任务解，未归一化 J 值）：
     归一化系数（ICGEM gl0660b.gfc）C̄20 = −9.0880834662e-5、
     C̄30 = −3.1974673493e-6、C̄40 = +3.2348047451e-6；
     未归一化换算 C_n0 = √(2n+1)·C̄_n0，J_n = −C_n0。 */
  var J2_MOON = Math.sqrt(5) * 9.088083466222001e-5;   // ≈ 2.03216e-4
  var J3_MOON = -Math.sqrt(7) * (-3.1974673492600000e-6); // ≈ +8.4597e-6
  var J4_MOON = -Math.sqrt(9) * 3.2348047451430000e-6;    // ≈ −9.7044e-6

  var MU_EARTH = 398600.4418;       // km³/s²，与 elements_perturb 一致（页面基准）
  var MU_EM_TOTAL = MU_EARTH + MU_MOON;
  var MU_RATIO = MU_MOON / MU_EM_TOTAL;   // ≈ 0.0121506
  var A_EM = 384400;                // km，地月特征距离（圆化模型半径）
  var OMEGA_EM = Math.sqrt(MU_EM_TOTAL / (A_EM * A_EM * A_EM)); // rad/s
  var T_EM = 2 * Math.PI / OMEGA_EM;      // s，模型恒星月
  var YEAR_S = 365.25 * 86400;
  var P_SUN = 4.56e-6;              // N/m²，与 elements_perturb 一致

  var WINDOW_S = 90 * 86400;        // 传播窗口：90 天

  /* ══ SPICE 历元锚定（spice_probe.py @ 2026-10-13 00:00 UTC，J2000） ══
   * kernel：data/lunar/ephemeris/{de442s.bsp, moon_pa_de421_1900-2050.bpc, ...}
   * 显示系基（J2000 表达，仅存档注释用；运行时一律用显示系内参数）：
   *   ẑ(极轴) = [0.015586, -0.376493, 0.926288]
   *   x̂(交线) = [0.844452, -0.491103, -0.213819]
   *   ŷ       = [0.535405,  0.785539,  0.310276]
   */
  var CASSINI_DEG = 6.687;          // 极轴 vs 轨道面法向夹角（SPICE 实测）
  var THETA_E0 = 1.382572;          // rad，历元地球方位（自 x̂、沿轨道面）
  /* 轨道面第二基矢 v2 = 轨道法向×x̂ 在显示系中的表达（单位矢） */
  var V2_DISP = [0, 0.993196, 0.116451];
  /* 会合系基（e1=地球→月球, e2, e3=轨道法向）在显示系中的表达（行主序：第 i 行 = e_i
     的显示系分量，三行两两正交）。轨道面惯性系 → 显示系 d = Σ_i o_i·e_i_disp = Mᵀ·o，
     须用转置乘（见 matTVec3 / orbitToDisp）。 */
  var M_ORBIT2DISP = [
    -0.187115, -0.975655, -0.114394,   // e1_disp
     0.982338, -0.185842, -0.021790,   // e2_disp
     0.000000, -0.116451,  0.993196,   // e3_disp
  ];
  var LAM0_SUN = -2.239594;         // rad，历元太阳方位（自 x̂、沿月赤道面）
  /* 太阳离月赤道面仅 ±1.5°（历元 −1.28°），光压模型取面内圆运动（忽略项） */

  /* ══ 基础矢量工具 ═════════════════════════════════════════════════ */

  function dot3(a, b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
  function cross3(a, b) {
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  }
  function matVec3(m, v) {
    return [
      m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
      m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
      m[6] * v[0] + m[7] * v[1] + m[8] * v[2]
    ];
  }
  /** 转置乘 Mᵀ·v（M 按行存放基矢时，Mᵀ·o 才是「轨道面惯性系 → 显示系」）。 */
  function matTVec3(m, v) {
    return [
      m[0] * v[0] + m[3] * v[1] + m[6] * v[2],
      m[1] * v[0] + m[4] * v[1] + m[7] * v[2],
      m[2] * v[0] + m[5] * v[1] + m[8] * v[2]
    ];
  }
  function clamp1(x) { return Math.min(1, Math.max(-1, x)); }

  /* ══ 月心环境模型（显示系） ═══════════════════════════════════════ */

  /** 地月轨道面惯性系（基 e1,e2,e3）→ 显示系：v_disp = Σ_i v_i·e_i_disp = Mᵀ·v_orbit。
      曾误用 matVec3（= M·v，实为「显示系 → 轨道面系」），使会合系轨迹 / 平动点 /
      航天器整体相对真实地球方位整体旋转 ~157°：L2 与地球方位夹角仅 22.6°（应 ≈180°），
      也与 earthPosMoonTab 给出的地球位置不符（夹角 157°）。 */
  function orbitToDisp(v) { return matTVec3(M_ORBIT2DISP, v); }

  /** 地球月心位置（显示系 km），t 自历元起算（s）。圆化模型（半径 a_EM）。 */
  function earthPosMoonTab(tSec) {
    var th = THETA_E0 + OMEGA_EM * tSec;
    var c = Math.cos(th), s = Math.sin(th);
    return [A_EM * c, A_EM * s * V2_DISP[1], A_EM * s * V2_DISP[2]];
  }

  /** 太阳方向单位矢（显示系），面内圆运动假设（黄道≈月赤道面，±1.5° 忽略）。 */
  function sunDirMoonTab(tSec) {
    var lam = LAM0_SUN + 2 * Math.PI * tSec / YEAR_S;
    return [Math.cos(lam), Math.sin(lam), 0];
  }

  /** 圆柱形月影受晒系数：本影内 0，外部 1（无半影；地球影忽略——月距处地影
      截面 ~4600 km 且仅当卫星-地球-太阳近似共线才涉及，对月轨道演示可忽略）。 */
  function moonShadowFactor(r, sHat) {
    var rdots = dot3(r, sHat);
    if (rdots >= 0) { return 1; }   // 向阳半球侧
    var r2 = dot3(r, r);
    var d2 = r2 - rdots * rdots;    // 到月-日轴垂距²
    return d2 < R_MOON * R_MOON ? 0 : 1;
  }

  /**
   * 月球带谐全引力（中心 + J2–J4，显示系 km/s²）。
   * 位函数 U = μ☾/r·[1 − Σ Jₙ(R_ref/r)ⁿ Pₙ(z/r)]，梯度与地球版同型
   * （elements_perturb.zonalAccel 换常数；R_ref 用 GL0660B 参考半径 1738.0 km）。
   */
  function moonZonalAccel(r, out) {
    var x = r[0], y = r[1], z = r[2];
    var r2 = x * x + y * y + z * z;
    var rm = Math.sqrt(r2);
    var s = z / rm;
    var P2 = (3 * s * s - 1) / 2;
    var P3 = (5 * s * s * s - 3 * s) / 2;
    var P4 = (35 * s * s * s * s - 30 * s * s + 3) / 8;
    var P2d = 3 * s;
    var P3d = (15 * s * s - 3) / 2;
    var P4d = (35 * s * s * s - 30 * s) / 4;
    var Re_r = R_MOON_GRAV / rm;
    var t2 = Re_r * Re_r, t3 = t2 * Re_r, t4 = t3 * Re_r;
    var sumP = J2_MOON * t2 * P2 + J3_MOON * t3 * P3 + J4_MOON * t4 * P4;
    var sumR = 2 * J2_MOON * t2 * P2 + 3 * J3_MOON * t3 * P3 + 4 * J4_MOON * t4 * P4;
    var sumS = J2_MOON * t2 * P2d + J3_MOON * t3 * P3d + J4_MOON * t4 * P4d;
    var mu_r = MU_MOON / rm;
    var dUdr = -mu_r / rm * (1 - sumP - sumR);
    var dUds = -mu_r * sumS;
    var r3 = r2 * rm;
    out[0] = dUdr * (x / rm) + dUds * (-x * z / r3);
    out[1] = dUdr * (y / rm) + dUds * (-y * z / r3);
    out[2] = dUdr * (z / rm) + dUds * ((x * x + y * y) / r3);
    return out;
  }

  /** 月心笛卡尔 → 经典根数（角度弧度；mu 参数化，供月球/地球复用）。 */
  function cartToClassicalMu(r, v, mu) {
    var rm = Math.hypot(r[0], r[1], r[2]);
    var vm2 = dot3(v, v);
    var h = cross3(r, v);
    var hm = Math.hypot(h[0], h[1], h[2]);
    var En = vm2 / 2 - mu / rm;
    var a = -mu / (2 * En);
    var i = Math.acos(clamp1(h[2] / hm));
    var n = [-h[1], h[0], 0];
    var nm = Math.hypot(n[0], n[1]);
    var Om = nm > 1e-12 ? Math.atan2(n[1], n[0]) : 0;
    var vh = cross3(v, h);
    var eVec = [vh[0] / mu - r[0] / rm, vh[1] / mu - r[1] / rm, vh[2] / mu - r[2] / rm];
    var e = Math.hypot(eVec[0], eVec[1], eVec[2]);
    var hunit = [h[0] / hm, h[1] / hm, h[2] / hm];
    function angAbout(u, w) {
      var cuw = cross3(u, w);
      return Math.atan2(dot3(cuw, hunit), dot3(u, w));
    }
    var w, nu;
    if (e > 1e-9) {
      w = nm > 1e-12 ? angAbout(n, eVec) : Math.atan2(eVec[1], eVec[0]);
      nu = angAbout(eVec, r);
    } else {
      w = 0;
      nu = nm > 1e-12 ? angAbout(n, r) : Math.atan2(r[1], r[0]);
    }
    return { a: a, e: e, i: i, Om: C.wrap2pi(Om), w: C.wrap2pi(w), nu: C.wrap2pi(nu) };
  }

  /* ══ CR3BP 核心（无量纲会合系） ═══════════════════════════════════ */

  var MU = MU_RATIO;

  /** 有效势 Ω(x,y,z)。 */
  function omegaPot(x, y, z) {
    var r1 = Math.sqrt((x + MU) * (x + MU) + y * y + z * z);
    var r2 = Math.sqrt((x - 1 + MU) * (x - 1 + MU) + y * y + z * z);
    return (1 - MU) / r1 + MU / r2 + 0.5 * (x * x + y * y);
  }

  /** 会合系运动方程右端（无量纲）。y = [x,y,z,vx,vy,vz]。 */
  function cr3bpDeriv(y, out) {
    var x = y[0], yy = y[1], z = y[2];
    var r1sq = (x + MU) * (x + MU) + yy * yy + z * z;
    var r2sq = (x - 1 + MU) * (x - 1 + MU) + yy * yy + z * z;
    var r1 = Math.sqrt(r1sq), r2 = Math.sqrt(r2sq);
    var r13 = r1sq * r1, r23 = r2sq * r2;
    var k1 = (1 - MU) / r13, k2 = MU / r23;
    out[0] = y[3]; out[1] = y[4]; out[2] = y[5];
    out[3] = 2 * y[4] + x - k1 * (x + MU) - k2 * (x - 1 + MU);
    out[4] = -2 * y[3] + yy - k1 * yy - k2 * yy;
    out[5] = -k1 * z - k2 * z;
    return out;
  }

  /** RK4 积分会合系方程：从 y0 积 nSteps 步（步长 h），返回末端状态。 */
  function integrateRot(y0, h, nSteps) {
    var y = Float64Array.from(y0);
    var k1 = new Float64Array(6), k2 = new Float64Array(6),
        k3 = new Float64Array(6), k4 = new Float64Array(6), yt = new Float64Array(6);
    for (var i = 0; i < nSteps; i++) {
      cr3bpDeriv(y, k1);
      for (var j = 0; j < 6; j++) { yt[j] = y[j] + 0.5 * h * k1[j]; }
      cr3bpDeriv(yt, k2);
      for (var j2 = 0; j2 < 6; j2++) { yt[j2] = y[j2] + 0.5 * h * k2[j2]; }
      cr3bpDeriv(yt, k3);
      for (var j3 = 0; j3 < 6; j3++) { yt[j3] = y[j3] + h * k3[j3]; }
      cr3bpDeriv(yt, k4);
      for (var j4 = 0; j4 < 6; j4++) {
        y[j4] += h / 6 * (k1[j4] + 2 * k2[j4] + 2 * k3[j4] + k4[j4]);
      }
    }
    return Array.from(y);
  }

  /** 雅可比积分 C = 2Ω − v²（无量纲）。 */
  function jacobi(y) {
    var v2 = y[3] * y[3] + y[4] * y[4] + y[5] * y[5];
    return 2 * omegaPot(y[0], y[1], y[2]) - v2;
  }

  /** CR3BP 状态方程雅可比 A = ∂f/∂y（中心差分，行主序 A[i·6+j] = ∂f_i/∂y_j）。 */
  function cr3bpJacobian(y, A) {
    var eps = 1e-7;
    var fp = new Float64Array(6), fm = new Float64Array(6);
    var yp = Float64Array.from(y), ym = Float64Array.from(y);
    for (var j = 0; j < 6; j++) {
      yp[j] += eps; ym[j] -= eps;
      cr3bpDeriv(yp, fp); cr3bpDeriv(ym, fm);
      for (var i = 0; i < 6; i++) { A[i * 6 + j] = (fp[i] - fm[i]) / (2 * eps); }
      yp[j] -= eps; ym[j] += eps;
    }
  }

  /** 单步 RK4 推进会合系状态 y（原地修改）。 */
  function rk4Step(y, h) {
    var k1 = new Float64Array(6), k2 = new Float64Array(6),
        k3 = new Float64Array(6), k4 = new Float64Array(6), yt = new Float64Array(6);
    cr3bpDeriv(y, k1);
    for (var j = 0; j < 6; j++) { yt[j] = y[j] + 0.5 * h * k1[j]; }
    cr3bpDeriv(yt, k2);
    for (var j2 = 0; j2 < 6; j2++) { yt[j2] = y[j2] + 0.5 * h * k2[j2]; }
    cr3bpDeriv(yt, k3);
    for (var j3 = 0; j3 < 6; j3++) { yt[j3] = y[j3] + h * k3[j3]; }
    cr3bpDeriv(yt, k4);
    for (var j4 = 0; j4 < 6; j4++) { y[j4] += h / 6 * (k1[j4] + 2 * k2[j4] + 2 * k3[j4] + k4[j4]); }
  }

  /** 单步 RK4 推进耦合系统 [状态(6) | STM(36)]（Φ' = A(y)·Φ，原地修改）。 */
  function rk4StepStm(s, h) {
    var N = 42;
    var k1 = new Float64Array(N), k2 = new Float64Array(N),
        k3 = new Float64Array(N), k4 = new Float64Array(N), st = new Float64Array(N);
    var A = new Float64Array(36);
    var f = function (src, out) {
      cr3bpDeriv(src, out);
      cr3bpJacobian(src, A);
      for (var i = 0; i < 6; i++) {
        for (var j = 0; j < 6; j++) {
          var sum = 0;
          for (var m = 0; m < 6; m++) { sum += A[i * 6 + m] * src[6 + m * 6 + j]; }
          out[6 + i * 6 + j] = sum;
        }
      }
    };
    f(s, k1);
    for (var j = 0; j < N; j++) { st[j] = s[j] + 0.5 * h * k1[j]; }
    f(st, k2);
    for (var j2 = 0; j2 < N; j2++) { st[j2] = s[j2] + 0.5 * h * k2[j2]; }
    f(st, k3);
    for (var j3 = 0; j3 < N; j3++) { st[j3] = s[j3] + h * k3[j3]; }
    f(st, k4);
    for (var j4 = 0; j4 < N; j4++) { s[j4] += h / 6 * (k1[j4] + 2 * k2[j4] + 2 * k3[j4] + k4[j4]); }
  }

  /** 共线平动点 x 坐标（无量纲、质心系）：解 ∂Ω/∂x=0（五次方程）。 */
  function collinearX(which) {
    function g(x) {
      var r1 = Math.abs(x + MU), r2 = Math.abs(x - 1 + MU);
      return x - (1 - MU) * (x + MU) / (r1 * r1 * r1) - MU * (x - 1 + MU) / (r2 * r2 * r2);
    }
    /* g 单调区间内二分：L1 ∈ (月球左, 月球)、L2 ∈ (月球, 月球右)、L3 ∈ (地球左, 地球)
       （括号端点取在奇点两侧，g 在奇点附近趋于 ±∞，符号必异） */
    var lo, hi;
    if (which === 'L1') { lo = 1 - MU - 0.3; hi = 1 - MU - 1e-9; }
    else if (which === 'L2') { lo = 1 - MU + 1e-9; hi = 1 - MU + 0.3; }
    else { lo = -1.5; hi = -MU - 1e-9; }   // L3 在地球左侧 ≈ −1.005
    var glo = g(lo), ghi = g(hi);
    if (glo * ghi > 0) { return NaN; }
    for (var k = 0; k < 200; k++) {
      var mid = 0.5 * (lo + hi);
      var gm = g(mid);
      if (gm * glo <= 0) { hi = mid; } else { lo = mid; glo = gm; }
    }
    return 0.5 * (lo + hi);
  }
  var X_L1 = collinearX('L1');
  var X_L2 = collinearX('L2');
  var X_L3 = collinearX('L3');
  var C_L1 = 2 * omegaPot(X_L1, 0, 0);   // 平动点处 v=0，C = 2Ω
  var C_L2 = 2 * omegaPot(X_L2, 0, 0);
  var C_L3 = 2 * omegaPot(X_L3, 0, 0);

  /* ── 坐标换算：会合系(无量纲) ↔ 显示系(km) ─────────────────────────
   * 轨道面惯性系（基 e1,e2,e3）：θ = ω_EM·t
   *   r_bary = l·[x cosθ − y sinθ, x sinθ + y cosθ, z]
   *   r_moon = (1−μ)·l·[cosθ, sinθ, 0]
   *   r_relMoon = r_bary − r_moon（轨道面系）
   *   r_disp = M·r_relMoon（显示系）
   * 速度对 t 求导后再经同一旋转。 */

  function rotStateToDisplay(y, tSec) {
    var th = OMEGA_EM * tSec;
    var c = Math.cos(th), s = Math.sin(th);
    var x = y[0], yy = y[1], z = y[2];
    var vx = y[3], vy = y[4], vz = y[5];
    var l = A_EM, mu = MU;
    /* 轨道面系位置（月心相对） */
    var ro = [
      l * (x * c - yy * s) - (1 - mu) * l * c,
      l * (x * s + yy * c) - (1 - mu) * l * s,
      l * z
    ];
    /* 轨道面系速度：d/dt[r_bary] − d/dt[r_moon]。
       v_bary = l·ω·R′(θ)·ρ + l·ω·R(θ)·ρ̇（ρ̇ 对无量纲时间 τ=ωt 求导，
       故 R·ρ̇ 项必须带 OMEGA_EM 因子——曾漏乘导致显示速度虚高 ~1/ω 倍，
       Halo/DRO 速度卡显示数万千米每秒）。 */
    var vo = [
      OMEGA_EM * l * (-(x * s + yy * c)) + OMEGA_EM * l * (vx * c - vy * s) + (1 - mu) * OMEGA_EM * l * s,
      OMEGA_EM * l * (x * c - yy * s) + OMEGA_EM * l * (vx * s + vy * c) - (1 - mu) * OMEGA_EM * l * c,
      l * OMEGA_EM * vz
    ];
    return { rE: orbitToDisp(ro), vE: orbitToDisp(vo) };
  }

  /** 显示系(月心) → 会合系(无量纲位置)。 */
  function displayToRotPos(rDisp, tSec) {
    var th = OMEGA_EM * tSec;
    var c = Math.cos(th), s = Math.sin(th);
    /* 显示系 → 轨道面系：v_orbit = M·v_disp（M 行存放基矢，M·v 即逐行点乘得基矢分量） */
    var ro = matVec3(M_ORBIT2DISP, rDisp);
    /* 轨道面系(月心相对) → 质心系无量纲 */
    var X = ro[0] + (1 - MU) * A_EM * c;
    var Y = ro[1] + (1 - MU) * A_EM * s;
    return [(X * c + Y * s) / A_EM, (-X * s + Y * c) / A_EM, ro[2] / A_EM];
  }

  /** 显示系(月心惯性) 完整状态 → 会合系无量纲状态 [x,y,z,vx,vy,vz]（速度对 τ=ωt）：
      rotStateToDisplay 的解析逆变换，用于把开普勒根数初值送入三体传播。
      推导：ro = l·R(θ)ρ − (1−μ)l·[c,s,0]；vo = lω·[−w+a+(1−μ)s, u+b−(1−μ)c, vz]，
      其中 [u,w] = R(θ)ρ 的位置分量、[a,b] = R(−θ)ρ̇。 */
  function displayToRotState(rDisp, vDisp, tSec) {
    var th = OMEGA_EM * tSec;
    var c = Math.cos(th), s = Math.sin(th);
    var l = A_EM, mu = MU, w = OMEGA_EM * l;
    var m = M_ORBIT2DISP;
    /* 显示系 → 轨道面系（M 行存放基矢，M·v 逐行点乘得基矢分量） */
    var ro = matVec3(m, rDisp);
    var vo = matVec3(m, vDisp);
    var u = ro[0] / l + (1 - mu) * c;        /* = x c − y s */
    var ww = ro[1] / l + (1 - mu) * s;       /* = x s + y c */
    var x = u * c + ww * s, y = -u * s + ww * c, z = ro[2] / l;
    var a = vo[0] / w + ww - (1 - mu) * s;   /* R(−θ)ρ̇ 分量 */
    var b = vo[1] / w - u + (1 - mu) * c;
    var vx = a * c + b * s, vy = -a * s + b * c;
    var vz = vo[2] / w;
    return [x, y, z, vx, vy, vz];
  }

  /** 平动点 L1/L2/L3 的显示系位置（km，随 θ 旋转）。 */
  function lpointDisplay(which, tSec) {
    var x = which === 'L1' ? X_L1 : which === 'L2' ? X_L2 : X_L3;
    return rotStateToDisplay([x, 0, 0, 0, 0, 0], tSec).rE;
  }

  /* ══ 微分修正（周期轨道求解） ═════════════════════════════════════ */

  function solveLin(A, n, b) {
    /* 高斯消元（带部分主元），A 行主序 n×n，b n；返回解数组 */
    var M = A.slice(), x = b.slice();
    for (var col = 0; col < n; col++) {
      var piv = col, mx = Math.abs(M[col * n + col]);
      for (var r = col + 1; r < n; r++) {
        if (Math.abs(M[r * n + col]) > mx) { mx = Math.abs(M[r * n + col]); piv = r; }
      }
      if (mx < 1e-14) { return null; }
      if (piv !== col) {
        for (var c = 0; c < n; c++) {
          var t = M[col * n + c]; M[col * n + c] = M[piv * n + c]; M[piv * n + c] = t;
        }
        var tb = x[col]; x[col] = x[piv]; x[piv] = tb;
      }
      for (var r2 = col + 1; r2 < n; r2++) {
        var f = M[r2 * n + col] / M[col * n + col];
        if (f === 0) { continue; }
        for (var c2 = col; c2 < n; c2++) { M[r2 * n + c2] -= f * M[col * n + c2]; }
        x[r2] -= f * x[col];
      }
    }
    for (var r3 = n - 1; r3 >= 0; r3--) {
      for (var c3 = r3 + 1; c3 < n; c3++) { x[r3] -= M[r3 * n + c3] * x[c3]; }
      x[r3] /= M[r3 * n + r3];
    }
    return x;
  }

  /**
   * Halo 族微分修正：固定 z0（族参数），变元 (x0, vy0, T/2)，
   * 约束半周期后 y = vx = vz = 0（xz 平面垂直穿越）。
   * seed = {x0, vy0, Thalf}（无量纲）。返回 {ok, x0, vy0, Thalf, iters}。
   *
   * 注 1：N 固定（h = Thalf/N 随变元连续变化），保证约束 F(v) 连续可微。
   * 注 2：带 Armijo 回溯线搜索——L2 邻域不稳定模 e^{λT} 放大残差，
   *   纯牛顿步带阻尼也常过冲跳进另一拓扑域（调试教训见算法文档）。
   */
  function dcHalo(z0, seed) {
    var v = [seed.x0, seed.vy0, seed.Thalf];
    var EPS = 1e-7;
    var N = 600;
    function F(v) {
      var yEnd = integrateRot([v[0], 0, z0, 0, v[1], 0], v[2] / N, N);
      return [yEnd[1], yEnd[3], yEnd[5]];
    }
    for (var it = 0; it < 60; it++) {
      var f = F(v);
      var nrm = Math.hypot(f[0], f[1], f[2]);
      if (nrm < 1e-11) { return { ok: true, x0: v[0], vy0: v[1], Thalf: v[2], iters: it }; }
      /* 数值雅可比：J[i][k] = ∂F_i/∂v_k（行=约束，列=变元） */
      var J = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
      for (var k = 0; k < 3; k++) {
        var vp = v.slice(); vp[k] += EPS;
        var fp = F(vp);
        for (var i = 0; i < 3; i++) { J[i][k] = (fp[i] - f[i]) / EPS; }
      }
      var dx = solveLin([J[0][0], J[0][1], J[0][2], J[1][0], J[1][1], J[1][2], J[2][0], J[2][1], J[2][2]], 3, [-f[0], -f[1], -f[2]]);
      if (!dx) { return { ok: false, iters: it }; }
      /* 回溯线搜索：α 减半直到 |F| 严格下降 */
      var alpha = 1, moved = false;
      for (var ls = 0; ls < 30; ls++) {
        var vn = [v[0] + alpha * dx[0], v[1] + alpha * dx[1], v[2] + alpha * dx[2]];
        if (vn[2] > 0.05 && isFinite(vn[0]) && isFinite(vn[1]) && isFinite(vn[2])) {
          var fn = F(vn);
          if (Math.hypot(fn[0], fn[1], fn[2]) < nrm) { v = vn; moved = true; break; }
        }
        alpha /= 2;
      }
      if (!moved) { return { ok: false, iters: it }; }
    }
    return { ok: false, iters: 60 };
  }

  /**
   * DRO 平面微分修正：固定半周期 Thalf，变元 (x0, vy0)，
   * 约束 y(Th)=0, vx(Th)=0（x 轴垂直穿越）。seed = {x0, vy0}。
   * 同 dcHalo：固定 N + 回溯线搜索。
   */
  function dcDro(Thalf, seed) {
    var v = [seed.x0, seed.vy0];
    var EPS = 1e-7;
    var N = 600;
    var h = Thalf / N;
    function F(v) {
      var yEnd = integrateRot([v[0], 0, 0, 0, v[1], 0], h, N);
      return [yEnd[1], yEnd[3]];
    }
    for (var it = 0; it < 60; it++) {
      var f = F(v);
      if (Math.hypot(f[0], f[1]) < 1e-11) {
        return { ok: true, x0: v[0], vy0: v[1], Thalf: Thalf, iters: it };
      }
      var J = [[0, 0], [0, 0]];
      for (var k = 0; k < 2; k++) {
        var vp = v.slice(); vp[k] += EPS;
        var fp = F(vp);
        for (var i = 0; i < 2; i++) { J[i][k] = (fp[i] - f[i]) / EPS; }
      }
      var dx = solveLin([J[0][0], J[0][1], J[1][0], J[1][1]], 2, [-f[0], -f[1]]);
      if (!dx) { return { ok: false, iters: it }; }
      var alpha = 1, moved = false;
      for (var ls = 0; ls < 30; ls++) {
        var vn = [v[0] + alpha * dx[0], v[1] + alpha * dx[1]];
        if (isFinite(vn[0]) && isFinite(vn[1])) {
          var fn = F(vn);
          if (Math.hypot(fn[0], fn[1]) < Math.hypot(f[0], f[1])) { v = vn; moved = true; break; }
        }
        alpha /= 2;
      }
      if (!moved) { return { ok: false, iters: it }; }
    }
    return { ok: false, iters: 60 };
  }

  /**
   * 平面 Lyapunov 轨道（L2 邻域，halo 族的分岔母族）：
   * 固定 x0，变元 (vy0, Thalf)，约束 y(Th)=0, vx(Th)=0。
   * 用于族延拓起点，返回 {ok, vy0, Thalf, iters}。
   */
  function dcPlanarLyapunov(x0, seed) {
    var v = [seed.vy0, seed.Thalf];
    var EPS = 1e-7;
    var N = 1200;
    function F(v) {
      var yEnd = integrateRot([x0, 0, 0, 0, v[0], 0], v[1] / N, N);
      return [yEnd[1], yEnd[3]];
    }
    for (var it = 0; it < 60; it++) {
      var f = F(v);
      if (Math.hypot(f[0], f[1]) < 1e-11) {
        return { ok: true, vy0: v[0], Thalf: v[1], iters: it };
      }
      var Jv = [0, 0], Jt = [0, 0];
      var fp1 = F([v[0] + EPS, v[1]]);
      Jv[0] = (fp1[0] - f[0]) / EPS; Jv[1] = (fp1[1] - f[1]) / EPS;
      var fp2 = F([v[0], v[1] + EPS]);
      Jt[0] = (fp2[0] - f[0]) / EPS; Jt[1] = (fp2[1] - f[1]) / EPS;
      var dx = solveLin([Jv[0], Jt[0], Jv[1], Jt[1]], 2, [-f[0], -f[1]]);
      if (!dx) { return { ok: false, iters: it }; }
      var alpha = 1, moved = false;
      for (var ls = 0; ls < 30; ls++) {
        var vn = [v[0] + alpha * dx[0], v[1] + alpha * dx[1]];
        if (vn[1] > 0.05 && isFinite(vn[0]) && isFinite(vn[1])) {
          var fn = F(vn);
          if (Math.hypot(fn[0], fn[1]) < Math.hypot(f[0], f[1])) { v = vn; moved = true; break; }
        }
        alpha /= 2;
      }
      if (!moved) { return { ok: false, iters: it }; }
    }
    return { ok: false, iters: 60 };
  }

  /* ══ 固化预设初值（tune_cr3bp_presets.cjs 离线求解，见算法文档） ══
   * 全部无量纲会合系；T 为整周期（无量纲时间 τ）。
   * 求解链：平面 Lyapunov 种子 → z 延拓进 Halo 族 → 双向族步进。
   * 来源锚点与校验：
   *   haloQQ1 —— 鹊桥一号使命轨道：地月 L2 南族 Halo，Az ≈ 13000 km，
   *              T = 14.79 d（B2 高珊等 2019：Ax≈12500 / Ay≈38000 / T≈14 d；
   *              解得 Ax=12223 / Ay=35311，圆化 CR3BP 与真实系统的正常差异）；
   *   nrho92  —— Gateway 类 L2 南族 9:2 共振 NRHO，T = 6.5644 d
   *              （目标 2×29.5306/9 = 6.5624 d；近月 1517 km，远月 69498 km）；
   *   dro21   —— 2:1 共振 DRO，T = 13.6608 d（恒星月/2），
   *              月心距 68858–92984 km，与月面平均距离 80589 km
   *              （B10 §2.1.3：80592 km；Ax=70523 / Ay=92917）；
   *   dro31   —— 3:1 共振 DRO，T = 恒星月/3 = 9.1069 d（τ 单位 2π/3），
   *              月心距 49248–58066 km，DRO 族延拓打靶解
   *              （dro21 出发沿族 T 递减，半周期正交穿越牛顿打靶，
   *              整周期闭合 ~4e-8；近月小环本征不稳定 ×17/周期 →
   *              传播器按理想周期参考轨道折叠显示，见 foldT）。
   *              注：按角速度比约定（ω_syn = -3 ω_EM）；会合系绕月 1 圈、
   *              惯性系绕月 2 圈/恒星月。无 4:1 预设（2026-10-15 用户移除）。 */
  var PRESET_IC = {
    haloQQ1: { x0: 1.17973572, z0: 0.03381894, vy0: -0.16214952, T: 3.40613902, azKm: 13000 },
    nrho92:  { x0: 1.02205629, z0: 0.18213144, vy0: -0.10333296, T: 1.51167822 },
    dro21:   { x0: 1.17564228, z0: 0, vy0: -0.49445949, T: 3.14585927 },
    dro31:   { x0: 1.12017939, z0: 0, vy0: -0.46178078, T: 2 * Math.PI / 3 },
  };

  /* ══ 传播器一：开普勒二体（月球） ═════════════════════════════════ */

  function MoonKeplerPropagator(el) {
    C.KeplerPropagator.call(this, el, MU_MOON, R_MOON);
  }
  MoonKeplerPropagator.prototype = Object.create(C.KeplerPropagator.prototype);
  MoonKeplerPropagator.id = 'kepler-moon';
  MoonKeplerPropagator.label = '开普勒二体（月球）';
  MoonKeplerPropagator.hasPerturbation = false;

  /* ══ 传播器二：月心 Cowell 数值摄动 ═══════════════════════════════ */

  var MOON_DEFAULTS = {
    grav: true,      // 月球非球形引力（J2–J4，GL0660B）
    earth3B: true,   // 地球第三体（直接项+间接项）
    srp: true,       // 太阳光压（圆柱月影）
    CR: 1.5,         // 光压系数
    SM: 0.02,        // 面质比 m²/kg
    windowS: WINDOW_S,
    hStep: 0,        // 0 = T₀/120 自动（钳 [60, 600] s）
  };

  /** 月心摄动加速度分解（显示系 km/s²）。 */
  function moonAccelBreakdown(y, tSec, o) {
    var r = [y[0], y[1], y[2]], v = [y[3], y[4], y[5]];
    var rm = Math.hypot(r[0], r[1], r[2]);
    var grav = [
      -MU_MOON * r[0] / (rm * rm * rm),
      -MU_MOON * r[1] / (rm * rm * rm),
      -MU_MOON * r[2] / (rm * rm * rm),
    ];
    var zonal = [0, 0, 0];
    if (o.grav) {
      moonZonalAccel(r, zonal);
      zonal[0] -= grav[0]; zonal[1] -= grav[1]; zonal[2] -= grav[2];
    }
    var earth3B = [0, 0, 0];
    if (o.earth3B) {
      var rE = earthPosMoonTab(tSec);
      var dx = r[0] - rE[0], dy = r[1] - rE[1], dz = r[2] - rE[2];
      var d2 = dx * dx + dy * dy + dz * dz, dm = Math.sqrt(d2);
      var e2 = dot3(rE, rE), em = Math.sqrt(e2);
      var k1 = -MU_EARTH / (d2 * dm), k2 = -MU_EARTH / (e2 * em);
      earth3B = [k1 * dx + k2 * rE[0], k1 * dy + k2 * rE[1], k1 * dz + k2 * rE[2]];
    }
    var srp = [0, 0, 0];
    if (o.srp) {
      var sHat = sunDirMoonTab(tSec);
      var nu = moonShadowFactor(r, sHat);
      if (nu > 0) {
        var ks = -nu * o.CR * o.SM * P_SUN / 1000;
        srp = [ks * sHat[0], ks * sHat[1], ks * sHat[2]];
      }
    }
    var total = [0, 0, 0];
    for (var i = 0; i < 3; i++) {
      total[i] = grav[i] + zonal[i] + earth3B[i] + srp[i];
    }
    return { grav: grav, zonal: zonal, earth3B: earth3B, srp: srp, total: total };
  }

  function moonDeriv(y, tSec, o, out) {
    var b = moonAccelBreakdown(y, tSec, o);
    out[0] = y[3]; out[1] = y[4]; out[2] = y[5];
    out[3] = b.total[0]; out[4] = b.total[1]; out[5] = b.total[2];
    return out;
  }

  function MoonCowellPropagator(el, opts) {
    this.isPerturbed = true;
    this.isMoonCowell = true;
    this.o = {};
    for (var k in MOON_DEFAULTS) { this.o[k] = MOON_DEFAULTS[k]; }
    if (opts) { for (var k2 in opts) { this.o[k2] = opts[k2]; } }

    this.el0 = { a: el.a, e: el.e, i: el.i, Om: el.Om, w: el.w };
    /* 积分历元恒锚近月点（ν=0）：内部 t=0 ↔ 近月点，与 timeAtTrueAnomaly
       及页面 tSincePeri（自近月点起算）同基准。el.nu 仅由页面 rebuild 换算成
       tSincePeri 初值，不参与积分初值——否则 y0@ν₀ 与近月点基准错位，
       rebuild（拖 Ω/ω/a/e/i）后相位跳变。 */
    this.nu0 = 0;

    var st0 = C.classicalToState(this.el0, 0, MU_MOON, R_MOON);
    this.y0 = [st0.rE[0], st0.rE[1], st0.rE[2], st0.vE[0], st0.vE[1], st0.vE[2]];

    var T0 = C.period(this.el0.a, MU_MOON);
    /* 步长：T₀/240，钳 [30, 300] s。大偏心率 ELFO 近月点通过时标 r/v ≈ 1100 s，
     * 600 s 级步长会产生系统性 e 衰减（调试教训见算法文档）；30–300 s
     * 覆盖 12 h 圆轨道（T/240=180 s）到 24 h 大椭圆（T/240=360→300 s）。 */
    this.h = this.o.hStep > 0 ? this.o.hStep
      : Math.min(300, Math.max(30, T0 / 240));

    /* RK4 全窗口落盘（与 elements_perturb.propagate 同型，含显式 t） */
    var win = this.o.windowS;
    var n = Math.max(1, Math.ceil(win / this.h));
    var h = win / n;                    // 均分覆盖整窗
    this.h = h;
    var t = new Float64Array(n + 1);
    var y = new Float64Array((n + 1) * 6);
    var yc = Float64Array.from(this.y0);
    y.set(yc, 0); t[0] = 0;
    var k1 = new Float64Array(6), k2 = new Float64Array(6),
        k3 = new Float64Array(6), k4 = new Float64Array(6), yt = new Float64Array(6);
    for (var i = 0; i < n; i++) {
      var ti = i * h;
      moonDeriv(yc, ti, this.o, k1);
      for (var j = 0; j < 6; j++) { yt[j] = yc[j] + 0.5 * h * k1[j]; }
      moonDeriv(yt, ti + 0.5 * h, this.o, k2);
      for (var j2 = 0; j2 < 6; j2++) { yt[j2] = yc[j2] + 0.5 * h * k2[j2]; }
      moonDeriv(yt, ti + 0.5 * h, this.o, k3);
      for (var j3 = 0; j3 < 6; j3++) { yt[j3] = yc[j3] + h * k3[j3]; }
      moonDeriv(yt, ti + h, this.o, k4);
      for (var j4 = 0; j4 < 6; j4++) {
        yc[j4] += h / 6 * (k1[j4] + 2 * k2[j4] + 2 * k3[j4] + k4[j4]);
      }
      t[i + 1] = (i + 1) * h;
      y.set(yc, (i + 1) * 6);
    }
    this.nSteps = n; this.tArr = t; this.yArr = y;
    this.window = n * h;
    this.tNow = 0;
    this._osc = cartToClassicalMu(this.y0.slice(0, 3), this.y0.slice(3, 6), MU_MOON);
  }

  MoonCowellPropagator.prototype.stateVecAt = function (tSec) {
    var t = Math.min(Math.max(tSec, 0), this.window);
    var x = t / this.h;
    var i = Math.min(Math.floor(x), this.nSteps - 1);
    var tau = (t - i * this.h) / this.h;
    var y = this.yArr, o0 = i * 6, o1 = (i + 1) * 6, h = this.h;
    var t2 = tau * tau, t3 = t2 * tau;
    var h00 = 2 * t3 - 3 * t2 + 1, h10 = t3 - 2 * t2 + tau;
    var h01 = -2 * t3 + 3 * t2, h11 = t3 - t2;
    var dh00 = 6 * t2 - 6 * tau, dh10 = 3 * t2 - 4 * tau + 1;
    var dh01 = -6 * t2 + 6 * tau, dh11 = 3 * t2 - 2 * tau;
    var out = new Array(6);
    for (var j = 0; j < 3; j++) {
      var p0 = y[o0 + j], v0 = y[o0 + 3 + j], p1 = y[o1 + j], v1 = y[o1 + 3 + j];
      out[j] = h00 * p0 + h10 * h * v0 + h01 * p1 + h11 * h * v1;
      out[3 + j] = (dh00 * p0 + dh10 * h * v0 + dh01 * p1 + dh11 * h * v1) / h;
    }
    return out;
  };

  MoonCowellPropagator.prototype.stateAtTime = function (tSec) {
    var y = this.stateVecAt(tSec);
    var rE = [y[0], y[1], y[2]], vE = [y[3], y[4], y[5]];
    var osc = cartToClassicalMu(rE, vE, MU_MOON);
    var rm = Math.hypot(rE[0], rE[1], rE[2]);
    var vm = Math.hypot(vE[0], vE[1], vE[2]);
    var rHat = [rE[0] / rm, rE[1] / rm, rE[2] / rm];
    var vr = dot3(vE, rHat);
    var h = cross3(rE, vE), hm = Math.hypot(h[0], h[1], h[2]);
    var hunit = [h[0] / hm, h[1] / hm, h[2] / hm];
    var tHat = cross3(hunit, rHat);
    var vt = dot3(vE, tHat);
    var E = C.nuToE(osc.nu, Math.max(osc.e, 1e-12));
    return {
      rE: rE, vE: vE, r: rm, v: vm, alt: rm - R_MOON,
      vr: vr, vt: vt, gamma: Math.atan2(vr, vt),
      osc: osc, E: E, M: E - osc.e * Math.sin(E), tSec: tSec,
    };
  };

  MoonCowellPropagator.prototype.trueAnomalyAtTime = function (tSec) {
    this.tNow = Math.min(Math.max(tSec, 0), this.window);
    this._osc = this.stateAtTime(this.tNow).osc;
    return this._osc.nu;
  };
  MoonCowellPropagator.prototype.timeAtTrueAnomaly = function (nuRad) {
    var M = C.nuToM(nuRad, this.el0.e).M;
    return M / C.meanMotion(this.el0.a, MU_MOON);
  };
  MoonCowellPropagator.prototype.stateAtTrueAnomaly = function () {
    return this.stateAtTime(this.tNow);
  };
  MoonCowellPropagator.prototype.periodSeconds = function () {
    return C.period(this._osc.a, MU_MOON);
  };
  MoonCowellPropagator.prototype.osculatingElements = function () { return this._osc; };
  MoonCowellPropagator.prototype.windowSeconds = function () { return this.window; };
  MoonCowellPropagator.prototype.stepSeconds = function () { return this.h; };

  MoonCowellPropagator.prototype.forceMagnitudes = function () {
    var y = this.stateVecAt(this.tNow);
    var b = moonAccelBreakdown(y, this.tNow, this.o);
    var mag = function (a) { return Math.hypot(a[0], a[1], a[2]); };
    return { grav: mag(b.grav), zonal: mag(b.zonal), earth3B: mag(b.earth3B), srp: mag(b.srp) };
  };

  MoonCowellPropagator.prototype.trailPoints = function (maxPts) {
    var stride = Math.max(1, Math.ceil((this.nSteps + 1) / (maxPts || 3000)));
    var x = [], y = [], z = [];
    for (var i = 0; i <= this.nSteps; i += stride) {
      var o = i * 6;
      x.push(this.yArr[o]); y.push(this.yArr[o + 1]); z.push(this.yArr[o + 2]);
    }
    var last = this.nSteps * 6;
    if (x[x.length - 1] !== this.yArr[last]) {
      x.push(this.yArr[last]); y.push(this.yArr[last + 1]); z.push(this.yArr[last + 2]);
    }
    return { x: x, y: y, z: z, stride: stride };
  };

  MoonCowellPropagator.id = 'perturbed-moon';
  MoonCowellPropagator.label = '数值摄动（月球 J2–J4 + 地球三体 + 光压，90 d）';
  MoonCowellPropagator.hasPerturbation = true;

  /* ══ 传播器三：CR3BP 会合系 ═══════════════════════════════════════ */

  /**
   * constructor(opts)：
   *   preset = 'halo' | 'nrho' | 'dro' | 'elements'
   *   azKm   = Halo 族 z 振幅（km，仅 halo；从预设解续修正）
   *   phase  = 0..1 起始相位（整周期内偏移）
   *   elements = {a,e,i,Om,w,nu}（弧度；仅 elements 模式：由根数构建历元初值，
   *            如 ELFO —— S13 冻结设计核心即地球三体，三体模式下直接传播）
   *   windowS = 传播窗口（默认 90 d）
   * 构造时完成全窗口 RK4（会合系无量纲），帧间 Hermite 插值，
   * 输出换算到显示系（月心 km）。
   */
  function Cr3bpPropagator(opts) {
    this.isPerturbed = true;
    this.isCr3bp = true;
    var o = opts || {};
    this.preset = o.preset || 'nrho';
    this.phase = (o.phase !== undefined) ? Math.min(1, Math.max(0, o.phase)) : 0;
    var win = (o.windowS !== undefined) ? o.windowS : WINDOW_S;

    var ic, y0;
    if (this.preset === 'elements') {
      /* 根数初值模式：开普勒根数 → 显示系惯性状态 → 会合系无量纲状态 */
      var el = o.elements || {};
      var st0 = C.classicalToState(
        { a: el.a, e: el.e, i: el.i, Om: el.Om, w: el.w },
        (el.nu !== undefined) ? el.nu : 0, MU_MOON, R_MOON);
      y0 = displayToRotState(st0.rE, st0.vE, 0);
      ic = { T: C.period(el.a, MU_MOON) * OMEGA_EM };   /* τ 单位的开普勒周期 */
    } else {
      if (this.preset === 'halo') {
        var az = (o.azKm !== undefined) ? o.azKm : PRESET_IC.haloQQ1.azKm;
        ic = this._solveHaloByAz(az);
        if (!ic) { ic = cloneIc(PRESET_IC.haloQQ1); }
      } else if (this.preset === 'dro' || this.preset === 'dro31') {
        var droIc = { dro: PRESET_IC.dro21, dro31: PRESET_IC.dro31 }[this.preset];
        ic = cloneIc(droIc);
      } else {
        ic = cloneIc(PRESET_IC.nrho92);
      }
      y0 = [ic.x0, 0, ic.z0, 0, ic.vy0, 0];
    }
    this.ic = ic;
    /* 理想周期参考轨道折叠（dro31）：近月小环 DRO 本征不稳定
       （×17 每周期），全窗口 90 d 预积分会被数值噪声指数放大而飞离
       （dro21 稳定漂移仅 5.7e-5 无此问题）。改为只预积分一个周期，
       stateVecAt 按 t mod T 折叠采样——显示恒为精确周期参考轨道。
       orbitPath 同样走 stateVecAt 折叠采样（见该函数注释）；
       时间轴（windowSeconds）仍为 90 d。 */
    this.foldT = (this.preset === 'dro31') ? ic.T : 0;
    /* 相位偏移前的完整会合系初值：elements 模式 ic 仅含 T（无 x0/z0/vy0），
       orbitPath/missionParams 必须从本状态出发积分，否则初值全 NaN */
    this.y0base = y0.slice();
    this.T = ic.T / OMEGA_EM;            // s
    this.periodT = this.T;

    /* 相位偏移后的起始状态 */
    if (this.phase > 0) {
      var nPh = Math.max(50, Math.ceil(ic.T * this.phase / 0.005));
      y0 = integrateRot(y0, ic.T * this.phase / nPh, nPh);
    }

    /* 全窗口积分（落盘会合系状态）。步长按预设定：halo/dro T/600、
     * nrho T/1200（近月点 ~1500 km，局部动力学快，需更高分辨率；
     * 90 d 实测雅可比漂移 halo 6.7e-11 / nrho 3.5e-4 / dro 5.7e-5）；
     * elements T/2400（ELFO 近月点高度 ~440 km：T/600 的 CJ 漂移 1.9e-3、
     * 近月点偏低 ~37 km 未收敛；T/2400 漂移 1.9e-6、与 T/38400 差 ~6 km，
     * 视为收敛。代价 216k 步 ≈ Cowell 摄动 26k 步的耗时量级）。 */
    /* 轨道修正（定轨保持）：仅 Halo 周期轨道。沿参考轨道（当前 IC 的周期解）
     * 预计算一个周期 8 段的参考相位点 skRef[k] 与段 STM skStm[k]（T/4800
     * 耦合积分 Φ'=AΦ）；全窗口积分时每当轨道相位越过边界 k·T/8，即在边界
     * 时刻（步内子步精确到达，避免沿轨伪差）施加单脉冲死拍控制
     *   Δv = −Φrv⁻¹(Φrr·δr + Φrv·δv)   （δ 为相对参考相位点的偏差）
     * 把下一段末位置偏差线性归零。背景：Halo 本征指数不稳定（λ≈0.0182/h，
     * 每周期 ×655），无修正时 ~1200 h 必然漂离（见漂移归因实验）；修正后
     * 90 d 保持在轨道上（段边界残差 km 量级，Δv 毫米/秒级）。 */
    this.stationKeeping = !!o.stationKeeping && this.preset === 'halo';
    this.corrections = [];
    var skRef = null, skStm = null, segTau = ic.T / 8;
    if (this.stationKeeping) {
      skRef = []; skStm = [];
      var s42 = new Float64Array(42);
      s42.set(this.y0base, 0);
      for (var q = 0; q < 6; q++) { s42[6 + q * 6 + q] = 1; }
      var hf = ic.T / 4800;
      skRef.push(Array.from(s42.slice(0, 6)));
      for (var sg = 0; sg < 8; sg++) {
        for (var si = 0; si < 600; si++) { rk4StepStm(s42, hf); }
        skRef.push(Array.from(s42.slice(0, 6)));
        skStm.push(Array.from(s42.slice(6, 42)));
      }
    }

    var winTau = this.foldT || (win * OMEGA_EM);   /* 折叠模式只预积分一个周期 */
    var stepsPerPeriod = { nrho: 1200, halo: 600, dro: 600, dro31: 600 }[this.preset] || 2400;
    var hTau = ic.T / stepsPerPeriod;
    var n = Math.max(1, Math.ceil(winTau / hTau));
    hTau = winTau / n;
    var t = new Float64Array(n + 1);
    var y = new Float64Array((n + 1) * 6);
    var yc = Float64Array.from(y0);
    y.set(yc, 0); t[0] = 0;
    var vUnitMs = A_EM * OMEGA_EM * 1000;   /* 无量纲速度 → m/s */
    var nextB = this.stationKeeping ? Math.floor(this.phase * 8 + 1e-9) + 1 : Infinity;
    for (var i = 0; i < n; i++) {
      var tau0 = i * hTau, tau1 = (i + 1) * hTau;
      /* 边界仿真时刻：轨道相位 = phase·T + t 越过 m·T/8 处（t_b = (m−8φ)·T/8） */
      var tb = (nextB - this.phase * 8) * segTau;
      if (this.stationKeeping && tb > tau0 + 1e-9 && tb <= tau1 + 1e-9) {
        /* 边界落在本步内：子步精确积到边界 → 施加修正 → 积完本步 */
        if (tb - tau0 > 1e-12) { rk4Step(yc, tb - tau0); }
        var seg = ((nextB % 8) + 8) % 8;
        var ref = skRef[seg], P = skStm[seg];
        var b3 = [0, 0, 0];
        for (var r = 0; r < 3; r++) {
          var sb = 0;
          for (var c = 0; c < 3; c++) {
            sb += P[r * 6 + c] * (yc[c] - ref[c]) + P[r * 6 + 3 + c] * (yc[3 + c] - ref[3 + c]);
          }
          b3[r] = sb;
        }
        var dv = solveLin([P[3], P[4], P[5], P[9], P[10], P[11], P[15], P[16], P[17]],
                          3, [-b3[0], -b3[1], -b3[2]]);
        if (dv) {
          yc[3] += dv[0]; yc[4] += dv[1]; yc[5] += dv[2];
          this.corrections.push({
            tSec: tb / OMEGA_EM,
            dvMps: Math.hypot(dv[0], dv[1], dv[2]) * vUnitMs,
            devKm: Math.hypot(yc[0] - ref[0], yc[1] - ref[1], yc[2] - ref[2]) * A_EM,
          });
        }
        nextB++;
        if (tau1 - tb > 1e-12) { rk4Step(yc, tau1 - tb); }
      } else {
        rk4Step(yc, hTau);
      }
      t[i + 1] = tau1 / OMEGA_EM;
      y.set(yc, (i + 1) * 6);
    }
    this.nSteps = n; this.hTau = hTau; this.h = hTau / OMEGA_EM;
    this.tArr = t; this.yArr = y;
    this.window = n * this.h;
    this.tNow = 0;
    this._mp = null;
  }

  function cloneIc(ic) {
    return { x0: ic.x0, z0: ic.z0, vy0: ic.vy0, T: ic.T, azKm: ic.azKm };
  }

  /** Az 续修正：从预设解出发，按 500 km 步长走向目标 Az。 */
  Cr3bpPropagator.prototype._solveHaloByAz = function (azKm) {
    var base = PRESET_IC.haloQQ1;
    var z0 = -Math.abs(azKm) / A_EM;          // 南族 z0 < 0（远端穿越）
    var cur = { x0: base.x0, vy0: base.vy0, Thalf: base.T / 2 };
    var azCur = base.azKm;
    var nSteps = Math.max(1, Math.ceil(Math.abs(azKm - azCur) / 500));
    for (var k = 1; k <= nSteps; k++) {
      var azTry = azCur + (azKm - azCur) * k / nSteps;
      var zTry = -Math.abs(azTry) / A_EM;
      var sol = dcHalo(zTry, cur);
      if (!sol.ok) { return null; }
      cur = { x0: sol.x0, vy0: sol.vy0, Thalf: sol.Thalf };
    }
    return { x0: cur.x0, z0: z0, vy0: cur.vy0, T: 2 * cur.Thalf, azKm: azKm };
  };

  Cr3bpPropagator.prototype.stateVecAt = function (tSec) {
    var t = Math.min(Math.max(tSec, 0), this.window);
    /* 理想周期轨道：相位折叠采样。量纲：foldT 是 τ 单位（ic.T），t 是秒——
       必须换算回秒再取模（2026-10-15 修复：原 t % foldT 秒对 τ 取模，余数恒
       <2.1 s → 任意时刻都返回 y(0)，dro31 惯性线退化为正圆、会合系环退化为
       点、旋转系视角下航天器静止如「卡死」） */
    if (this.foldT) { t = t % (this.foldT / OMEGA_EM); }
    var x = t / this.h;
    var i = Math.min(Math.floor(x), this.nSteps - 1);
    var tau = (t - i * this.h) / this.h;
    var y = this.yArr, o0 = i * 6, o1 = (i + 1) * 6, h = this.h;
    var t2 = tau * tau, t3 = t2 * tau;
    var h00 = 2 * t3 - 3 * t2 + 1, h10 = t3 - 2 * t2 + tau;
    var h01 = -2 * t3 + 3 * t2, h11 = t3 - t2;
    var dh00 = 6 * t2 - 6 * tau, dh10 = 3 * t2 - 4 * tau + 1;
    var dh01 = -6 * t2 + 6 * tau, dh11 = 3 * t2 - 2 * tau;
    var out = new Array(6);
    /* 会合系速度（对 τ 求导）→ 对 t 的导数 = ω·(d/dτ) */
    for (var j = 0; j < 3; j++) {
      var p0 = y[o0 + j], v0 = y[o0 + 3 + j] * OMEGA_EM, p1 = y[o1 + j], v1 = y[o1 + 3 + j] * OMEGA_EM;
      out[j] = h00 * p0 + h10 * h * v0 + h01 * p1 + h11 * h * v1;
      out[3 + j] = (dh00 * p0 + dh10 * h * v0 + dh01 * p1 + dh11 * h * v1) / h;
    }
    return out;   // [x,y,z, ẋ,ẏ,ż] 会合系无量纲（速度已对 t）
  };

  /** tSec 处完整状态（显示系 km/km/s + 会合系读数 + 雅可比积分）。 */
  Cr3bpPropagator.prototype.stateAtTime = function (tSec) {
    var yv = this.stateVecAt(tSec);
    var yRot = [yv[0], yv[1], yv[2], yv[3] / OMEGA_EM, yv[4] / OMEGA_EM, yv[5] / OMEGA_EM];
    var disp = rotStateToDisplay(yRot, tSec);
    var rE = disp.rE, vE = disp.vE;
    var osc = cartToClassicalMu(rE, vE, MU_MOON);
    var rm = Math.hypot(rE[0], rE[1], rE[2]);
    var vm = Math.hypot(vE[0], vE[1], vE[2]);
    var rHat = [rE[0] / rm, rE[1] / rm, rE[2] / rm];
    var vr = dot3(vE, rHat);
    var h = cross3(rE, vE), hm = Math.hypot(h[0], h[1], h[2]);
    var hunit = [h[0] / hm, h[1] / hm, h[2] / hm];
    var tHat = cross3(hunit, rHat);
    var vt = dot3(vE, tHat);
    var r1 = Math.hypot(yRot[0] + MU, yRot[1], yRot[2]);
    var r2 = Math.hypot(yRot[0] - 1 + MU, yRot[1], yRot[2]);
    return {
      rE: rE, vE: vE, r: rm, v: vm, alt: rm - R_MOON,
      vr: vr, vt: vt, gamma: Math.atan2(vr, vt),
      osc: osc, tSec: tSec,
      rot: { x: yRot[0], y: yRot[1], z: yRot[2], vx: yRot[3], vy: yRot[4], vz: yRot[5] },
      cj: jacobi(yRot), r1: r1, r2: r2,
    };
  };

  Cr3bpPropagator.prototype.trueAnomalyAtTime = function (tSec) {
    this.tNow = Math.min(Math.max(tSec, 0), this.window);
    return (this.tNow / this.T) * C.TWO_PI;   // 相位代理（供 UI 锚定）
  };
  Cr3bpPropagator.prototype.timeAtTrueAnomaly = function (nuRad) {
    return (C.wrap2pi(nuRad) / C.TWO_PI) * this.T;
  };
  Cr3bpPropagator.prototype.stateAtTrueAnomaly = function () {
    return this.stateAtTime(this.tNow);
  };
  Cr3bpPropagator.prototype.periodSeconds = function () { return this.T; };
  /** 惯性系轨迹图案的完整跨度（秒）：DRO n:1 共振在惯性系 nT 闭合
      （会合系周期 T = 恒星月/n，惯性图案 = 1 个恒星月 = nT：ρ(t) 周期 T、
      旋转变换 R(ωt) 周期 nT，乘积周期 lcm = nT）。
      Halo/NRHO 惯性系缓慢进动、无短闭合，恒取 T（单圈弧）。 */
  Cr3bpPropagator.prototype.inertialSpanSeconds = function () {
    if (this.preset === 'dro') return 2 * this.T;      /* 2:1 共振 */
    if (this.preset === 'dro31') return 3 * this.T;    /* 3:1 共振 */
    return this.T;
  };
  Cr3bpPropagator.prototype.osculatingElements = function () {
    return this.stateAtTime(this.tNow).osc;
  };
  Cr3bpPropagator.prototype.windowSeconds = function () { return this.window; };
  Cr3bpPropagator.prototype.stepSeconds = function () { return this.h; };

  /** 一个整周期的显示系轨迹（轨道线绘制用）。
      syn=true 时输出会合系（旋转坐标系）静态轨迹：所有采样点统一用
      θ=0 的变换，地月连线固定——L2 视角下卫星绕 L2 转圈的闭合环。
      spanSec：可选采样跨度（秒），默认 syn ? T : inertialSpanSeconds()。
      页面地心轨迹线（geoLine）传恒星月 TWO_PI/OMEGA_EM——地心视角绕地球
      一圈的自然周期（2026-10-15 修复 halo 绿线画半圈即停：原沿用
      inertialSpanSeconds= T = 54% 恒星月，nrho 仅 24%）。 */
  Cr3bpPropagator.prototype.orbitPath = function (nPts, syn, spanSec) {
    var n = nPts || 256;
    var ic = this.ic;
    var pts = [];
    /* 周期轨道（dro/dro31）：直接按折叠采样 yArr 输出，与航天器（stateVecAt）
       同一数据源。原先的二次大步长 RK4 积分对本征不稳定环（dro31 ×17/周期）
       会把步长误差指数放大，画出的环变形错位——页面上表现为航天器与轨道线
       方向/位置对不上（dro21 稳定故此前未暴露）。 */
    if (this.preset === 'dro' || this.preset === 'dro31') {
      /* 惯性轨迹采样完整图案跨度（DRO 2:1 → 2T、3:1 → 3T，惯性系 nT 闭合），
         否则页面只显示前半个惯性周期（2026-10-15 用户反馈「固定残存半个周期」）；
         会合系环仍以 T 闭合，不变。 */
      var Tp = syn ? this.periodSeconds() : (spanSec || this.inertialSpanSeconds());
      for (var i = 0; i <= n; i++) {
        var tSec = i / n * Tp;
        var yv = this.stateVecAt(tSec);
        var yR = [yv[0], yv[1], yv[2], yv[3] / OMEGA_EM, yv[4] / OMEGA_EM, yv[5] / OMEGA_EM];
        var dd = syn ? rotStateToDisplay(yR, 0) : rotStateToDisplay(yR, tSec);
        pts.push(dd.rE);
      }
      return pts;
    }
    /* 从相位偏移前的完整会合系初值出发（elements 模式 ic 仅含 T）。
       spanSec 可超出 ic.T（halo/nrho 地心轨迹需跨恒星月，积分窗 90 d 足够）：
       步长 h = span/nSteps，密度与预积分（T/600）同量级，RK4 精度不受损。 */
    var y = this.y0base.slice();
    var total = spanSec ? spanSec * OMEGA_EM : ic.T;   /* τ 单位总跨度 */
    var nSteps = Math.max(240, n * 2);
    var h = total / nSteps;
    for (var i = 0; i <= nSteps; i++) {
      var d = syn ? rotStateToDisplay(y, 0) : rotStateToDisplay(y, i * h / OMEGA_EM);
      pts.push(d.rE);
      if (i < nSteps) {
        /* 单步 RK4 推进（逐步落点，保持与 integrateRot 一致精度） */
        var k1 = cr3bpDeriv(y, [0,0,0,0,0,0]);
        var yt = y.map((v, k) => v + 0.5 * h * k1[k]);
        var k2 = cr3bpDeriv(yt, [0,0,0,0,0,0]);
        var yt2 = y.map((v, k) => v + 0.5 * h * k2[k]);
        var k3 = cr3bpDeriv(yt2, [0,0,0,0,0,0]);
        var yt3 = y.map((v, k) => v + h * k3[k]);
        var k4 = cr3bpDeriv(yt3, [0,0,0,0,0,0]);
        y = y.map((v, k) => v + h / 6 * (k1[k] + 2 * k2[k] + 2 * k3[k] + k4[k]));
      }
    }
    return pts;
  };

  /** 使命轨道参数（周期/近月点/远月点/振幅/雅可比）。 */
  Cr3bpPropagator.prototype.missionParams = function () {
    if (this._mp) { return this._mp; }
    var ic = this.ic;
    /* 从相位偏移前的完整会合系初值出发（elements 模式 ic 仅含 T） */
    var y = this.y0base.slice();
    var nSteps = 480;
    var h = ic.T / nSteps;
    var rp = Infinity, ra = -Infinity, xMin = Infinity, xMax = -Infinity,
        yMax = 0, zMax = 0;
    for (var i = 0; i <= nSteps; i++) {
      var r2 = Math.hypot(y[0] - 1 + MU, y[1], y[2]);
      var rpKm = r2 * A_EM - R_MOON;
      if (rpKm < rp) { rp = rpKm; }
      var raKm = r2 * A_EM - R_MOON;
      if (raKm > ra) { ra = raKm; }
      if (y[0] < xMin) { xMin = y[0]; }
      if (y[0] > xMax) { xMax = y[0]; }
      if (Math.abs(y[1]) > yMax) { yMax = Math.abs(y[1]); }
      if (Math.abs(y[2]) > zMax) { zMax = Math.abs(y[2]); }
      if (i < nSteps) {
        var k1 = cr3bpDeriv(y, [0,0,0,0,0,0]);
        var yt = y.map((v, k) => v + 0.5 * h * k1[k]);
        var k2 = cr3bpDeriv(yt, [0,0,0,0,0,0]);
        var yt2 = y.map((v, k) => v + 0.5 * h * k2[k]);
        var k3 = cr3bpDeriv(yt2, [0,0,0,0,0,0]);
        var yt3 = y.map((v, k) => v + h * k3[k]);
        var k4 = cr3bpDeriv(yt3, [0,0,0,0,0,0]);
        y = y.map((v, k) => v + h / 6 * (k1[k] + 2 * k2[k] + 2 * k3[k] + k4[k]));
      }
    }
    this._mp = {
      T_s: ic.T / OMEGA_EM,
      periluneKm: rp, apoluneKm: ra,
      azKm: zMax * A_EM, axKm: (xMax - xMin) / 2 * A_EM, ayKm: yMax * A_EM,
      cj: jacobi(this.y0base),
      cjL1: C_L1, cjL2: C_L2, cjL3: C_L3,
      xL1: X_L1, xL2: X_L2, xL3: X_L3,
    };
    return this._mp;
  };

  /** 拖尾点列（显示系 km）。 */
  Cr3bpPropagator.prototype.trailPoints = function (maxPts) {
    var stride = Math.max(1, Math.ceil((this.nSteps + 1) / (maxPts || 3000)));
    var x = [], y = [], z = [];
    for (var i = 0; i <= this.nSteps; i += stride) {
      var o = i * 6;
      var d = rotStateToDisplay(
        [this.yArr[o], this.yArr[o + 1], this.yArr[o + 2],
         this.yArr[o + 3], this.yArr[o + 4], this.yArr[o + 5]],
        this.tArr[i]);
      x.push(d.rE[0]); y.push(d.rE[1]); z.push(d.rE[2]);
    }
    return { x: x, y: y, z: z, stride: stride };
  };

  Cr3bpPropagator.id = 'cr3bp';
  Cr3bpPropagator.label = '三体会合系（CR3BP：Halo / NRHO / DRO，90 d）';
  Cr3bpPropagator.hasPerturbation = true;

  /* 注册 */
  C.registerPropagator(MoonKeplerPropagator);
  C.registerPropagator(MoonCowellPropagator);
  C.registerPropagator(Cr3bpPropagator);

  return {
    /* 常数 */
    MU_MOON: MU_MOON, R_MOON: R_MOON, R_MOON_GRAV: R_MOON_GRAV,
    J2_MOON: J2_MOON, J3_MOON: J3_MOON, J4_MOON: J4_MOON,
    MU_EARTH: MU_EARTH, MU_EM_TOTAL: MU_EM_TOTAL, MU_RATIO: MU_RATIO,
    A_EM: A_EM, OMEGA_EM: OMEGA_EM, T_EM: T_EM,
    CASSINI_DEG: CASSINI_DEG, THETA_E0: THETA_E0, LAM0_SUN: LAM0_SUN,
    WINDOW_S: WINDOW_S, P_SUN: P_SUN, YEAR_S: YEAR_S,
    X_L1: X_L1, X_L2: X_L2, X_L3: X_L3, C_L1: C_L1, C_L2: C_L2, C_L3: C_L3,
    PRESET_IC: PRESET_IC,
    /* 模型函数 */
    earthPosMoonTab: earthPosMoonTab,
    sunDirMoonTab: sunDirMoonTab,
    moonShadowFactor: moonShadowFactor,
    moonZonalAccel: moonZonalAccel,
    moonAccelBreakdown: moonAccelBreakdown,
    cartToClassicalMu: cartToClassicalMu,
    orbitToDisp: orbitToDisp,
    rotStateToDisplay: rotStateToDisplay,
    displayToRotPos: displayToRotPos,
    displayToRotState: displayToRotState,
    lpointDisplay: lpointDisplay,
    omegaPot: omegaPot,
    cr3bpDeriv: cr3bpDeriv,
    integrateRot: integrateRot,
    jacobi: jacobi,
    collinearX: collinearX,
    dcHalo: dcHalo,
    dcDro: dcDro,
    dcPlanarLyapunov: dcPlanarLyapunov,
    solveLin: solveLin,
    /* 传播器 */
    MoonKeplerPropagator: MoonKeplerPropagator,
    MoonCowellPropagator: MoonCowellPropagator,
    Cr3bpPropagator: Cr3bpPropagator,
    MOON_DEFAULTS: MOON_DEFAULTS,
  };
}));
