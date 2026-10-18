/*!
 * elements_perturb.js — 近地轨道数值摄动传播器（Cowell 特殊摄动法）。
 *
 * 与 elements_core.js 同一 UMD 双端约定（浏览器挂 window.ElementsPerturb，
 * node 走 module.exports），并在加载时自动向 ElementsCore 注册传播器，
 * elements_demo.html 的「传播模型」下拉框零改动即可识别。
 *
 * 力学模型（公式编号均指国防科技大学《轨道动力学》教材，下称"教材"）：
 *   - 中心引力 + 地球非球形带谐项 J2/J3/J4（教材 7.4 节，位函数
 *     U = μ/r·[1 − Σ J_n·(R⊕/r)^n·P_n(sinφ)]，加速度取其梯度；
 *     J2 分量与教材式(7-4-33)闭式等价——node 自检中有逐项比对）
 *   - 大气阻力：a = −½·C_D·(S/m)·ρ·|v_rel|·v_rel（教材式 7-5-3/7-5-14），
 *     v_rel = v − ω⊕×r（大气随地球刚性自转，教材 7.5.3 简化假设），
 *     密度 ρ 用分段指数模型 ρ = ρ0·exp(−(h−h0)/H)（教材式 7-5-5），
 *     节点间标高按教材式(7-5-9) H = Δh/ln(ρ_i/ρ_{i+1}) 由标准大气表均值定出
 *   - 太阳光压：a = −ν·(S/m)·C_R·p⊙·î⊙（教材式 7-6-1），
 *     p⊙ = 4.56×10⁻⁶ N/m²，C_R = 1 + η（教材式 7-6-2），
 *     ν 为圆柱形地影受晒系数（本影内 ν=0，无半影过渡——假设项）
 *
 * 积分器：定步长 RK4（经典四阶 Runge-Kutta），步长 h = clamp(T₀/120, 10, 60) s，
 * 全窗口状态落盘，帧间用三次 Hermite（p,v 端点）插值——位置插值误差 ~ O(Δt⁴)，
 * 对 60 s 采样间距的 LEO 轨道为厘米级。
 *
 * 单位约定：与 core 一致——距离 km、时间 s、μ = 398600.4418 km³/s²；
 * 阻力/光压的 SI 计算（kg/m³、m/s、N/m²）结果统一 ÷1000 化为 km/s²。
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) {
    module.exports = factory(require('./elements_core.js'));
  } else {
    root.ElementsPerturb = factory(root.ElementsCore);
  }
}(typeof globalThis !== 'undefined' ? globalThis : this, function (C) {
  'use strict';

  if (!C) { throw new Error('ElementsPerturb 需要 ElementsCore 先加载'); }

  var MU = C.MU;            // 398600.4418 km³/s²
  var RE = C.R_EARTH;       // 6378.137 km
  var TWO_PI = C.TWO_PI;

  /* 地球非球形带谐系数（EGM2008，未归一化 zonal 值）。
     J3/J4 用模块内可变变量存储，并暴露 getter/setter——
     node 自检需要临时置零做闭合 J2-only 对拍（教材式 7-4-41 严格只含 J2）。 */
  var J2 = 1.0826267e-3;
  var J3 = -2.5327e-6;
  var J4 = -1.6199e-6;

  var OMEGA_E = 7.2921159e-5;              // rad/s，地球自转角速度
  var P_SUN = 4.56e-6;                     // N/m²，教材式(7-6-1) 光压常数
  var OBLIQUITY = 23.4393 * C.D2R;         // 黄赤交角
  var YEAR_S = 365.25 * 86400;             // 太阳视运动圆轨道周期

  var WINDOW_S = 90 * 86400;               // 默认传播窗口：90 天

  /* ── 月球三体摄动模型（地球选项卡可选开关） ──────────────────────────
   *
   * 月球位置用圆化解析模型（非 SPICE 星历）：历元几何由 SPICE 探针
   * （tools/interactive_viewer/spice_probe.py，kernel = data/lunar/ephemeris/
   * de442s.bsp + latest_leapseconds.tls）在 2026-10-13 00:00 UTC 锁定：
   *   - 历元地月距离 396501.7 km（近拱点附近）、瞬时角速度 2.505983e-6 rad/s
   *   - 轨道面极轴/节线/历元相位收进 J2000（≈ICRF）单位矢量
   * 圆化假设（偏心率 0.0549、黄道倾角 5.14° 忽略）在 90 d 窗口内引入的
   * 月球位置误差约 ±21000 km（相位漂移主项），对 LEO/GEO 量级演示可接受。
   * 摄动加速度 = 直接项 + 间接项（教材式 13-2-4）。 */
  var MU_MOON = 4902.800066;               // km³/s²，DE440 月心引力常数
  var MOON_ORBIT = {
    a: 396501.7,                           // km，历元地心距（圆化半径）
    omega: 2.505983e-6,                    // rad/s，历元瞬时角速度
    pole: [-0.046868, -0.465408, 0.883855],   // 轨道面极轴（J2000 单位矢）
    node: [0.994968, -0.100197, 0.0],         // 升交点方向（J2000 单位矢）
    theta0: -2.233765,                     // rad，历元月球自节线起算方位角
  };
  (function buildMoonBasis() {
    var p = MOON_ORBIT.pole, n = MOON_ORBIT.node;
    /* 面内第二基矢 = pole × node（右手系） */
    MOON_ORBIT.eperp = [
      p[1] * n[2] - p[2] * n[1],
      p[2] * n[0] - p[0] * n[2],
      p[0] * n[1] - p[1] * n[0]
    ];
  })();
  /** 月球地心位置（ECI km），t 自历元起算（s）。圆化解析模型。 */
  function moonPosEci(tSec) {
    var th = MOON_ORBIT.theta0 + MOON_ORBIT.omega * tSec;
    var c = Math.cos(th), s = Math.sin(th);
    var a = MOON_ORBIT.a, n = MOON_ORBIT.node, e = MOON_ORBIT.eperp;
    return [a * (c * n[0] + s * e[0]), a * (c * n[1] + s * e[1]), a * (c * n[2] + s * e[2])];
  }

  /* 默认卫星假设参数（页面左栏「摄动力与假设参数」可改，均注明） */
  var DEFAULTS = {
    grav: true,        // 地球非球形引力（J2+J3+J4）
    drag: true,        // 大气阻力
    srp: true,         // 太阳光压
    moon3B: false,     // 月球第三体摄动（圆化解析模型，默认关）
    Cd: 2.2,           // 阻力系数（自由分子流，教材表 7-1 取值）
    CR: 1.5,           // 光压系数 = 1+η，η=0.5（教材式 7-6-2）
    SM: 0.02,          // 面质比 S/m，m²/kg（≈ 100 kg 卫星、2 m² 迎波面；光压与阻力共用）
    lam0: 0.0,         // 历元太阳黄经（rad），0 = 春分（假设项）
    windowS: WINDOW_S, // 传播窗口 s
    hStep: 0,          // 0 = 按 T₀/120 自动
  };

  /* 分段指数大气：标准大气表均值（h km, ρ kg/m³），标高由式(7-5-9)算好 */
  var ATMO_NODES = [
    [100, 5.25e-7], [150, 1.73e-9], [200, 2.41e-10], [250, 5.97e-11],
    [300, 1.87e-11], [350, 6.66e-12], [400, 2.62e-12], [450, 1.09e-12],
    [500, 4.76e-13], [550, 2.14e-13], [600, 9.89e-14], [650, 4.73e-14],
    [700, 2.36e-14], [750, 1.24e-14], [800, 6.95e-15], [850, 4.22e-15],
    [900, 2.78e-15], [950, 1.98e-15], [1000, 1.49e-15], [1250, 5.70e-16],
    [1500, 2.79e-16],
  ];
  var ATMO_TOP = 1500;   // 此高度以上密度取 0
  var atmoBands = [];
  (function buildAtmo() {
    for (var i = 0; i < ATMO_NODES.length - 1; i++) {
      var h0 = ATMO_NODES[i][0], r0 = ATMO_NODES[i][1];
      var h1 = ATMO_NODES[i + 1][0], r1 = ATMO_NODES[i + 1][1];
      var H = (h1 - h0) / Math.log(r0 / r1);   // 教材式(7-5-9)
      atmoBands.push([h0, r0, H]);
    }
  })();

  /** 大气密度 kg/m³。h<100 km 按首段向上外推（衰减演示会触地，仅保数值稳定）。 */
  function density(hKm) {
    if (hKm >= ATMO_TOP) { return 0; }
    if (hKm <= ATMO_NODES[0][0]) {
      var b0 = atmoBands[0];
      return b0[1] * Math.exp((b0[0] - hKm) / b0[2]);
    }
    for (var i = atmoBands.length - 1; i >= 0; i--) {
      if (hKm >= atmoBands[i][0]) {
        var b = atmoBands[i];
        return b[1] * Math.exp(-(hKm - b[0]) / b[2]);
      }
    }
    return 0;
  }

  /** 太阳方向单位矢量（ECI）：太阳视运动取匀速圆运动，黄经 λ = λ0 + 2πt/年。 */
  function sunDir(tSec, lam0) {
    var lam = (lam0 || 0) + TWO_PI * tSec / YEAR_S;
    var ce = Math.cos(OBLIQUITY), se = Math.sin(OBLIQUITY);
    return [Math.cos(lam), Math.sin(lam) * ce, Math.sin(lam) * se];
  }

  /** 圆柱形地影受晒系数：本影内 0，外部 1（无半影——假设项）。 */
  function shadowFactor(r, sHat) {
    var rdots = r[0] * sHat[0] + r[1] * sHat[1] + r[2] * sHat[2];
    if (rdots >= 0) { return 1; }
    var r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
    var d2 = r2 - rdots * rdots;      // 到日地轴的垂距平方
    return d2 < RE * RE ? 0 : 1;
  }

  /**
   * 带谐引力摄动加速度（J2+J3+J4，ECI）。
   * 位函数 U = μ/r·[1 − Σ J_n(Re/r)ⁿ P_n(s)]，s = z/r；
   * a = ∇U，对 r 与 s 的偏导分别按求和式展开（见文件头推导，J2 与教材闭式一致）。
   */
  function zonalAccel(r, out) {
    var x = r[0], y = r[1], z = r[2];
    var r2 = x * x + y * y + z * z;
    var rm = Math.sqrt(r2);
    var s = z / rm;
    /* Legendre 多项式及导数（递推，n≤4 展开，避免每步分配数组） */
    var P0 = 1, P1 = s;
    var P2 = (3 * s * P1 - P0) / 2;
    var P3 = (5 * s * P2 - 2 * P1) / 3;
    var P4 = (7 * s * P3 - 3 * P2) / 4;
    var P1d = 1;
    var P2d = (3 * (P1 + s * P1d) - P1d * 0) / 2;   // = 3s
    var P3d = (5 * (P2 + s * P2d) - 2 * P1d) / 3;   // = (15s²−3)/2
    var P4d = (7 * (P3 + s * P3d) - 3 * P2d) / 4;
    /* Σ J_n (Re/r)ⁿ P_n 及其对 r、s 的导数组合 */
    var Re_r = RE / rm;
    var t2 = Re_r * Re_r;      // (Re/r)²
    var t3 = t2 * Re_r;        // (Re/r)³
    var t4 = t3 * Re_r;        // (Re/r)⁴
    var sumP = J2 * t2 * P2 + J3 * t3 * P3 + J4 * t4 * P4;
    var sumR = 2 * J2 * t2 * P2 + 3 * J3 * t3 * P3 + 4 * J4 * t4 * P4; // Σ n·J_n(Re/r)ⁿ P_n
    var sumS = J2 * t2 * P2d + J3 * t3 * P3d + J4 * t4 * P4d;          // Σ J_n(Re/r)ⁿ P'_n
    var mu_r = MU / rm;
    var dUdr = -mu_r / rm * (1 - sumP - sumR);   // −μ/r²·[1 − Σ(1+n)J_n(Re/r)ⁿP_n]
    var dUds = -mu_r * sumS;                     // −μ/r·Σ J_n(Re/r)ⁿP'_n
    var r3 = r2 * rm;
    /* 全引力（中心 + 带谐摄动）= ∇[μ/r(1−Σ)]。
       注意 U 中谐项为 −μ/r·Σ，故 ∇U 的摄动部分与 Σ 的梯度反号；
       node 自检 FD 对拍确认下式即 ∇U_full。 */
    var ax = dUdr * (x / rm) + dUds * (-x * z / r3);
    var ay = dUdr * (y / rm) + dUds * (-y * z / r3);
    var az = dUdr * (z / rm) + dUds * ((x * x + y * y) / r3);
    out[0] = ax; out[1] = ay; out[2] = az;
    return out;
  }

  /**
   * 摄动加速度分解（km/s²）。y = [rx,ry,rz,vx,vy,vz]（km, km/s）。
   * 返回 { grav, zonal, drag, srp, moon3B, total }，每项为 [ax,ay,az]。
   * 阻力/光压在 SI 下计算后 ÷1000 转 km/s²。
   */
  function accelBreakdown(y, tSec, o) {
    o = o || DEFAULTS;
    var r = [y[0], y[1], y[2]], v = [y[3], y[4], y[5]];
    var rm = Math.hypot(r[0], r[1], r[2]);
    var grav = [
      -MU * r[0] / (rm * rm * rm),
      -MU * r[1] / (rm * rm * rm),
      -MU * r[2] / (rm * rm * rm),
    ];
    var zonal = [0, 0, 0];
    if (o.grav) {
      /* zonalAccel 返回全引力（中心 + 带谐），减去中心引力得纯摄动项 */
      zonalAccel(r, zonal);
      zonal[0] -= grav[0]; zonal[1] -= grav[1]; zonal[2] -= grav[2];
    }
    var drag = [0, 0, 0];
    if (o.drag) {
      var rho = density(rm - RE);
      if (rho > 0) {
        /* 大气随地球刚性自转：v_rel = v − ω⊕×r，ω⊕ = (0,0,ω) */
        var wxr = [-OMEGA_E * y[1], OMEGA_E * y[0], 0];
        var vrx = (v[0] - wxr[0]) * 1000, vry = (v[1] - wxr[1]) * 1000, vrz = (v[2] - wxr[2]) * 1000;
        var vrm = Math.hypot(vrx, vry, vrz);
        var k = -0.5 * o.Cd * o.SM * rho * vrm / 1000;   // km/s² per (m/s) 速度分量
        drag = [k * vrx, k * vry, k * vrz];
      }
    }
    var srp = [0, 0, 0];
    if (o.srp) {
      var sHat = sunDir(tSec, o.lam0);
      var nu = shadowFactor(r, sHat);
      if (nu > 0) {
        var ks = -nu * o.CR * o.SM * P_SUN / 1000;     // km/s²
        srp = [ks * sHat[0], ks * sHat[1], ks * sHat[2]];
      }
    }
    var moon3B = [0, 0, 0];
    if (o.moon3B) {
      /* 月球第三体：a = −μ☾[(r−r☾)/d³ + r☾/|r☾|³]（直接项+间接项） */
      var rmn = moonPosEci(tSec);
      var dx = r[0] - rmn[0], dy = r[1] - rmn[1], dz = r[2] - rmn[2];
      var d2m = dx * dx + dy * dy + dz * dz;
      var dm = Math.sqrt(d2m);
      var rn2 = rmn[0] * rmn[0] + rmn[1] * rmn[1] + rmn[2] * rmn[2];
      var rn = Math.sqrt(rn2);
      var k1 = -MU_MOON / (d2m * dm);
      var k2 = -MU_MOON / (rn2 * rn);
      moon3B = [k1 * dx + k2 * rmn[0], k1 * dy + k2 * rmn[1], k1 * dz + k2 * rmn[2]];
    }
    var total = [0, 0, 0];
    for (var i = 0; i < 3; i++) {
      total[i] = grav[i] + zonal[i] + drag[i] + srp[i] + moon3B[i];
    }
    return { grav: grav, zonal: zonal, drag: drag, srp: srp, moon3B: moon3B, total: total };
  }

  /** 动力学右端 f(y)（仅加速度部分；自引力 + 开关摄动）。 */
  function deriv(y, tSec, o, out) {
    var b = accelBreakdown(y, tSec, o);
    out[0] = y[3]; out[1] = y[4]; out[2] = y[5];
    out[3] = b.total[0]; out[4] = b.total[1]; out[5] = b.total[2];
    return out;
  }

  /**
   * 定步长 RK4 全程积分。y0 六个分量；返回 {n, h, t, y}（每步落盘，含起点）。
   */
  function propagate(y0, windowS, h, o) {
    var n = Math.max(1, Math.ceil(windowS / h));
    var t = new Float64Array(n + 1);
    var y = new Float64Array((n + 1) * 6);
    var yc = Float64Array.from(y0);
    y.set(yc, 0); t[0] = 0;
    var k1 = new Float64Array(6), k2 = new Float64Array(6),
        k3 = new Float64Array(6), k4 = new Float64Array(6),
        yt = new Float64Array(6);
    for (var i = 0; i < n; i++) {
      var ti = i * h;
      deriv(yc, ti, o, k1);
      for (var j = 0; j < 6; j++) { yt[j] = yc[j] + 0.5 * h * k1[j]; }
      deriv(yt, ti + 0.5 * h, o, k2);
      for (var j = 0; j < 6; j++) { yt[j] = yc[j] + 0.5 * h * k2[j]; }
      deriv(yt, ti + 0.5 * h, o, k3);
      for (var j = 0; j < 6; j++) { yt[j] = yc[j] + h * k3[j]; }
      deriv(yt, ti + h, o, k4);
      for (var j = 0; j < 6; j++) {
        yc[j] += h / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
      }
      t[i + 1] = (i + 1) * h;
      y.set(yc, (i + 1) * 6);
    }
    return { n: n, h: h, t: t, y: y };
  }

  /* ── 笛卡尔 ↔ 经典根数 ─────────────────────────────────────────── */

  function clamp1(x) { return Math.min(1, Math.max(-1, x)); }
  function dot3(a, b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
  function cross3(a, b) {
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  }

  /**
   * ECI 笛卡尔状态 → 经典密切根数（弧度）。
   * 奇异情形按工程惯例退化处理：圆轨道 e<1e-9 时 ω=0、ν 取幅角；
   * 赤道轨道 i<~1e-7 时 Ω=0、ω 自 x 轴量起。
   */
  function cartToClassical(r, v) {
    var rm = Math.hypot(r[0], r[1], r[2]);
    var vm2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    var h = cross3(r, v);
    var hm = Math.hypot(h[0], h[1], h[2]);
    var En = vm2 / 2 - MU / rm;
    var a = -MU / (2 * En);
    var i = Math.acos(clamp1(h[2] / hm));
    var n = [-h[1], h[0], 0];
    var nm = Math.hypot(n[0], n[1]);
    var Om = nm > 1e-12 ? Math.atan2(n[1], n[0]) : 0;
    var vh = cross3(v, h);
    var eVec = [vh[0] / MU - r[0] / rm, vh[1] / MU - r[1] / rm, vh[2] / MU - r[2] / rm];
    var e = Math.hypot(eVec[0], eVec[1], eVec[2]);
    var hunit = [h[0] / hm, h[1] / hm, h[2] / hm];
    var angAbout = function (u, w) {
      var cuw = cross3(u, w);
      return Math.atan2(dot3(cuw, hunit), dot3(u, w));
    };
    var w, nu;
    if (e > 1e-9) {
      if (nm > 1e-12) { w = angAbout(n, eVec); }
      else { w = Math.atan2(eVec[1], eVec[0]); }   // 赤道：ω 自春分点方向量起
      nu = angAbout(eVec, r);
    } else {
      w = 0;
      nu = nm > 1e-12 ? angAbout(n, r) : Math.atan2(r[1], r[0]);
    }
    return { a: a, e: e, i: i, Om: C.wrap2pi(Om), w: C.wrap2pi(w), nu: C.wrap2pi(nu) };
  }

  /* ── 数值摄动传播器（契约同 KeplerPropagator，见 core 注释） ─────── */

  /**
   * constructor(el, opts)：
   *   el   = {a, e, i, Om, w}（a km；角度弧度；积分历元恒锚近地点 ν=0，el.nu 忽略）
   *   opts = DEFAULTS 的覆盖子集（开关、Cd、CR、SM、lam0、windowS、hStep）
   *
   * 构造时同步完成全窗口 RK4 积分（LEO 约 1–3 s）；页面请先显示
   * 「传播计算中」再实例化。之后任意时刻状态走三次 Hermite 插值。
   */
  function PerturbedPropagator(el, opts) {
    this.isPerturbed = true;
    this.o = {};
    for (var k in DEFAULTS) { this.o[k] = DEFAULTS[k]; }
    if (opts) { for (var k2 in opts) { this.o[k2] = opts[k2]; } }

    this.el0 = { a: el.a, e: el.e, i: el.i, Om: el.Om, w: el.w };
    /* 积分历元恒锚近地点（ν=0）：内部 t=0 ↔ 近地点，与 timeAtTrueAnomaly
       及页面 tSincePeri（自近地点起算）同基准。el.nu 仅由页面 rebuild 换算成
       tSincePeri 初值，不参与积分初值——否则 y0@ν₀ 与近地点基准错位，
       rebuild（拖 Ω/ω/a/e/i）后相位跳变。 */
    this.nu0 = 0;

    var st0 = C.classicalToState(this.el0, 0);
    this.y0 = [st0.rE[0], st0.rE[1], st0.rE[2], st0.vE[0], st0.vE[1], st0.vE[2]];

    var T0 = C.period(this.el0.a);
    this.h = this.o.hStep > 0 ? this.o.hStep : Math.min(60, Math.max(10, T0 / 120));
    var sol = propagate(this.y0, this.o.windowS, this.h, this.o);
    this.nSteps = sol.n; this.h = sol.h;
    this.tArr = sol.t; this.yArr = sol.y;
    this.window = this.nSteps * this.h;
    this.tNow = 0;
    this._osc = cartToClassical(this.y0.slice(0, 3), this.y0.slice(3, 6));
  }

  /** 三次 Hermite 插值：tSec 处状态 → [rx,ry,rz,vx,vy,vz]（不做边界外推）。 */
  PerturbedPropagator.prototype.stateVecAt = function (tSec) {
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

  /** tSec 处完整状态对象（与 classicalToState 同构的常用字段）。 */
  PerturbedPropagator.prototype.stateAtTime = function (tSec) {
    var y = this.stateVecAt(tSec);
    var rE = [y[0], y[1], y[2]], vE = [y[3], y[4], y[5]];
    var osc = cartToClassical(rE, vE);
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
      rE: rE, vE: vE, r: rm, v: vm, alt: rm - RE,
      vr: vr, vt: vt, gamma: Math.atan2(vr, vt),
      osc: osc, E: E, M: E - osc.e * Math.sin(E), tSec: tSec,
    };
  };

  /** 契约：t（s，历元起算）→ 瞬时密切真近点角；同时把 tNow 锚到 t。 */
  PerturbedPropagator.prototype.trueAnomalyAtTime = function (tSec) {
    this.tNow = Math.min(Math.max(tSec, 0), this.window);
    this._osc = this.stateAtTime(this.tNow).osc;
    return this._osc.nu;
  };

  /** 契约：ν → 历元密切根数下的等价二体时间（相位锚定用）。 */
  PerturbedPropagator.prototype.timeAtTrueAnomaly = function (nuRad) {
    var M = C.nuToM(nuRad, this.el0.e).M;
    return M / C.meanMotion(this.el0.a);
  };

  /**
   * 契约：返回 tNow 处的积分状态（摄动轨道 ν 不唯一决定状态，
   * 故忽略 nu 参数；页面卫星标记/面板必须走本方法）。
   */
  PerturbedPropagator.prototype.stateAtTrueAnomaly = function (nuRad) {
    return this.stateAtTime(this.tNow);
  };

  /** 契约：当前密切轨道周期（s）——播放节奏与相位显示用。 */
  PerturbedPropagator.prototype.periodSeconds = function () {
    return C.period(this._osc.a);
  };

  /** tNow 处的密切根数 {a,e,i,Om,w,nu}（弧度）——瞬时椭圆绘制用。 */
  PerturbedPropagator.prototype.osculatingElements = function () {
    return this._osc;
  };

  PerturbedPropagator.prototype.windowSeconds = function () { return this.window; };
  PerturbedPropagator.prototype.stepSeconds = function () { return this.h; };

  /** 当前 tNow 处各摄动加速度大小（km/s²）。 */
  PerturbedPropagator.prototype.forceMagnitudes = function () {
    var y = this.stateVecAt(this.tNow);
    var b = accelBreakdown(y, this.tNow, this.o);
    var mag = function (a) { return Math.hypot(a[0], a[1], a[2]); };
    return { grav: mag(b.grav), zonal: mag(b.zonal), drag: mag(b.drag), srp: mag(b.srp), moon3B: mag(b.moon3B) };
  };

  /**
   * 全窗口拖尾点列（ECI km，降采样到 ≤ maxPts 点）。
   * 返回 {x,y,z} 三个普通数组（three.js BufferAttribute 直接可用）。
   */
  PerturbedPropagator.prototype.trailPoints = function (maxPts) {
    var stride = Math.max(1, Math.ceil((this.nSteps + 1) / (maxPts || 3000)));
    var x = [], y = [], z = [];
    for (var i = 0; i <= this.nSteps; i += stride) {
      var o = i * 6;
      x.push(this.yArr[o]); y.push(this.yArr[o + 1]); z.push(this.yArr[o + 2]);
    }
    var last = this.nSteps * 6;
    var lx = this.yArr[last], ly = this.yArr[last + 1], lz = this.yArr[last + 2];
    if (x[x.length - 1] !== lx) { x.push(lx); y.push(ly); z.push(lz); }
    return { x: x, y: y, z: z, stride: stride };
  };

  PerturbedPropagator.id = 'perturbed';
  PerturbedPropagator.label = '数值摄动（J2–J4 + 大气阻力 + 太阳光压，90 d）';
  PerturbedPropagator.hasPerturbation = true;

  /* 自动注册：elements_demo.html 传播模型下拉框零改动接入 */
  C.registerPropagator(PerturbedPropagator);

  return {
    J2: J2,
    get J3() { return J3; }, set J3(v) { J3 = v; },
    get J4() { return J4; }, set J4(v) { J4 = v; },
    OMEGA_E: OMEGA_E, P_SUN: P_SUN, OBLIQUITY: OBLIQUITY, YEAR_S: YEAR_S,
    WINDOW_S: WINDOW_S, DEFAULTS: DEFAULTS,
    ATMO_NODES: ATMO_NODES,
    MU_MOON: MU_MOON, MOON_ORBIT: MOON_ORBIT,
    density: density,
    sunDir: sunDir,
    shadowFactor: shadowFactor,
    zonalAccel: zonalAccel,
    accelBreakdown: accelBreakdown,
    moonPosEci: moonPosEci,
    propagate: propagate,
    cartToClassical: cartToClassical,
    PerturbedPropagator: PerturbedPropagator,
  };
}));
