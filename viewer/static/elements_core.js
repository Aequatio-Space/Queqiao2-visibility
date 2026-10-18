/*!
 * elements_core.js — 轨道六根数二体核心数学（浏览器 / node 双端复用）。
 *
 * 与 tools/interactive_viewer/gen_elements_selftest.py 中的手写 numpy 路径逐公式同源
 * （Vallado《Fundamentals of Astrodynamics》经典根数 → PQW → ECI 链路），
 * 并由 poliastro 0.17 Orbit.from_classical 交叉验证（见 static/elements_selftest.js）。
 *
 * 约定：
 *   - 角度在 core 内部一律用弧度；UI 层负责度 ↔ 弧度换算。
 *   - el = {a, e, i, Om, w}：a 单位 km，e 无量纲，i/Om/w 弧度。
 *   - 坐标系：ECI 惯性系（x → 春分点 ♈，z → 北极），单位 km。
 *   - μ = 398600.4418 km³/s²（EGM2008，与 poliastro/astropy Earth GM 一致）
 *   - R⊕ = 6378.137 km（WGS-84 半长轴，球体近似）
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) { module.exports = factory(); }
  else { root.ElementsCore = factory(); }
}(typeof globalThis !== 'undefined' ? globalThis : this, function () {
  'use strict';

  var MU = 398600.4418;    // km^3/s^2
  var R_EARTH = 6378.137;  // km
  var D2R = Math.PI / 180;
  var R2D = 180 / Math.PI;
  var TWO_PI = Math.PI * 2;

  function wrap2pi(x) {
    x = x % TWO_PI;
    if (x < 0) { x += TWO_PI; }
    return x;
  }

  /** 两个角量（弧度）之间的最短角距离。 */
  function angDist(a, b) {
    var d = Math.abs(wrap2pi(a) - wrap2pi(b));
    return Math.min(d, TWO_PI - d);
  }

  /**
   * PQW → ECI 旋转矩阵 R = R₃(−Ω)·R₁(−i)·R₃(−ω)（行主序 9 元素数组）。
   * 即经典根数到惯性系的标准变换（Vallado Eq. 3-29）。
   */
  function pqwToEciMatrix(iRad, OmRad, wRad) {
    var cO = Math.cos(OmRad), sO = Math.sin(OmRad);
    var ci = Math.cos(iRad),  si = Math.sin(iRad);
    var cw = Math.cos(wRad), sw = Math.sin(wRad);
    return [
      cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci,  sO * si,
      sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si,
      si * sw,                 si * cw,                 ci
    ];
  }

  function matVec(m, v) {
    return [
      m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
      m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
      m[6] * v[0] + m[7] * v[1] + m[8] * v[2]
    ];
  }

  function norm3(v) { return Math.hypot(v[0], v[1], v[2]); }

  /** 真近点角 → 偏近点角（半角公式，全象限安全）。 */
  function nuToE(nuRad, e) {
    return 2 * Math.atan2(Math.sqrt(1 - e) * Math.sin(nuRad / 2),
                          Math.sqrt(1 + e) * Math.cos(nuRad / 2));
  }

  /** 偏近点角 → 真近点角。 */
  function eToNu(ERad, e) {
    return 2 * Math.atan2(Math.sqrt(1 + e) * Math.sin(ERad / 2),
                          Math.sqrt(1 - e) * Math.cos(ERad / 2));
  }

  /**
   * 开普勒方程 M = E − e·sinE 的牛顿迭代求解。
   * 返回 {E, iters}；E 已归约到 [0, 2π)。
   * 初值策略：e < 0.8 取 E₀ = M，否则取 E₀ = π（Vallado 建议，e ≤ 0.98 收敛可靠）。
   */
  function mToE(MRad, e) {
    var M = wrap2pi(MRad);
    if (e < 1e-12) { return { E: M, iters: 0 }; }
    var E = (e < 0.8) ? M : Math.PI;
    var it = 0;
    for (; it < 100; it++) {
      var f = E - e * Math.sin(E) - M;
      var fp = 1 - e * Math.cos(E);
      var d = f / fp;
      E -= d;
      if (Math.abs(d) < 1e-13) { it += 1; break; }
    }
    return { E: wrap2pi(E), iters: it };
  }

  /** 真近点角 → 平近点角（附带 E）。 */
  function nuToM(nuRad, e) {
    var E = nuToE(nuRad, e);
    return { M: E - e * Math.sin(E), E: E };
  }

  /*
   * 经典六根数（ν 为自变量）→ 完整二体状态。
   * 可选 mu / bodyR：中心天体引力常数（km³/s²）与参考半径（km），
   * 缺省为地球（MU / R_EARTH）。月球传 MU_MOON / R_MOON
   * （定义在 elements_cr3bp.js，页面层传入），核心不硬编码月球常数。
   * 返回：{p, r, alt, rP, vP, matrix, rE, vE, v, vr, vt, gamma, E, M, hmag}
   *   - rE/vE：惯性系位置/速度（km, km/s）
   *   - vr/vt：径向/横向速度分量（km/s）
   *   - gamma：飞行路径角（速度与当地水平面的夹角，弧度，远离中心天体为正）
   *   - E/M：该 ν 对应的偏近点角/平近点角（弧度，M 已归约）
   */
  function classicalToState(el, nuRad, mu, bodyR) {
    var MU_ = mu || MU;
    var R_ = bodyR || R_EARTH;
    var a = el.a, e = el.e;
    var p = a * (1 - e * e);
    var r = p / (1 + e * Math.cos(nuRad));
    var rP = [r * Math.cos(nuRad), r * Math.sin(nuRad), 0];
    var smp = Math.sqrt(MU_ / p);
    var vP = [-smp * Math.sin(nuRad), smp * (e + Math.cos(nuRad)), 0];
    var R = pqwToEciMatrix(el.i, el.Om, el.w);
    var rE = matVec(R, rP);
    var vE = matVec(R, vP);
    var vr = smp * e * Math.sin(nuRad);
    var vt = smp * (1 + e * Math.cos(nuRad));
    var gamma = Math.atan2(e * Math.sin(nuRad), 1 + e * Math.cos(nuRad));
    var E = nuToE(nuRad, e);
    var M = E - e * Math.sin(E);
    return {
      p: p, r: r, alt: r - R_,
      rP: rP, vP: vP, matrix: R, rE: rE, vE: vE,
      v: norm3(vE), vr: vr, vt: vt, gamma: gamma,
      E: E, M: M, hmag: Math.sqrt(MU_ * p)
    };
  }

  /** 平均角速度 n = √(μ/a³)（rad/s）。mu 缺省地球。 */
  function meanMotion(a, mu) { return Math.sqrt((mu || MU) / (a * a * a)); }

  /** 轨道周期 T = 2π/n（s）。mu 缺省地球。 */
  function period(a, mu) { return TWO_PI / meanMotion(a, mu); }

  /** 自最近过近拱点起算的时间 t(ν)（s）。mu 缺省地球。 */
  function timeAtTrueAnomaly(el, nuRad, mu) {
    var E = nuToE(nuRad, el.e);
    return (E - el.e * Math.sin(E)) / meanMotion(el.a, mu);
  }

  /** t → ν（s，经 M → 开普勒方程 → ν 链路）。mu 缺省地球。 */
  function trueAnomalyAtTime(el, tSec, mu) {
    var M = wrap2pi(meanMotion(el.a, mu) * tSec);
    return eToNu(mToE(M, el.e).E, el.e);
  }

  /* ── 传播器（Propagator）契约 — 摄动接口预留 ────────────────────────────
   *
   * 当前仅注册开普勒二体。未来引入轨道摄动（J2、日月三体、大气阻力、数值外推…）时，
   * 实现同一契约并调用 registerPropagator(类) 即可被 UI 下拉框自动识别，页面与
   * 3D 场景零改动：
   *
   *   static id                — 唯一标识（<select> 的 value）
   *   static label             — 下拉框显示名
   *   static hasPerturbation   — 是否含摄动（UI 据此提示“含摄动修正”）
   *   constructor(elements)    — elements: {a, e, i, Om, w}（角度弧度）
   *   stateAtTrueAnomaly(nuRad) → 与 classicalToState 同构的状态对象
   *                              （至少含 rE/vE/r/v/alt/gamma；摄动模型应在此叠加 Δr/Δv）
   *   trueAnomalyAtTime(tSec)  → nuRad（tSec 自最近一次过近地点起算）
   *   timeAtTrueAnomaly(nuRad) → tSec
   *   periodSeconds()          → 轨道周期（播放动画用）
   */
  function KeplerPropagator(el, mu, bodyR) {
    this.el = { a: el.a, e: el.e, i: el.i, Om: el.Om, w: el.w };
    this.mu = mu || MU;          // 中心天体引力常数（缺省地球）
    this.bodyR = bodyR || R_EARTH;
  }
  KeplerPropagator.id = 'kepler';
  KeplerPropagator.label = '开普勒二体（无摄动）';
  KeplerPropagator.hasPerturbation = false;
  KeplerPropagator.prototype.stateAtTrueAnomaly = function (nuRad) {
    return classicalToState(this.el, nuRad, this.mu, this.bodyR);
  };
  KeplerPropagator.prototype.trueAnomalyAtTime = function (tSec) {
    return trueAnomalyAtTime(this.el, tSec, this.mu);
  };
  KeplerPropagator.prototype.timeAtTrueAnomaly = function (nuRad) {
    return timeAtTrueAnomaly(this.el, nuRad, this.mu);
  };
  KeplerPropagator.prototype.periodSeconds = function () { return period(this.el.a, this.mu); };

  var registry = {};
  function registerPropagator(cls) {
    if (!cls || !cls.id || !cls.label) {
      throw new Error('Propagator 契约不符：需要 static id / label');
    }
    registry[cls.id] = cls;
    return cls;
  }
  registerPropagator(KeplerPropagator);

  return {
    MU: MU, R_EARTH: R_EARTH, D2R: D2R, R2D: R2D, TWO_PI: TWO_PI,
    wrap2pi: wrap2pi, angDist: angDist,
    pqwToEciMatrix: pqwToEciMatrix, matVec: matVec, norm3: norm3,
    nuToE: nuToE, eToNu: eToNu, mToE: mToE, nuToM: nuToM,
    classicalToState: classicalToState,
    meanMotion: meanMotion, period: period,
    timeAtTrueAnomaly: timeAtTrueAnomaly,
    trueAnomalyAtTime: trueAnomalyAtTime,
    KeplerPropagator: KeplerPropagator,
    registerPropagator: registerPropagator,
    propagators: registry
  };
}));
