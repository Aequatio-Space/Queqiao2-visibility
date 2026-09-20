/* coast_tool/app.js - 月球转移发射窗口计算器（停泊轨道滑行方案）
 *
 * 由 coast_tool/index.html 的 4 段内联脚本按原顺序提取（脚本内容未改动），
 * 以满足全站 CSP `script-src 'self'`（禁止内联脚本执行）。
 * 来源：lunarlaunchwindow 项目 coast_tool/
 * （MIT License, Copyright (c) 2026 Zhuoxiao）。
 */

/* ===== 原 index.html 内联脚本块 1/4（未改动） ===== */
/* =====================================================================
 *  第一部分: 天文历表与轨道力学  (equatorial-of-date 惯性系)
 * ===================================================================== */
const DEG = Math.PI / 180;
const MU = 398600.4418;      // km^3/s^2
const MU_MOON = 4902.800;
const R_EARTH = 6378.137;    // km equatorial
const R_MEAN  = 6371.0;

const LAUNCH_SITES = {
    wenchang: { name: '文昌', lat: 19.6145, lon: 110.951, tz: 8, tzName: 'UTC+8',
                azMin: 90, azMax: 110 },
    xichang:  { name: '西昌', lat: 28.2466, lon: 102.027, tz: 8, tzName: 'UTC+8',
                azMin: 95, azMax: 105 },
    ksc:      { name: 'KSC',  lat: 28.5729, lon: -80.649, tz: -5, tzName: 'UTC-5',
                azMin: 37, azMax: 114 }
};

function mod360(x) { return ((x % 360) + 360) % 360; }

function julianDateUTC(ms) {           // ms = Date.getTime()
    return ms / 86400000 + 2440587.5;
}
function jdTT(jdUTC) { return jdUTC + 69.2 / 86400; }  // ΔT ~ 69 s (2020s)

function gmstDeg(jdUTC) {
    const d = jdUTC - 2451545.0, T = d / 36525;
    return mod360(280.46061837 + 360.98564736629 * d + 0.000387933 * T * T - T * T * T / 38710000);
}
function obliquityDeg(jdTTv) {
    const T = (jdTTv - 2451545.0) / 36525;
    return 23.43929111 - 0.01300417 * T - 1.64e-7 * T * T;
}

/* ---- 月球位置: Meeus ch.47 截断级数 (~0.01° 经度, ~0.5 km 距离主项) ---- */
// [D, M, M', F, Σl(1e-6 deg), Σr(1e-3 km)]
const LR_TERMS = [
[0,0,1,0, 6288774,-20905355],[2,0,-1,0, 1274027,-3699111],[2,0,0,0, 658314,-2955968],
[0,0,2,0, 213618,-569925],[0,1,0,0, -185116,48888],[0,0,0,2, -114332,-3149],
[2,0,-2,0, 58793,246158],[2,-1,-1,0, 57066,-152138],[2,0,1,0, 53322,-170733],
[2,-1,0,0, 45758,-204586],[0,1,-1,0, -40923,-129620],[1,0,0,0, -34720,108743],
[0,1,1,0, -30383,104755],[2,0,0,-2, 15327,10321],[0,0,1,2, -12528,0],
[0,0,1,-2, 10980,79661],[4,0,-1,0, 10675,-34782],[0,0,3,0, 10034,-23210],
[4,0,-2,0, 8548,-21636],[2,1,-1,0, -7888,24208],[2,1,0,0, -6766,30824],
[1,0,-1,0, -5163,-8379],[1,1,0,0, 4987,-16675],[2,-1,1,0, 4036,-12831],
[2,0,2,0, 3994,-10445],[4,0,0,0, 3861,-11650],[2,0,-3,0, 3665,14403],
[0,1,-2,0, -2689,-7003],[2,0,-1,2, -2602,0],[2,-1,-2,0, 2390,10056],
[1,0,1,0, -2348,6322],[2,-2,0,0, 2236,-9884]
];
// [D, M, M', F, Σb(1e-6 deg)]
const B_TERMS = [
[0,0,0,1, 5128122],[0,0,1,1, 280602],[0,0,1,-1, 277693],[2,0,0,-1, 173237],
[2,0,-1,1, 55413],[2,0,-1,-1, 46271],[2,0,0,1, 32573],[0,0,2,1, 17198],
[2,0,1,-1, 9266],[0,0,2,-1, 8822],[2,-1,0,-1, 8216],[2,0,-2,-1, 4324],
[2,0,1,1, 4200],[2,1,0,-1, -3359],[2,-1,-1,1, 2463],[2,-1,0,1, 2211],
[2,-1,-1,-1, 2065],[0,1,-1,-1, -1870],[4,0,-1,-1, 1828],[0,1,0,1, -1794]
];

// 返回黄道 of date: {lonDeg, latDeg, distKm}
function moonEcliptic(jdUTC) {
    const T = (jdTT(jdUTC) - 2451545.0) / 36525;
    const Lp = mod360(218.3164477 + 481267.88123421*T - 0.0015786*T*T + T*T*T/538841);
    const D  = mod360(297.8501921 + 445267.1114034*T - 0.0018819*T*T + T*T*T/545868);
    const M  = mod360(357.5291092 + 35999.0502909*T - 0.0001536*T*T);
    const Mp = mod360(134.9633964 + 477198.8675055*T + 0.0087414*T*T + T*T*T/69699);
    const F  = mod360(93.2720950 + 483202.0175233*T - 0.0036539*T*T);
    const E  = 1 - 0.002516*T - 0.0000074*T*T;
    const A1 = mod360(119.75 + 131.849*T);
    const A2 = mod360(53.09 + 479264.290*T);
    const A3 = mod360(313.45 + 481266.484*T);

    let sl = 0, sr = 0, sb = 0;
    for (const [d, m, mp, f, cl, cr] of LR_TERMS) {
        const arg = (d*D + m*M + mp*Mp + f*F) * DEG;
        const ef = (m === 1 || m === -1) ? E : (m === 2 || m === -2) ? E*E : 1;
        sl += cl * ef * Math.sin(arg);
        sr += cr * ef * Math.cos(arg);
    }
    for (const [d, m, mp, f, cb] of B_TERMS) {
        const arg = (d*D + m*M + mp*Mp + f*F) * DEG;
        const ef = (m === 1 || m === -1) ? E : (m === 2 || m === -2) ? E*E : 1;
        sb += cb * ef * Math.sin(arg);
    }
    sl += 3958*Math.sin(A1*DEG) + 1962*Math.sin((Lp - F)*DEG) + 318*Math.sin(A2*DEG);
    sb += -2235*Math.sin(Lp*DEG) + 382*Math.sin(A3*DEG)
        + 175*Math.sin((A1 - F)*DEG) + 175*Math.sin((A1 + F)*DEG)
        + 127*Math.sin((Lp - Mp)*DEG) - 115*Math.sin((Lp + Mp)*DEG);

    return {
        lonDeg: mod360(Lp + sl / 1e6),
        latDeg: sb / 1e6,
        distKm: 385000.56 + sr / 1000
    };
}

function sunEclipticLon(jdUTC) {
    const n = jdTT(jdUTC) - 2451545.0;
    const L = mod360(280.460 + 0.9856474 * n);
    const g = mod360(357.528 + 0.9856003 * n);
    return mod360(L + 1.915 * Math.sin(g*DEG) + 0.020 * Math.sin(2*g*DEG));
}

/* ---- 向量工具 ---- */
const V = {
    dot: (a,b) => a[0]*b[0]+a[1]*b[1]+a[2]*b[2],
    cross: (a,b) => [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]],
    norm: a => Math.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]),
    unit: a => { const n = V.norm(a); return [a[0]/n, a[1]/n, a[2]/n]; },
    scale: (a,s) => [a[0]*s, a[1]*s, a[2]*s],
    add: (a,b) => [a[0]+b[0], a[1]+b[1], a[2]+b[2]],
    sub: (a,b) => [a[0]-b[0], a[1]-b[1], a[2]-b[2]]
};

// 黄道(of date) → 赤道(of date)
function eclToEq(v, jdUTC) {
    const e = obliquityDeg(jdTT(jdUTC)) * DEG, c = Math.cos(e), s = Math.sin(e);
    return [v[0], v[1]*c - v[2]*s, v[1]*s + v[2]*c];
}
function moonEq(jdUTC) {
    const m = moonEcliptic(jdUTC);
    const lr = m.lonDeg*DEG, br = m.latDeg*DEG;
    const ec = [Math.cos(br)*Math.cos(lr), Math.cos(br)*Math.sin(lr), Math.sin(br)];
    const eq = eclToEq(ec, jdUTC);
    return { u: eq, dist: m.distKm, ecl: m };
}
// 月球轨道面法向 (数值差分)
function moonPlaneNormal(jdUTC) {
    const dt = 0.04; // days
    const a = moonEq(jdUTC - dt), b = moonEq(jdUTC + dt);
    const n = V.unit(V.cross(V.scale(a.u, a.dist), V.scale(b.u, b.dist)));
    return n[2] < 0 ? V.scale(n, -1) : n;
}

// 发射场地心方向 (惯性系, 赤道 of date)
function sitePosUnit(site, jdUTC) {
    // 地心纬度 (WGS84 修正)
    const f = 1/298.257223563;
    const phi_gc = Math.atan((1-f)*(1-f) * Math.tan(site.lat*DEG));
    const ra = (gmstDeg(jdUTC) + site.lon) * DEG;
    return { u: [Math.cos(phi_gc)*Math.cos(ra), Math.cos(phi_gc)*Math.sin(ra), Math.sin(phi_gc)],
             phi_gc };
}

/* ---- Kepler: 近地点(ν=0)到 ν 的飞行时间 ---- */
function kepTOF(rp, e, nuDeg) {
    const a = rp / (1 - e);
    const nu = nuDeg * DEG;
    let E = 2 * Math.atan2(Math.sqrt(1-e) * Math.sin(nu/2), Math.sqrt(1+e) * Math.cos(nu/2));
    let Mm = E - e * Math.sin(E);
    if (Mm < 0) Mm += 2 * Math.PI;
    return Mm / Math.sqrt(MU / (a*a*a));  // s
}
// 反解: 从 TLI 起 t 秒后的真近点角 (deg)
function kepNuAt(rp, e, tSec) {
    const a = rp / (1 - e);
    const n = Math.sqrt(MU / (a*a*a));
    const Mm = n * tSec;
    // 高偏心率近近地点用立方根初值, 否则用 M
    let E = (e > 0.8 && Mm < Math.PI) ? Math.cbrt(6 * Mm) : Mm;
    for (let i = 0; i < 40; i++) {
        const dE = (E - e*Math.sin(E) - Mm) / (1 - e*Math.cos(E));
        E -= dE;
        if (Math.abs(dE) < 1e-10) break;
    }
    const nu = 2 * Math.atan2(Math.sqrt(1+e) * Math.sin(E/2), Math.sqrt(1-e) * Math.cos(E/2));
    return mod360(nu / DEG);
}

/* =====================================================================
 *  第二部分: 发射窗口求解 (核心)
 *
 *  给定发射时刻 tL:
 *    LTO 平面 = span{ r_site(tL), r_moon(tA) }  (tA = tL + 转移时间)
 *    → 唯一确定 倾角 i、发射方位角 Az(顺行切向)、升交点 Ω
 *  滑行闭合: TLI 在近地点、切向点火;
 *    航程角闭合  θ_total = θ_ascent + ω_park·t_coast + Δν
 *    时间闭合    T_transfer = T_ascent + t_coast + TOF_kepler(Δν)
 *    两式联立 → 解 t_coast (扫描找根)
 * ===================================================================== */
function solveGeometry(jdL, site, transferDays) {
    const jdA = jdL + transferDays;
    const sp = sitePosUnit(site, jdL);
    const mn = moonEq(jdA);
    const u1 = sp.u, u2 = mn.u;

    let n = V.cross(u1, u2);
    let theta = Math.atan2(V.norm(n), V.dot(u1, u2)) / DEG;
    n = V.unit(n);
    if (n[2] < 0) { n = V.scale(n, -1); theta = 360 - theta; }   // 顺行

    const inc = Math.acos(Math.max(-1, Math.min(1, n[2]))) / DEG;
    const raan = mod360(Math.atan2(n[0], -n[1]) / DEG);

    // 唯一方位角: 速度方向 v = n × u1, 在当地地平系分解
    const v = V.cross(n, u1);
    const east = V.unit([-u1[1], u1[0], 0]);
    const north = V.unit(V.cross(u1, east));
    const az = mod360(Math.atan2(V.dot(v, east), V.dot(v, north)) / DEG);

    return { jdL, jdA, u1, u2, n, theta, inc, raan, az,
             moonDist: mn.dist, moonDec: Math.asin(u2[2]) / DEG, phi_gc: sp.phi_gc };
}

function coastSolutions(geo, params) {
    const { rp, tAscent, thAscent, maxCoastSec, transferSec } = params;
    const omPark = Math.sqrt(MU / (rp*rp*rp)) / DEG;   // deg/s
    const d = geo.moonDist;
    const sols = [];
    const f = (tc) => {
        const dnu = mod360(geo.theta - thAscent - omPark * tc);
        if (dnu < 95 || dnu > 250) return null;
        const e = (d - rp) / (rp - d * Math.cos(dnu * DEG));
        if (e <= 0 || e >= 1) return null;
        return { val: kepTOF(rp, e, dnu) - (transferSec - tAscent - tc), dnu, e };
    };
    let prev = null, prevTc = 0;
    for (let tc = 0; tc <= maxCoastSec; tc += 5) {
        const cur = f(tc);
        if (cur && prev && prev.val * cur.val < 0) {
            const tcr = prevTc + 5 * prev.val / (prev.val - cur.val);
            const r = f(tcr);
            if (r) {
                const a = rp / (1 - r.e);
                const vPer = Math.sqrt(MU * (2/rp - 1/a));
                sols.push({
                    coast: tcr, dnu: r.dnu, e: r.e, a,
                    apoAlt: a * (1 + r.e) - R_MEAN,
                    dv: vPer - Math.sqrt(MU / rp),
                    c3: vPer*vPer - 2*MU/rp,
                    tofSec: transferSec - params.tAscent - tcr
                });
            }
        }
        prev = cur; prevTc = tc;
    }
    return sols;
}

function readParams() {
    const rpAlt = parseFloat(document.getElementById('parking-alt').value) || 200;
    return {
        rp: R_MEAN + rpAlt,
        tAscent: parseFloat(document.getElementById('ascent-time').value) || 800,
        thAscent: parseFloat(document.getElementById('ascent-arc').value) || 20,
        maxCoastSec: (parseFloat(document.getElementById('max-coast').value) || 20) * 60,
        azMin: parseFloat(document.getElementById('az-min').value),
        azMax: parseFloat(document.getElementById('az-max').value),
        periTarget: parseFloat(document.getElementById('peri-target').value) || 150,
        lunarGrav: document.getElementById('lunar-grav').checked,
        biasSide: parseInt(document.getElementById('bias-side').value) || 0,
        polarLock: document.getElementById('polar-lock').checked
    };
}

// 扫描一天: 返回每个格点 {hUTC, az(实飞=二体+偏置), az2b, feasible, sol, geo}
function scanDay(jd0, site, transferDays, params, azBias = 0, stepMin = 2) {
    const pts = [];
    for (let m = 0; m <= 1440; m += stepMin) {
        const jdL = jd0 + m / 1440;
        const geo = solveGeometry(jdL, site, transferDays);
        const azReal = mod360(geo.az + azBias);
        const inCorr = azReal >= params.azMin && azReal <= params.azMax;
        let sol = null;
        if (inCorr) {
            const sols = coastSolutions(geo, { ...params, transferSec: transferDays * 86400 });
            sol = sols.length ? sols[0] : null;
        }
        pts.push({ hUTC: m / 60, jdL, az: azReal, az2b: geo.az, geo, sol,
                   feasible: inCorr && !!sol });
    }
    // 提取连续窗口
    const windows = [];
    let cur = null;
    for (const p of pts) {
        if (p.feasible) {
            if (!cur) cur = { open: p.hUTC, close: p.hUTC, pts: [] };
            cur.close = p.hUTC; cur.pts.push(p);
        } else if (cur) { windows.push(cur); cur = null; }
    }
    if (cur) windows.push(cur);
    return { pts, windows };
}

/* =====================================================================
 *  第二部分补充: 月球引力修正 (四体: 地 + 月 + 日 点质量, 地心系含间接项)
 *
 *  二体闭合把 LTO 瞄准 "月球@到达" 的质心方向, 忽略了月球引力对末段的
 *  拉弯加速 — 实飞会提前数小时到达且近月点大幅偏低 (常见直接撞月)。
 *  修正: 在同一发射几何上加偏置, 用 RK4 四体积分打靶。
 *  基本约束 (2×2): 近月点高度 = 目标 且 近月点时刻 = 设定到达时刻,
 *    控制 (射向偏置 dAz, 滑行偏置 dtc)。该问题有 ± 两支解 (B 面偏置在
 *    月球两侧), 捕获轨道面取向不受控 (是自由落出的量)。
 *  极轨捕获 (3×3): 追加约束 近月点轨道面法向 ⊥ 月自转轴 (i_M = 90°),
 *    控制量增加 TLI 速度偏置 dvp — 与全数值参考解 (DE440s/DOP853 +
 *    LM 微分修正) 同构; 到达贴近远地点时 TOF 对能量极端敏感, 亚 m/s 的
 *    dvp 即可与 dtc 联动转动 B 矢量方向数十度, 故该约束可控。
 *  窗口扫描采用 "窗口中心标定一次, 全天套用偏置" (偏置随发射时刻近似
 *  常数, 表现为窗口整体平移, 宽度不变)。
 * ===================================================================== */
const MU_SUN = 1.32712440018e11;   // km^3/s^2
const R_MOON = 1737.4;             // km
const AU_KM = 1.495978707e8;
// IAU 月自转轴 (J2000 赤道系; 与 of-date 差 <0.3°)
const MOON_POLE = (() => {
    const ra = 269.9949 * DEG, de = 66.5392 * DEG;
    return [Math.cos(de) * Math.cos(ra), Math.cos(de) * Math.sin(ra),
            Math.sin(de)];
})();

function sunEqVec(jdUTC) {         // 太阳位置矢量 (km, 赤道系 of date)
    const n = jdTT(jdUTC) - 2451545.0;
    const L = mod360(280.460 + 0.9856474 * n);
    const g = mod360(357.528 + 0.9856003 * n) * DEG;
    const lam = mod360(L + 1.915 * Math.sin(g) + 0.020 * Math.sin(2 * g)) * DEG;
    const r = (1.00014 - 0.01671 * Math.cos(g) - 0.00014 * Math.cos(2 * g)) * AU_KM;
    return V.scale(eclToEq([Math.cos(lam), Math.sin(lam), 0], jdUTC), r);
}

// 历表: 均匀采样月/日位置 + Catmull-Rom 三次插值 (积分内上万次调用)
function buildEph(jd0, jd1) {
    const step = 0.02;
    const n = Math.ceil((jd1 - jd0) / step) + 4;
    const moonP = [], sunP = [];
    for (let k = 0; k < n; k++) {
        const jd = jd0 - step + k * step;
        const m = moonEq(jd);
        moonP.push(V.scale(m.u, m.dist));
        sunP.push(sunEqVec(jd));
    }
    const interp = (arr, jd) => {
        const x = (jd - (jd0 - step)) / step;
        const i = Math.max(1, Math.min(arr.length - 3, Math.floor(x)));
        const t = x - i, out = [0, 0, 0];
        for (let c = 0; c < 3; c++) {
            const p0 = arr[i - 1][c], p1 = arr[i][c],
                  p2 = arr[i + 1][c], p3 = arr[i + 2][c];
            out[c] = p1 + 0.5 * t * (p2 - p0 + t * (2*p0 - 5*p1 + 4*p2 - p3
                     + t * (3 * (p1 - p2) + p3 - p0)));
        }
        return out;
    };
    const moonF = jd => interp(moonP, jd);
    const hv = 0.004;
    return { moon: moonF, sun: jd => interp(sunP, jd),
             moonVel: jd => {
                 const a = moonF(jd - hv), b = moonF(jd + hv);
                 return [(b[0] - a[0]) / (2 * hv * 86400),
                         (b[1] - a[1]) / (2 * hv * 86400),
                         (b[2] - a[2]) / (2 * hv * 86400)];
             } };
}

function accel4b(jdRef, tSec, r, eph) {
    const jd = jdRef + tSec / 86400;
    const rm = eph.moon(jd), rs = eph.sun(jd);
    const dm = V.sub(r, rm), ds = V.sub(r, rs);
    const dR = V.norm(r), nm = V.norm(dm), nrm = V.norm(rm),
          ns = V.norm(ds), nrs = V.norm(rs);
    const cE = -MU / (dR*dR*dR), cM1 = -MU_MOON / (nm*nm*nm),
          cM2 = -MU_MOON / (nrm*nrm*nrm), cS1 = -MU_SUN / (ns*ns*ns),
          cS2 = -MU_SUN / (nrs*nrs*nrs);
    return [cE*r[0] + cM1*dm[0] + cM2*rm[0] + cS1*ds[0] + cS2*rs[0],
            cE*r[1] + cM1*dm[1] + cM2*rm[1] + cS1*ds[1] + cS2*rs[1],
            cE*r[2] + cM1*dm[2] + cM2*rm[2] + cS1*ds[2] + cS2*rs[2]];
}

function rk4step(jdRef, t, y, h, eph) {
    const f = (tt, yy) => {
        const a = accel4b(jdRef, tt, [yy[0], yy[1], yy[2]], eph);
        return [yy[3], yy[4], yy[5], a[0], a[1], a[2]];
    };
    const k1 = f(t, y);
    const k2 = f(t + h/2, y.map((v, i) => v + h/2 * k1[i]));
    const k3 = f(t + h/2, y.map((v, i) => v + h/2 * k2[i]));
    const k4 = f(t + h, y.map((v, i) => v + h * k3[i]));
    return y.map((v, i) => v + h/6 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]));
}

function parabolaMinT(p0, p1, p2) {   // 三点 (t,d) 抛物线极小值位置
    const a = ((p0[1] - p1[1]) / (p0[0] - p1[0])
             - (p1[1] - p2[1]) / (p1[0] - p2[0])) / (p0[0] - p2[0]);
    if (Math.abs(a) < 1e-15) return p1[0];
    const b = (p0[1] - p1[1]) / (p0[0] - p1[0]) - a * (p0[0] + p1[0]);
    return -b / (2 * a);
}

// 从 TLI 状态积分, 返回近月点 { tPeri, dPeri, y (近月点附近状态) }
function propagateToPerilune(jdRef, y0, tEnd, eph) {
    let t = 0, y = y0.slice();
    const rec = [];
    let best = [0, Infinity, y0.slice()];
    const dMoon = (tt, yy) => {
        const rm = eph.moon(jdRef + tt / 86400);
        return Math.hypot(yy[0]-rm[0], yy[1]-rm[1], yy[2]-rm[2]);
    };
    while (t < tEnd) {
        const d = dMoon(t, y);
        rec.push([t, d, y]);
        if (rec.length > 3) rec.shift();
        if (d < best[1]) best = [t, d, y];
        if (rec.length === 3 && rec[1][1] < rec[0][1] && rec[1][1] <= rec[2][1]
            && rec[1][1] < 60000) {
            const tv = parabolaMinT(rec[0], rec[1], rec[2]);
            const a = ((rec[0][1] - rec[1][1]) / (rec[0][0] - rec[1][0])
                     - (rec[1][1] - rec[2][1]) / (rec[1][0] - rec[2][0]))
                     / (rec[0][0] - rec[2][0]);
            const dmin = rec[1][1] - a * (rec[1][0] - tv) * (rec[1][0] - tv);
            return { tPeri: tv, dPeri: dmin, y: rec[1][2] };
        }
        const rE = Math.hypot(y[0], y[1], y[2]);
        let h = rE < 20000 ? 10 : rE < 60000 ? 40 : 240;
        if (d < 60000) h = Math.min(h, 30);
        if (d < 25000) h = 15;
        if (d < 8000) h = 5;
        if (d < 2500) h = 2;
        h = Math.min(h, tEnd - t);
        y = rk4step(jdRef, t, y, h, eph);
        t += h;
    }
    return { tPeri: best[0], dPeri: best[1], y: best[2] };
}

// 与 solveGeometry/coastSolutions 相同几何下的 TLI 状态 (射向加偏置 dAz)
function tliStateBiased(geo, params, tcSec, vp, dAzDeg) {
    const u1 = geo.u1;
    const east = V.unit([-u1[1], u1[0], 0]);
    const north = V.cross(u1, east);
    const azr = (geo.az + dAzDeg) * DEG;
    const vdir = V.add(V.scale(east, Math.sin(azr)), V.scale(north, Math.cos(azr)));
    const n = V.unit(V.cross(u1, vdir));
    const omPark = Math.sqrt(MU / (params.rp * params.rp * params.rp)) / DEG;
    const ang = (params.thAscent + omPark * tcSec) * DEG;
    const ehat = V.add(V.scale(u1, Math.cos(ang)),
                       V.scale(V.cross(n, u1), Math.sin(ang)));
    const phat = V.cross(n, ehat);
    return [...V.scale(ehat, params.rp), ...V.scale(phat, vp)];
}

function solve3(J, b) {              // 3×3 Cramer
    const det = J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])
              - J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])
              + J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
    if (!isFinite(det) || Math.abs(det) < 1e-15) return null;
    const rep = (M, c, v) => M.map((row, i) => row.map((x, j) => j === c ? v[i] : x));
    const d3 = M => M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
               - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
               + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);
    return [d3(rep(J, 0, b))/det, d3(rep(J, 1, b))/det, d3(rep(J, 2, b))/det];
}

// 打靶主函数: 见文件头注释。params.biasSide ∈ {0 自动, ±1}; params.polarLock
function lunarCorrection(geo, sol, params) {
    const vp0 = Math.sqrt(MU * (2 / params.rp - 1 / sol.a));
    const eph = buildEph(geo.jdL - 0.05, geo.jdA + 1.6);
    const rTarget = R_MOON + (params.periTarget || 200);
    const tcLo = -sol.coast, tcHi = params.maxCoastSec - sol.coast;
    const F = (dAz, dtc, dvp = 0) => {
        const jdTLI = geo.jdL + (params.tAscent + sol.coast + dtc) / 86400;
        const tofSec = (geo.jdA - jdTLI) * 86400;
        const y0 = tliStateBiased(geo, params, sol.coast + dtc, vp0 + dvp, dAz);
        const p = propagateToPerilune(jdTLI, y0, tofSec + 1.5 * 86400, eph);
        const jdP = jdTLI + p.tPeri / 86400;
        const rm = eph.moon(jdP), vm = eph.moonVel(jdP);
        const rr = [p.y[0]-rm[0], p.y[1]-rm[1], p.y[2]-rm[2]];
        const vr = [p.y[3]-vm[0], p.y[4]-vm[1], p.y[5]-vm[2]];
        const pol = V.dot(V.unit(V.cross(rr, vr)), MOON_POLE);
        return { f0: p.dPeri - rTarget, f1: (p.tPeri - tofSec) / 3600,
                 pol, iM: Math.acos(Math.max(-1, Math.min(1, pol))) / DEG, p };
    };
    const nrm = c => Math.hypot(c.f0 / 80, c.f1);   // 80 km ≈ 1 h 等权
    const raw = F(0, 0);

    // 时刻求解: 固定 dAz, 对 dtc 括号+二分 (f1 随 dtc 单调增, ~0.1 h/s)
    const solveTime = (azFix, d0, c0) => {
        if (Math.abs(c0.f1) < 0.04) return [d0, c0];
        const dir = c0.f1 < 0 ? 1 : -1;          // 提前到达 → 加长滑行
        let lo = d0, glo = c0.f1, hi = null, cl = c0;
        for (const st of [30, 60, 120, 240, 480]) {
            const d = Math.max(tcLo, Math.min(tcHi, d0 + dir * st));
            const c = F(azFix, d);
            if ((c.f1 < 0) !== (c0.f1 < 0)) { hi = d; cl = c; break; }
            lo = d; glo = c.f1; cl = c;
            if (d === tcLo || d === tcHi) break;
        }
        if (hi === null) return [lo, cl];
        let bestD = hi, best = cl;
        for (let k = 0; k < 9 && Math.abs(best.f1) > 0.04; k++) {
            const dm = (lo + hi) / 2;
            const cm = F(azFix, dm);
            if ((cm.f1 < 0) === (glo < 0)) { lo = dm; glo = cm.f1; }
            else hi = dm;
            bestD = dm; best = cm;
        }
        return [bestD, best];
    };

    // 高度求解: 固定 dtc, 横向扫 dAz — "近月距-侧偏" 是 V 型褶皱, 不能盲用
    // 多维牛顿; 从中心向外找变号区间, 弦截细化; 返回两侧候选 (含 i_M)
    const solveAlt = (dtcFix, sidePref) => {
        const cBase = F(0, dtcFix);
        const refine = (aLo, cLo, aHi, cHi) => {
            let lo = [aLo, cLo.f0], hi = [aHi, cHi.f0];
            let cBest = cHi; cBest.dAz = aHi;
            for (let k = 0; k < 7 && Math.abs(cBest.f0) > 20; k++) {
                const am = lo[0] + (hi[0] - lo[0]) * lo[1] / (lo[1] - hi[1]);
                const cm = F(am, dtcFix); cm.dAz = am;
                cBest = cm;
                if ((cm.f0 < 0) === (lo[1] < 0)) lo = [am, cm.f0];
                else hi = [am, cm.f0];
            }
            return cBest;
        };
        const cands = [];
        for (const sgn of (sidePref ? [sidePref] : [1, -1])) {
            let prevA = 0, prevC = cBase;
            for (const step of [0.35, 0.7, 1.1, 1.6, 2.3, 3.2]) {
                const a = sgn * step;
                const c = F(a, dtcFix);
                if ((prevC.f0 < 0) !== (c.f0 < 0)) {
                    const cb = refine(prevA, prevC, a, c);
                    cb.side = sgn;
                    cands.push(cb);
                    break;
                }
                prevA = a; prevC = c;
            }
        }
        if (!cands.length) { cBase.dAz = 0; cBase.side = 0; return [cBase]; }
        cands.sort((p, q) => Math.abs(p.dAz) - Math.abs(q.dAz));  // 自动 = 最小修正
        return cands;
    };

    // 两遍交替 (时刻↔高度弱耦合, 两遍即入抛光半径)
    let [dtc, cT] = solveTime(0, 0, raw);
    let cands = solveAlt(dtc, params.biasSide || 0);
    let cA = cands[0];
    let dAz = cA.dAz;
    [dtc, cT] = solveTime(dAz, dtc, cA);
    cands = solveAlt(dtc, 0);                    // 双侧都解 (供展示/换支)
    const want = params.biasSide || cA.side || 0;
    cA = cands.find(c => c.side === want) || cands[0];
    dAz = cA.dAz;
    let dvp = 0, cur = cA, ok = false;
    const branches = cands.map(c => ({ side: c.side, dAz: c.dAz, iM: c.iM }));

    if (params.polarLock) {
        // 3×3 抛光: (dAz, dtc, dvp) ↔ (高度, 时刻, 极轨面 n̂·p̂=0)
        // 能量杠杆极陡 (贴远地点) → 小差分步 + 强钳制
        let hA = 0.05, hT = 10, hV = 1.5e-4;
        const nrm3 = c => Math.hypot(c.f0 / 80, c.f1,
            Math.asin(Math.max(-1, Math.min(1, c.pol))) / DEG / 1.5);
        cur = F(dAz, dtc, dvp);
        for (let it = 0; it < 12; it++) {
            if (Math.abs(cur.f0) < 3 && Math.abs(cur.f1) < 0.02
                && Math.abs(cur.pol) < 8.7e-4) { ok = true; break; }
            const fA = F(dAz + hA, dtc, dvp), fT = F(dAz, dtc + hT, dvp),
                  fV = F(dAz, dtc, dvp + hV);
            const J = [
                [(fA.f0-cur.f0)/hA, (fT.f0-cur.f0)/hT, (fV.f0-cur.f0)/hV],
                [(fA.f1-cur.f1)/hA, (fT.f1-cur.f1)/hT, (fV.f1-cur.f1)/hV],
                [(fA.pol-cur.pol)/hA, (fT.pol-cur.pol)/hT, (fV.pol-cur.pol)/hV]];
            const dx = solve3(J, [-cur.f0, -cur.f1, -cur.pol]);
            if (!dx) break;
            dx[0] = Math.max(-0.5, Math.min(0.5, dx[0]));
            dx[1] = Math.max(-120, Math.min(120, dx[1]));
            dx[2] = Math.max(-5e-4, Math.min(5e-4, dx[2]));
            let lam = 1, nxt = null;
            for (let k = 0; k < 4; k++) {
                nxt = F(dAz + lam*dx[0], dtc + lam*dx[1], dvp + lam*dx[2]);
                if (nrm3(nxt) < nrm3(cur)) break;
                lam /= 2;
            }
            dAz += lam*dx[0]; dtc += lam*dx[1]; dvp += lam*dx[2];
            cur = nxt;
            hA = Math.max(0.012, hA * 0.6);
            hT = Math.max(3, hT * 0.6);
            hV = Math.max(5e-5, hV * 0.6);
        }
        if (!ok) ok = Math.abs(cur.f0) < 25 && Math.abs(cur.f1) < 0.15
                      && Math.abs(cur.pol) < 0.02;
    } else {
        // 2×2 抛光 (捕获面不约束)
        let hA = 0.05, hT = 12;
        let x = [dAz, dtc];
        for (let it = 0; it < 5; it++) {
            if (Math.abs(cur.f0) < 3 && Math.abs(cur.f1) < 0.02) { ok = true; break; }
            const fA = F(x[0] + hA, x[1]), fT = F(x[0], x[1] + hT);
            const J00 = (fA.f0 - cur.f0) / hA, J01 = (fT.f0 - cur.f0) / hT,
                  J10 = (fA.f1 - cur.f1) / hA, J11 = (fT.f1 - cur.f1) / hT;
            const det = J00 * J11 - J01 * J10;
            if (!isFinite(det) || Math.abs(det) < 1e-12) break;
            let dxA = -( J11 * cur.f0 - J01 * cur.f1) / det;
            let dxT = -(-J10 * cur.f0 + J00 * cur.f1) / det;
            dxA = Math.max(-0.4, Math.min(0.4, dxA));
            dxT = Math.max(-90, Math.min(90, dxT));
            let lam = 1, nxt = null;
            for (let k = 0; k < 3; k++) {
                nxt = F(x[0] + lam * dxA, x[1] + lam * dxT);
                if (nrm(nxt) < nrm(cur)) break;
                lam /= 2;
            }
            x = [x[0] + lam * dxA, x[1] + lam * dxT];
            cur = nxt;
            hA = Math.max(0.015, hA * 0.5);
            hT = Math.max(4, hT * 0.5);
        }
        if (!ok) ok = Math.abs(cur.f0) < 25 && Math.abs(cur.f1) < 0.15;
        dAz = x[0]; dtc = x[1];
    }
    return { dAz, dtc, dvp, side: cA.side || 0, iM: cur.iM, branches,
             polar: !!params.polarLock,
             periAlt: cur.p.dPeri - R_MOON, periDtH: cur.f1,
             rawAlt: raw.p.dPeri - R_MOON, rawDtH: raw.f1,
             converged: ok };
}

/* ===== 原 index.html 内联脚本块 2/4（未改动） ===== */
/* =====================================================================
 *  第三部分: Three.js 场景 (r144, 借用 satobserver 的地球昼夜着色器)
 *  场景轴: sceneX = 赤道系X, sceneY = 赤道系Z(北极), sceneZ = -赤道系Y
 * ===================================================================== */
const SCALE = 26 / 384400;      // 月距 ~26 单位
const EARTH_VIS_R = 1.0;
const MOON_VIS_R = 0.30;

function eqToScene(v) { return new THREE.Vector3(v[0], v[2], -v[1]); }

const viewport = document.getElementById('viewport');
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(50, 1, 0.05, 20000);
camera.position.set(24, 22, 34);

const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
renderer.setClearColor(0x05070c);
renderer.outputEncoding = THREE.sRGBEncoding;
viewport.appendChild(renderer.domElement);

const controls = new THREE.OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
controls.dampingFactor = 0.08;
controls.minDistance = 1.6;
controls.maxDistance = 300;

scene.add(new THREE.AmbientLight(0x223344, 0.7));
const sunLight = new THREE.DirectionalLight(0xffffff, 1.15);
sunLight.position.set(50, 0, 0);
scene.add(sunLight);

// 星空
{
    const g = new THREE.BufferGeometry(), pos = [];
    for (let i = 0; i < 2600; i++) {
        const th = Math.random() * 2 * Math.PI, ph = Math.acos(2 * Math.random() - 1);
        const r = 5000 + Math.random() * 2500;
        pos.push(r*Math.sin(ph)*Math.cos(th), r*Math.cos(ph), r*Math.sin(ph)*Math.sin(th));
    }
    g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
    scene.add(new THREE.Points(g, new THREE.PointsMaterial({
        color: 0xbbccdd, size: 1.6, sizeAttenuation: false })));
}

/* ---- 地球: 昼夜贴图着色器 (satobserver globe3d.js 方案) ---- */
const earthGroup = new THREE.Group();   // rotation.y = GMST → ECEF 内容转到惯性系
scene.add(earthGroup);

let dayNightMat = null, phongMat = null, earthMesh = null, texturesOK = false;
{
    const loader = new THREE.TextureLoader();
    const geo = new THREE.SphereGeometry(EARTH_VIS_R, 96, 64);
    phongMat = new THREE.MeshPhongMaterial({
        color: 0x2255aa, emissive: 0x101c30, specular: 0x222831, shininess: 12 });
    earthMesh = new THREE.Mesh(geo, phongMat);
    earthGroup.add(earthMesh);

    let loaded = 0;
    const dayTex = loader.load('assets/earth_day.jpg', onTex, undefined, onErr);
    const nightTex = loader.load('assets/earth_night.jpg', onTex, undefined, onErr);
    function onTex() {
        if (++loaded < 2) return;
        [dayTex, nightTex].forEach(t => {
            t.encoding = THREE.sRGBEncoding;
            t.anisotropy = renderer.capabilities.getMaxAnisotropy();
        });
        dayNightMat = new THREE.ShaderMaterial({
            uniforms: {
                dayMap: { value: dayTex }, nightMap: { value: nightTex },
                sunDir: { value: new THREE.Vector3(1, 0, 0) }   // earthGroup 局部系
            },
            vertexShader: `
                varying vec2 vUv; varying vec3 vNormalL;
                void main() {
                    vUv = uv;
                    vNormalL = normal;
                    gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
                }`,
            fragmentShader: `
                uniform sampler2D dayMap, nightMap; uniform vec3 sunDir;
                varying vec2 vUv; varying vec3 vNormalL;
                void main() {
                    float d = dot(normalize(vNormalL), normalize(sunDir));
                    float k = smoothstep(-0.05, 0.05, d);
                    vec3 day = texture2D(dayMap, vUv).rgb * (0.22 + 0.95 * clamp(d, 0.0, 1.0));
                    vec3 night = texture2D(nightMap, vUv).rgb * vec3(1.05, 0.95, 0.78);
                    gl_FragColor = vec4(pow(mix(night, day, k), vec3(1.0/2.2)), 1.0);
                }`
        });
        texturesOK = true;
        applyEarthMaterial();
        // 纹理晚于首次计算到位时, 重建场景以刷新晨昏线朝向
        if (typeof updateScene3D === 'function' && typeof selected !== 'undefined' && selected)
            updateScene3D();
    }
    function onErr() { /* 纹理加载失败: 保持纯色 (file:// 打开时属正常) */ }
}
function applyEarthMaterial() {
    const want = document.getElementById('show-daynight').checked && texturesOK;
    earthMesh.material = want && dayNightMat ? dayNightMat : phongMat;
}

// 大气辉光 (satobserver)
{
    const mat = new THREE.ShaderMaterial({
        vertexShader: `
            varying vec3 vN;
            void main() {
                vN = normalize(normalMatrix * normal);
                gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
            }`,
        fragmentShader: `
            varying vec3 vN;
            void main() {
                float i = pow(clamp(0.66 - dot(vN, vec3(0.0, 0.0, 1.0)), 0.0, 1.0), 5.0);
                gl_FragColor = vec4(0.30, 0.55, 1.0, 1.0) * i * 0.65;
            }`,
        side: THREE.BackSide, blending: THREE.AdditiveBlending,
        transparent: true, depthWrite: false
    });
    earthGroup.add(new THREE.Mesh(new THREE.SphereGeometry(EARTH_VIS_R * 1.045, 64, 48), mat));
}

// ECEF → earthGroup 局部坐标 (地固系, group 自转带动)
function geoToLocal(latDeg, lonDeg, r) {
    const la = latDeg * DEG, lo = lonDeg * DEG;
    return new THREE.Vector3(r*Math.cos(la)*Math.cos(lo), r*Math.sin(la), -r*Math.cos(la)*Math.sin(lo));
}

/* ---- 图标 marker: 恒定屏幕尺寸的小标识 (canvas glyph → sprite) ----
   depthTest 开启: 转到地球背面时被正确遮挡, 不穿透 */
const SCALED_GLYPHS = [];
function makeGlyphSprite(drawFn, sizeWorld, centerY = 0.5) {
    const c = document.createElement('canvas');
    c.width = 128; c.height = 128;
    drawFn(c.getContext('2d'));
    const mat = new THREE.SpriteMaterial({
        map: new THREE.CanvasTexture(c), transparent: true, depthTest: true });
    const s = new THREE.Sprite(mat);
    s.center.set(0.5, centerY);
    s.userData.glyphSize = sizeWorld;
    s.scale.set(sizeWorld, sizeWorld, 1);
    SCALED_GLYPHS.push(s);
    return s;
}
function drawPinGlyph(x) {   // 发射场: 定位销
    x.lineWidth = 7; x.strokeStyle = 'rgba(6,26,14,0.9)'; x.fillStyle = '#22cc66';
    x.beginPath();
    x.arc(64, 46, 30, Math.PI * 0.82, Math.PI * 0.18);
    x.lineTo(64, 118);
    x.closePath(); x.fill(); x.stroke();
    x.fillStyle = '#0d1117';
    x.beginPath(); x.arc(64, 46, 12, 0, 6.3); x.fill();
}
function drawBurstGlyph(x) {  // TLI 点火: 星芒
    x.fillStyle = '#00ff88'; x.strokeStyle = 'rgba(0,40,20,0.9)'; x.lineWidth = 5;
    x.beginPath();
    for (let i = 0; i < 16; i++) {
        const r = i % 2 === 0 ? 56 : 22;
        const a = i / 16 * 2 * Math.PI - Math.PI / 2;
        x[i ? 'lineTo' : 'moveTo'](64 + r * Math.cos(a), 64 + r * Math.sin(a));
    }
    x.closePath(); x.fill(); x.stroke();
    x.fillStyle = '#eafff4';
    x.beginPath(); x.arc(64, 64, 10, 0, 6.3); x.fill();
}
function drawShipGlyph(x) {   // 飞船: 小火箭
    x.strokeStyle = 'rgba(10,10,20,0.9)'; x.lineWidth = 6; x.fillStyle = '#f2f5f8';
    x.beginPath();
    x.moveTo(64, 10);
    x.quadraticCurveTo(88, 42, 84, 78);
    x.lineTo(44, 78);
    x.quadraticCurveTo(40, 42, 64, 10);
    x.closePath(); x.fill(); x.stroke();
    x.fillStyle = '#8b949e';
    x.beginPath(); x.moveTo(46, 58); x.lineTo(26, 86); x.lineTo(46, 80); x.closePath(); x.fill();
    x.beginPath(); x.moveTo(82, 58); x.lineTo(102, 86); x.lineTo(82, 80); x.closePath(); x.fill();
    x.fillStyle = '#d0782e';
    x.beginPath(); x.moveTo(52, 82); x.lineTo(64, 116); x.lineTo(76, 82); x.closePath(); x.fill();
    x.fillStyle = '#4d8fd0';
    x.beginPath(); x.arc(64, 50, 9, 0, 6.3); x.fill();
}

// 发射场标记 (针尖锚定地表) + 射向走廊扇面 (地固系子对象)
const siteMarker = makeGlyphSprite(drawPinGlyph, 0.5, 0.04);
earthGroup.add(siteMarker);
const corridorGroup = new THREE.Group();
earthGroup.add(corridorGroup);

function greatCircleFrom(latDeg, lonDeg, azDeg, arcDeg, n = 40) {
    // 从 (lat,lon) 沿方位角 az 的大圆弧, 返回 ECEF 局部坐标点列
    const la = latDeg * DEG, lo = lonDeg * DEG, az = azDeg * DEG;
    const pts = [];
    for (let i = 0; i <= n; i++) {
        const d = (i / n) * arcDeg * DEG;
        const lat2 = Math.asin(Math.sin(la)*Math.cos(d) + Math.cos(la)*Math.sin(d)*Math.cos(az));
        const lon2 = lo + Math.atan2(Math.sin(az)*Math.sin(d)*Math.cos(la),
                                     Math.cos(d) - Math.sin(la)*Math.sin(lat2));
        pts.push(geoToLocal(lat2/DEG, lon2/DEG, EARTH_VIS_R * 1.004));
    }
    return pts;
}
function rebuildCorridor(site, azMin, azMax) {
    while (corridorGroup.children.length) {
        const c = corridorGroup.children[0];
        c.geometry && c.geometry.dispose();
        corridorGroup.remove(c);
    }
    const mat = new THREE.LineBasicMaterial({ color: 0xd0782e, transparent: true, opacity: 0.55 });
    const matEdge = new THREE.LineBasicMaterial({ color: 0xd0782e, transparent: true, opacity: 0.9 });
    const nAz = 8;
    for (let i = 0; i <= nAz; i++) {
        const az = azMin + (i / nAz) * (azMax - azMin);
        const pts = greatCircleFrom(site.lat, site.lon, az, 22);
        const g = new THREE.BufferGeometry().setFromPoints(pts);
        corridorGroup.add(new THREE.Line(g, (i === 0 || i === nAz) ? matEdge : mat));
    }
    siteMarker.position.copy(geoToLocal(site.lat, site.lon, EARTH_VIS_R * 1.012));
    siteLabel.position.copy(geoToLocal(site.lat, site.lon, EARTH_VIS_R * 1.16));
}

/* ---- 标签 ---- */
function makeTextSprite(text, color = '#ffffff', scale = 1.0) {
    const canvas = document.createElement('canvas');
    canvas.width = 512; canvas.height = 128;
    const ctx = canvas.getContext('2d');
    ctx.font = '64px "PingFang SC", Arial';
    ctx.fillStyle = color;
    ctx.textAlign = 'center';
    ctx.strokeStyle = 'rgba(0,0,0,0.7)'; ctx.lineWidth = 6;
    ctx.strokeText(text, 256, 84); ctx.fillText(text, 256, 84);
    const mat = new THREE.SpriteMaterial({
        map: new THREE.CanvasTexture(canvas), transparent: true, depthTest: true });
    const s = new THREE.Sprite(mat);
    s.scale.set(2.6 * scale, 0.65 * scale, 1);
    return s;
}
const earthLabel = makeTextSprite('地球', '#6699cc', 0.6); earthLabel.position.set(0, 1.5, 0); scene.add(earthLabel);
const moonLabel = makeTextSprite('月球@发射', '#aaaaaa', 0.6); scene.add(moonLabel);
const moonArrLabel = makeTextSprite('月球@到达', '#f0883e', 0.6); scene.add(moonArrLabel);
const siteLabel = makeTextSprite('发射场', '#22cc66', 0.5); earthGroup.add(siteLabel);
const tliLabel = makeTextSprite('TLI', '#00ff88', 0.45); scene.add(tliLabel);
const sunLabelSpr = makeTextSprite('☀ 太阳方向', '#ffcc44', 0.55); scene.add(sunLabelSpr);
// 标签随相机距离缩放, 保持近似恒定屏幕大小 (satobserver rescaleScreenItems 思路)
const SCALED_LABELS = [[earthLabel, 0.6], [moonLabel, 0.6], [moonArrLabel, 0.6],
                       [siteLabel, 0.5], [tliLabel, 0.45], [sunLabelSpr, 0.55]];
const _lblWP = new THREE.Vector3();

/* ---- 示意位置矢量: 地心 → 发射场 (r_S), 地心 → 月心@到达 (r_M) ---- */
// r_S: earthGroup 子对象 (地固系), 随地球自转; 穿透地球绘制, 避免箭杆被地球遮挡
// r_M: 场景惯性系, 指向 solveGeometry 中到达时刻月心 u2 (与"月球@到达"一致); 同样穿透地球
const siteVecArrow = new THREE.ArrowHelper(
    new THREE.Vector3(1, 0, 0), new THREE.Vector3(0, 0, 0),
    EARTH_VIS_R * 1.18, 0x22cc66, 0.18, 0.10);
siteVecArrow.line.material.depthTest = false;
siteVecArrow.cone.material.depthTest = false;
siteVecArrow.line.renderOrder = 8;
siteVecArrow.cone.renderOrder = 8;
earthGroup.add(siteVecArrow);
const siteVecLabel = makeTextSprite('r_S 发射场位置', '#22cc66', 0.5);
siteVecLabel.material.depthTest = false;
siteVecLabel.renderOrder = 9;
earthGroup.add(siteVecLabel);
const moonVecArrow = new THREE.ArrowHelper(
    new THREE.Vector3(1, 0, 0), new THREE.Vector3(0, 0, 0),
    MOON_VIS_R * 2, 0xccaa55, 0.7, 0.32);
moonVecArrow.line.material.depthTest = false;
moonVecArrow.cone.material.depthTest = false;
moonVecArrow.line.renderOrder = 8;
moonVecArrow.cone.renderOrder = 8;
scene.add(moonVecArrow);
const moonVecLabel = makeTextSprite('r_M 月心@到达', '#ccaa55', 0.5);
scene.add(moonVecLabel);
SCALED_LABELS.push([siteVecLabel, 0.5], [moonVecLabel, 0.5]);
siteVecArrow.visible = siteVecLabel.visible = false;
moonVecArrow.visible = moonVecLabel.visible = false;
function rescaleLabels() {
    for (const [spr, base] of SCALED_LABELS) {
        if (!spr.visible) continue;
        spr.getWorldPosition(_lblWP);
        const k = base * Math.max(0.12, Math.min(2.2, camera.position.distanceTo(_lblWP) / 42));
        spr.scale.set(2.6 * k, 0.65 * k, 1);
    }
    // 图标 marker: 几乎完全恒定的屏幕尺寸
    for (const spr of SCALED_GLYPHS) {
        if (!spr.visible) continue;
        spr.getWorldPosition(_lblWP);
        const s = spr.userData.glyphSize *
            Math.max(0.02, Math.min(8, camera.position.distanceTo(_lblWP) / 42));
        spr.scale.set(s, s, 1);
    }
}

/* ---- 月球 ---- */
function makeMoonTexture() {
    const c = document.createElement('canvas'); c.width = 512; c.height = 256;
    const x = c.getContext('2d');
    x.fillStyle = '#9a9a9a'; x.fillRect(0, 0, 512, 256);
    let seed = 9973;
    const rnd = () => (seed = (seed * 16807) % 2147483647) / 2147483647;
    for (const [mx, my, rx, ry] of [[180,100,60,40],[280,120,45,35],[350,90,30,25],[120,140,50,30],[400,80,35,28]]) {
        x.beginPath(); x.ellipse(mx, my, rx, ry, rnd()*0.5, 0, 6.3);
        x.fillStyle = `rgba(70,70,75,${0.3+rnd()*0.2})`; x.fill();
    }
    for (let i = 0; i < 260; i++) {
        const b = Math.floor(rnd()*50 + 110);
        x.beginPath(); x.arc(rnd()*512, rnd()*256, rnd()*7 + 1, 0, 6.3);
        x.fillStyle = `rgb(${b},${b},${b})`; x.fill();
    }
    return new THREE.CanvasTexture(c);
}
const moonMesh = new THREE.Mesh(
    new THREE.SphereGeometry(MOON_VIS_R, 32, 32),
    new THREE.MeshPhongMaterial({ map: makeMoonTexture(), shininess: 4 }));
scene.add(moonMesh);
const moonGhost = new THREE.Mesh(
    new THREE.SphereGeometry(MOON_VIS_R * 0.85, 16, 16),
    new THREE.MeshBasicMaterial({ color: 0xf0883e, transparent: true, opacity: 0.28, wireframe: true }));
scene.add(moonGhost);

/* ---- 轨道线组 ---- */
const moonOrbitLine = new THREE.Line(new THREE.BufferGeometry(),
    new THREE.LineBasicMaterial({ color: 0x8899aa, transparent: true, opacity: 0.45 }));
scene.add(moonOrbitLine);

const lunarPlaneMesh = new THREE.Mesh(
    new THREE.RingGeometry(1.5, 30, 72),
    new THREE.MeshBasicMaterial({ color: 0xccaa44, side: THREE.DoubleSide, transparent: true, opacity: 0.05 }));
scene.add(lunarPlaneMesh);

const ltoPlaneMesh = new THREE.Mesh(
    new THREE.RingGeometry(1.2, 30, 72),
    new THREE.MeshBasicMaterial({ color: 0xd0782e, side: THREE.DoubleSide, transparent: true, opacity: 0.05 }));
ltoPlaneMesh.visible = false;
scene.add(ltoPlaneMesh);

const equatorPlaneMesh = new THREE.Mesh(
    new THREE.RingGeometry(1.2, 28, 72),
    new THREE.MeshBasicMaterial({ color: 0x44cc66, side: THREE.DoubleSide, transparent: true, opacity: 0.04 }));
equatorPlaneMesh.rotation.x = Math.PI / 2;   // 场景基面即赤道面
equatorPlaneMesh.visible = false;
scene.add(equatorPlaneMesh);

const eclipticPlaneMesh = new THREE.Mesh(
    new THREE.RingGeometry(1.2, 29, 72),
    new THREE.MeshBasicMaterial({ color: 0x4466cc, side: THREE.DoubleSide, transparent: true, opacity: 0.05 }));
eclipticPlaneMesh.visible = false;
scene.add(eclipticPlaneMesh);

const trajGroup = new THREE.Group();   // 上升段+停泊+LTO
scene.add(trajGroup);
const ascentGroup = new THREE.Group();   // 上升段
const coastGroup = new THREE.Group();    // 停泊轨道 + 滑行弧
const ltoGroup = new THREE.Group();      // LTO + 剩余轨道虚线
trajGroup.add(ascentGroup, coastGroup, ltoGroup);

// 航天器 与 TLI 点火 图标
const scMarker = makeGlyphSprite(drawShipGlyph, 0.55);
scMarker.visible = false;
scene.add(scMarker);
const tliMarker = makeGlyphSprite(drawBurstGlyph, 0.42);
tliMarker.visible = false;
scene.add(tliMarker);
const scTrail = new THREE.Line(new THREE.BufferGeometry(),
    new THREE.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.5 }));
scene.add(scTrail);

function setPlaneQuat(mesh, nEq) {
    const nS = eqToScene(nEq).normalize();
    mesh.quaternion.setFromUnitVectors(new THREE.Vector3(0, 0, 1), nS);
}
function disposeGroup(g) {
    while (g.children.length) {
        const c = g.children[0];
        c.geometry && c.geometry.dispose();
        g.remove(c);
    }
}

/* ===== 原 index.html 内联脚本块 3/4（未改动） ===== */
/* =====================================================================
 *  第四部分: 图表 (SVG, 深色主题; 调色板已验证 CVD 安全)
 *  蓝 #4d8fd0 方位角 | 橙 #d0782e 窗口/走廊 | 青 #2ba8a0 滑行
 * ===================================================================== */
const CH = {
    ink: '#e0e0e0', muted: '#8b949e', grid: '#1e2d3d',
    blue: '#4d8fd0', orange: '#d0782e', teal: '#2ba8a0'
};
function svgEl(tag, attrs) {
    const e = document.createElementNS('http://www.w3.org/2000/svg', tag);
    for (const k in attrs) e.setAttribute(k, attrs[k]);
    return e;
}
function clearSvg(svg) { while (svg.firstChild) svg.removeChild(svg.firstChild); }

// 方位角-时刻图: 曲线 + 走廊带 + 可行段
function drawAzChart(svg, pts, params, selHUTC, tz) {
    clearSvg(svg);
    const W = svg.clientWidth || 320, H = 150;
    svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
    const mL = 30, mR = 8, mT = 8, mB = 18;
    const x = h => mL + (h / 24) * (W - mL - mR);
    const azLo = 40, azHi = 180;
    const y = az => mT + (1 - (az - azLo) / (azHi - azLo)) * (H - mT - mB);

    // 走廊带
    svg.appendChild(svgEl('rect', {
        x: mL, y: y(params.azMax), width: W - mL - mR,
        height: Math.abs(y(params.azMin) - y(params.azMax)),
        fill: CH.orange, opacity: 0.13 }));
    // 网格 + 轴
    for (let az = 60; az <= 180; az += 30) {
        svg.appendChild(svgEl('line', { x1: mL, y1: y(az), x2: W - mR, y2: y(az), stroke: CH.grid, 'stroke-width': 1 }));
        const t = svgEl('text', { x: mL - 4, y: y(az) + 3, fill: CH.muted, 'font-size': 9, 'text-anchor': 'end' });
        t.textContent = az; svg.appendChild(t);
    }
    for (let h = 0; h <= 24; h += 6) {
        const t = svgEl('text', { x: x(h), y: H - 5, fill: CH.muted, 'font-size': 9, 'text-anchor': 'middle' });
        t.textContent = `${mod360(h + tz) % 24}时`; svg.appendChild(t);
    }
    // 方位角曲线 (蓝), 可行段覆盖橙色粗线
    let d = '', dF = '', pen = false, penF = false;
    for (const p of pts) {
        const az = p.az;
        if (az < azLo || az > azHi) { pen = false; penF = false; continue; }
        d += (pen ? 'L' : 'M') + x(p.hUTC).toFixed(1) + ' ' + y(az).toFixed(1); pen = true;
        if (p.feasible) { dF += (penF ? 'L' : 'M') + x(p.hUTC).toFixed(1) + ' ' + y(az).toFixed(1); penF = true; }
        else penF = false;
    }
    svg.appendChild(svgEl('path', { d, fill: 'none', stroke: CH.blue, 'stroke-width': 1.6 }));
    if (dF) svg.appendChild(svgEl('path', { d: dF, fill: 'none', stroke: CH.orange, 'stroke-width': 3, 'stroke-linecap': 'round' }));
    // 选中时刻
    if (selHUTC != null) {
        svg.appendChild(svgEl('line', { x1: x(selHUTC), y1: mT, x2: x(selHUTC), y2: H - mB,
            stroke: CH.orange, 'stroke-width': 1, 'stroke-dasharray': '3 3' }));
    }
    // 悬停
    addHover(svg, mL, W - mR, hFrac => {
        const h = hFrac * 24;
        const p = pts[Math.round(h * 60 / 2)];
        if (!p) return '';
        return `${fmtHM(h, tz)} 射向 ${p.az.toFixed(1)}°${p.feasible ? ' ✓可发射' : ''}`;
    }, 'chart-az-hint');
}

function drawCoastChart(svg, pts, tz) {
    clearSvg(svg);
    const W = svg.clientWidth || 320, H = 110;
    svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
    const mL = 30, mR = 8, mT = 8, mB = 18;
    const feas = pts.filter(p => p.feasible && p.sol);
    if (!feas.length) {
        const t = svgEl('text', { x: W/2, y: H/2, fill: CH.muted, 'font-size': 11, 'text-anchor': 'middle' });
        t.textContent = '当日无可行窗口'; svg.appendChild(t); return;
    }
    const cMax = Math.max(...feas.map(p => p.sol.coast / 60)) * 1.15 + 2;
    const hMin = Math.min(...feas.map(p => p.hUTC)) - 0.5;
    const hMax = Math.max(...feas.map(p => p.hUTC)) + 0.5;
    const x = h => mL + (h - hMin) / (hMax - hMin) * (W - mL - mR);
    const y = c => mT + (1 - c / cMax) * (H - mT - mB);
    for (let c = 0; c <= cMax; c += Math.ceil(cMax / 4 / 5) * 5 || 5) {
        svg.appendChild(svgEl('line', { x1: mL, y1: y(c), x2: W - mR, y2: y(c), stroke: CH.grid }));
        const t = svgEl('text', { x: mL - 4, y: y(c) + 3, fill: CH.muted, 'font-size': 9, 'text-anchor': 'end' });
        t.textContent = c; svg.appendChild(t);
    }
    let d = '', pen = false, lastH = null;
    for (const p of feas) {
        if (lastH != null && p.hUTC - lastH > 0.1) pen = false;
        d += (pen ? 'L' : 'M') + x(p.hUTC).toFixed(1) + ' ' + y(p.sol.coast / 60).toFixed(1);
        pen = true; lastH = p.hUTC;
    }
    svg.appendChild(svgEl('path', { d, fill: 'none', stroke: CH.teal, 'stroke-width': 2 }));
    for (let i = 0; i <= 2; i++) {
        const h = hMin + 0.5 + (i / 2) * (hMax - hMin - 1);
        const t = svgEl('text', { x: x(h), y: H - 5, fill: CH.muted, 'font-size': 9, 'text-anchor': 'middle' });
        t.textContent = fmtHM(h, tz); svg.appendChild(t);
    }
    const lbl = svgEl('text', { x: mL + 4, y: mT + 9, fill: CH.teal, 'font-size': 9 });
    lbl.textContent = '滑行 min'; svg.appendChild(lbl);
}

function drawMonthChart(svg, days, selIdx, onPick) {
    clearSvg(svg);
    const W = svg.clientWidth || 320, H = 110;
    svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
    const mL = 26, mR = 6, mT = 8, mB = 20;
    const wMax = Math.max(10, ...days.map(d => d.minutes));
    const bw = (W - mL - mR) / days.length;
    const y = v => mT + (1 - v / (wMax * 1.1)) * (H - mT - mB);
    for (let v = 0; v <= wMax; v += wMax > 120 ? 60 : 30) {
        svg.appendChild(svgEl('line', { x1: mL, y1: y(v), x2: W - mR, y2: y(v), stroke: CH.grid }));
        const t = svgEl('text', { x: mL - 3, y: y(v) + 3, fill: CH.muted, 'font-size': 8.5, 'text-anchor': 'end' });
        t.textContent = v; svg.appendChild(t);
    }
    days.forEach((d, i) => {
        const h = Math.max(0, y(0) - y(d.minutes));
        const r = svgEl('rect', {
            x: mL + i * bw + 1, y: d.minutes > 0 ? y(d.minutes) : y(0) - 1,
            width: Math.max(1, bw - 2), height: d.minutes > 0 ? h : 1,
            rx: 2, fill: i === selIdx ? CH.orange : CH.blue,
            opacity: d.minutes > 0 ? 0.9 : 0.25, cursor: 'pointer'
        });
        r.addEventListener('click', () => onPick(i));
        const tip = svgEl('title', {});
        tip.textContent = `${d.label}: ${d.minutes.toFixed(0)} min`;
        r.appendChild(tip);
        svg.appendChild(r);
        if (i % 5 === 0) {
            const t = svgEl('text', { x: mL + i * bw + bw/2, y: H - 6, fill: CH.muted, 'font-size': 8.5, 'text-anchor': 'middle' });
            t.textContent = d.label.slice(5); svg.appendChild(t);
        }
    });
}

function addHover(svg, x0, x1, fmt, hintId) {
    const hint = document.getElementById(hintId);
    const base = hint ? hint.textContent : '';
    svg.addEventListener('mousemove', ev => {
        const r = svg.getBoundingClientRect();
        const fx = (ev.clientX - r.left) * (svg.viewBox.baseVal.width / r.width);
        if (fx < x0 || fx > x1) return;
        const s = fmt((fx - x0) / (x1 - x0));
        if (hint && s) hint.textContent = s;
    });
    svg.addEventListener('mouseleave', () => { if (hint) hint.textContent = base; });
}

function fmtHM(hUTC, tz) {
    let h = hUTC + tz, d = '';
    if (h >= 24) { h -= 24; d = '+1'; } else if (h < 0) { h += 24; d = '-1'; }
    const hh = Math.floor(h), mm = Math.floor((h - hh) * 60);
    return `${String(hh).padStart(2,'0')}:${String(mm).padStart(2,'0')}${d}`;
}
function fmtJD(jd, tz, withDate = true) {
    const ms = (jd - 2440587.5) * 86400000 + tz * 3600000;
    const dt = new Date(ms);
    const s = dt.toISOString();
    return withDate ? s.slice(0, 16).replace('T', ' ') : s.slice(11, 16);
}

/* ---- 月相 ---- */
function drawMoonPhase(canvas, elong) {
    const ctx = canvas.getContext('2d');
    const w = canvas.width, h = canvas.height, cx = w/2, cy = h/2, r = w/2 - 2;
    ctx.clearRect(0, 0, w, h);
    ctx.beginPath(); ctx.arc(cx, cy, r, 0, 6.3);
    ctx.fillStyle = '#c8c4bc'; ctx.fill();
    const cp = Math.cos(elong * DEG);
    ctx.save(); ctx.beginPath();
    if (elong <= 180) {
        ctx.arc(cx, cy, r, Math.PI/2, -Math.PI/2, false);
        ctx.ellipse(cx, cy, r*Math.abs(cp), r, 0, -Math.PI/2, Math.PI/2, elong > 90);
    } else {
        ctx.arc(cx, cy, r, -Math.PI/2, Math.PI/2, false);
        ctx.ellipse(cx, cy, r*Math.abs(cp), r, 0, Math.PI/2, -Math.PI/2, elong <= 270);
    }
    ctx.closePath(); ctx.fillStyle = '#1a1a2e'; ctx.fill(); ctx.restore();
}
function phaseName(e) {
    const names = ['新月 (朔)','蛾眉月','上弦月','盈凸月','满月 (望)','亏凸月','下弦月','残月'];
    return names[Math.floor(mod360(e + 22.5) / 45)] || '新月 (朔)';
}

/* ===== 原 index.html 内联脚本块 4/4（未改动） ===== */
/* =====================================================================
 *  第五部分: 主流程 — 计算、场景更新、动画、UI
 * ===================================================================== */
let dayResult = null;        // scanDay 输出
let selected = null;         // { geo, sol, params, site } 选定发射时刻的完整解
let monthDays = null;
let visBasis = null;         // 3D 可视化基向量缓存 (动画用)

function currentSite() { return LAUNCH_SITES[document.getElementById('launch-site').value]; }

function calcDay() {
    const site = currentSite();
    const dateStr = document.getElementById('launch-date').value;
    if (!dateStr) return;
    const jd0 = julianDateUTC(Date.parse(dateStr + 'T00:00:00Z'));
    const transferDays = parseFloat(document.getElementById('transfer-time').value);
    const params = readParams();

    // 二体扫描 → (可选) 窗口中心四体标定射向偏置 → 带偏置重扫
    let scan = scanDay(jd0, site, transferDays, params);
    let corr = null, azBias = 0;
    if (params.lunarGrav) {
        let jdRefL;
        if (scan.windows.length) {
            const w = scan.windows.reduce((a, b) =>
                (b.close - b.open) > (a.close - a.open) ? b : a);
            jdRefL = jd0 + (w.open + w.close) / 2 / 24;
        } else {
            const mid = (params.azMin + params.azMax) / 2;
            const p = scan.pts.reduce((a, b) =>
                Math.abs(b.az2b - mid) < Math.abs(a.az2b - mid) ? b : a);
            jdRefL = p.jdL;
        }
        const geoR = solveGeometry(jdRefL, site, transferDays);
        const solsR = coastSolutions(geoR, { ...params, transferSec: transferDays * 86400 });
        if (solsR.length) {
            corr = lunarCorrection(geoR, solsR[0], params);
            if (corr.converged) {
                azBias = corr.dAz;
                scan = scanDay(jd0, site, transferDays, params, azBias);
            }
        }
    }
    dayResult = { ...scan, jd0, site, transferDays, params, corr, azBias };
    document.getElementById('chart-az-hint').textContent = (corr && corr.converged)
        ? `曲线: 实飞射向 = 二体共面射向 ${azBias >= 0 ? '+' : ''}${azBias.toFixed(2)}° (月引力偏置, 窗口中心标定)`
        : '曲线: 瞄准"月球@到达"共面所需射向';
    renderWindowList();
    if (dayResult.windows.length) selectInWindow(dayResult.windows[0], 0.5);
    else { selected = null; updateResultPanel(); updateScene3D(); }
    drawAzChart(document.getElementById('chart-az'), dayResult.pts, params,
        selected ? (selected.geo.jdL - jd0) * 24 : null, site.tz);
    drawCoastChart(document.getElementById('chart-coast'), dayResult.pts, site.tz);
}

function renderWindowList() {
    const box = document.getElementById('window-list');
    box.innerHTML = '';
    if (!dayResult || !dayResult.windows.length) {
        box.innerHTML = '<div class="no-window">当日射向走廊内无共面发射机会<br>(试试相邻日期或放宽走廊/滑行限制)</div>';
        document.getElementById('in-window-box').style.display = 'none';
        return;
    }
    const { site } = dayResult;
    dayResult.windows.forEach((w, i) => {
        const div = document.createElement('div');
        div.className = 'window-item' + (i === 0 ? ' active' : '');
        const width = (w.close - w.open) * 60;
        div.innerHTML = `<div class="w-time">${fmtHM(w.open, site.tz)} – ${fmtHM(w.close, site.tz)}</div>
            <div class="w-sub">宽度 ~${width.toFixed(0)} min · 射向 ${w.pts[0].az.toFixed(1)}°→${w.pts[w.pts.length-1].az.toFixed(1)}°</div>`;
        div.addEventListener('click', () => {
            document.querySelectorAll('.window-item').forEach(el => el.classList.remove('active'));
            div.classList.add('active');
            selectInWindow(w, 0.5);
        });
        box.appendChild(div);
    });
    document.getElementById('in-window-box').style.display = '';
}

function selectInWindow(w, frac) {
    const { site, transferDays, params, jd0 } = dayResult;
    const hUTC = w.open + frac * (w.close - w.open);
    const jdL = jd0 + hUTC / 24;
    const geo = solveGeometry(jdL, site, transferDays);
    const sols = coastSolutions(geo, { ...params, transferSec: transferDays * 86400 });
    if (!sols.length) return;
    selected = { geo, sol: sols[0], params, site, transferDays, window: w,
                 azBias: dayResult.azBias, corr: dayResult.corr };
    document.getElementById('in-window-value').textContent =
        `${fmtHM(hUTC, site.tz)} ${site.tzName} (射向 ${mod360(geo.az + dayResult.azBias).toFixed(1)}°)`;
    const slider = document.getElementById('in-window');
    slider.value = Math.round(frac * 100);
    updateResultPanel();
    updateScene3D();
    drawAzChart(document.getElementById('chart-az'), dayResult.pts, params, hUTC, site.tz);
    resetAnimation();
}

function scanMonth() {
    const site = currentSite();
    const dateStr = document.getElementById('launch-date').value;
    if (!dateStr) return;
    const jdBase = julianDateUTC(Date.parse(dateStr + 'T00:00:00Z'));
    const transferDays = parseFloat(document.getElementById('transfer-time').value);
    const params = readParams();
    // 月引力偏置: 沿用最近一次当日标定 (逐日重标定太慢; 偏置月内漂移 <0.5°)
    const azBias = (params.lunarGrav && dayResult && dayResult.corr
                    && dayResult.corr.converged) ? dayResult.azBias : 0;
    monthDays = [];
    for (let d = 0; d < 30; d++) {
        const jd0 = jdBase + d;
        // 粗扫: 5 min 步长, 只测走廊 (滑行解一般存在, 由 maxCoast 约束)
        let minutes = 0;
        for (let m = 0; m < 1440; m += 5) {
            const geo = solveGeometry(jd0 + m / 1440, site, transferDays);
            const azReal = mod360(geo.az + azBias);
            if (azReal >= params.azMin && azReal <= params.azMax) {
                const sols = coastSolutions(geo, { ...params, transferSec: transferDays * 86400 });
                if (sols.length) minutes += 5;
            }
        }
        const dt = new Date((jd0 - 2440587.5) * 86400000);
        monthDays.push({ jd0, minutes, label: dt.toISOString().slice(0, 10) });
    }
    document.getElementById('month-box').style.display = '';
    const pickDay = i => {
        document.getElementById('launch-date').value = monthDays[i].label;
        calcDay();
        drawMonthChart(document.getElementById('chart-month'), monthDays, i, pickDay);
    };
    drawMonthChart(document.getElementById('chart-month'), monthDays, 0, pickDay);
}

/* ---- 结果面板 ---- */
function updateResultPanel() {
    const q = id => document.getElementById(id);
    if (!selected) {
        ['r-launch','r-az','r-coast','r-tli','r-dnu','r-inc','r-raan','r-ecc',
         'r-apo','r-dv','r-c3','r-gap','r-arrival','r-mdec','r-mdist','r-minc',
         'r-bias','r-dvp','r-im','r-peri','r-perit','r-raw']
            .forEach(id => q(id).textContent = '--');
        return;
    }
    const { geo, sol, site, params } = selected;
    const corr = selected.corr, azBias = selected.azBias || 0;
    const corrOn = corr && corr.converged;
    const nMoon = moonPlaneNormal(geo.jdA);
    const gap = Math.acos(Math.max(-1, Math.min(1, V.dot(geo.n, nMoon)))) / DEG;
    const iMoon = Math.acos(nMoon[2]) / DEG;
    q('r-launch').textContent = fmtJD(geo.jdL, site.tz) + ' ' + site.tzName;
    q('r-az').textContent = mod360(geo.az + azBias).toFixed(2) + '°';
    const coastEff = sol.coast + (corrOn ? corr.dtc : 0);
    q('r-coast').textContent = (coastEff / 60).toFixed(1) + ' min';
    const jdTLI = geo.jdL + (params.tAscent + coastEff) / 86400;
    q('r-tli').textContent = fmtJD(jdTLI, site.tz) + ' ' + site.tzName;
    q('r-dnu').textContent = sol.dnu.toFixed(1) + '°';
    q('r-inc').textContent = geo.inc.toFixed(2) + '°';
    q('r-raan').textContent = geo.raan.toFixed(2) + '°';
    q('r-ecc').textContent = sol.e.toFixed(5);
    q('r-apo').textContent = Math.round(sol.apoAlt).toLocaleString() + ' km';
    q('r-dv').textContent = (sol.dv + (corrOn ? corr.dvp : 0)).toFixed(3) + ' km/s';
    q('r-c3').textContent = sol.c3.toFixed(2) + ' km²/s²';
    q('r-gap').textContent = gap.toFixed(2) + '°';
    q('r-arrival').textContent = fmtJD(geo.jdA, site.tz) + ' ' + site.tzName;
    q('r-mdec').textContent = geo.moonDec.toFixed(2) + '°';
    q('r-mdist').textContent = Math.round(geo.moonDist).toLocaleString() + ' km';
    q('r-minc').textContent = iMoon.toFixed(2) + '°';

    // 月球引力修正面板
    if (corrOn) {
        const sg = v => (v >= 0 ? '+' : '') + v;
        const other = (corr.branches || []).find(b => b.side !== corr.side);
        q('r-bias').textContent = `${sg(corr.dAz.toFixed(2))}° [${corr.side > 0 ? '+' : '−'}支]`
            + (other ? ` (另支 ${sg(other.dAz.toFixed(2))}°)` : '');
        q('r-dvp').textContent = `${sg(corr.dtc.toFixed(0))} s / `
            + `${sg((corr.dvp * 1000).toFixed(2))} m/s`;
        q('r-im').textContent = corr.polar
            ? `${corr.iM.toFixed(1)}° (约束 90°)`
            : `${corr.iM.toFixed(1)}° (未约束)`
              + (other ? `; 另支 ${other.iM.toFixed(1)}°` : '');
        q('r-peri').textContent = `${Math.round(corr.periAlt)} km (标定点)`;
        q('r-perit').textContent = `${(corr.periDtH * 60).toFixed(0)} min`;
        const rawA = corr.rawAlt < 0 ? '撞月'
            : `近月点 ${Math.round(corr.rawAlt).toLocaleString()} km`;
        const rawT = corr.rawDtH <= 0 ? `提前 ${(-corr.rawDtH).toFixed(1)} h`
            : `滞后 ${corr.rawDtH.toFixed(1)} h`;
        q('r-raw').textContent = `${rawA} · ${rawT}`;
    } else if (params.lunarGrav) {
        ['r-bias','r-dvp','r-im','r-peri','r-perit','r-raw']
            .forEach(id => q(id).textContent = '未收敛');
    } else {
        ['r-bias','r-dvp','r-im','r-peri','r-perit','r-raw']
            .forEach(id => q(id).textContent = '关');
    }

    // 月相 @到达
    const mA = moonEcliptic(geo.jdA);
    const elong = mod360(mA.lonDeg - sunEclipticLon(geo.jdA));
    drawMoonPhase(document.getElementById('moon-phase-canvas'), elong);
    document.getElementById('phase-name').textContent = phaseName(elong);
    document.getElementById('phase-detail').textContent =
        `到达时光照 ${((1 - Math.cos(elong*DEG)) / 2 * 100).toFixed(0)}% · 角距 ${elong.toFixed(0)}°`;
}

/* ---- 3D 场景更新 ---- */
function updateScene3D() {
    disposeGroup(ascentGroup);
    disposeGroup(coastGroup);
    disposeGroup(ltoGroup);
    scMarker.visible = false;
    tliMarker.visible = false;
    scTrail.geometry.dispose();
    scTrail.geometry = new THREE.BufferGeometry();
    visBasis = null;
    siteVecArrow.visible = siteVecLabel.visible = false;
    moonVecArrow.visible = moonVecLabel.visible = false;
    if (!selected) return;

    const { geo, sol, params, site } = selected;

    // 地球自转到发射时刻
    earthGroup.rotation.y = gmstDeg(geo.jdL) * DEG;
    rebuildCorridor(site, params.azMin, params.azMax);

    // 太阳方向 (赤道系)
    const sl = sunEclipticLon(geo.jdL) * DEG;
    const sunEq = eclToEq([Math.cos(sl), Math.sin(sl), 0], geo.jdL);
    const sunS = eqToScene(sunEq).normalize();
    sunLight.position.copy(sunS.clone().multiplyScalar(60));
    sunLabelSpr.position.copy(sunS.clone().multiplyScalar(31));
    if (dayNightMat) {
        // earthGroup 局部系 = 场景系绕 Y 转 -gmst
        const inv = sunS.clone().applyAxisAngle(new THREE.Vector3(0, 1, 0), -earthGroup.rotation.y);
        dayNightMat.uniforms.sunDir.value.copy(inv);
    }

    // 月球 @发射 与 @到达, 月球轨道 (±14 天真实历表)
    const mL = moonEq(geo.jdL);
    moonMesh.position.copy(eqToScene(V.scale(mL.u, mL.dist * SCALE)));
    moonLabel.position.copy(moonMesh.position).add(new THREE.Vector3(0, 0.75, 0));
    moonGhost.position.copy(eqToScene(V.scale(geo.u2, geo.moonDist * SCALE)));
    moonArrLabel.position.copy(moonGhost.position).add(new THREE.Vector3(0, 0.7, 0));

    // 位置矢量箭头: r_S 地心→发射场 (随地球自转); r_M 地心→到达时刻月心 u2
    siteVecArrow.setDirection(geoToLocal(site.lat, site.lon, 1).normalize());
    siteVecArrow.setLength(EARTH_VIS_R * 1.18, 0.18, 0.10);
    siteVecLabel.position.copy(geoToLocal(site.lat, site.lon, EARTH_VIS_R * 1.18 * 0.75));
    moonVecArrow.setDirection(moonGhost.position.clone().normalize());
    moonVecArrow.setLength(Math.max(0.2, moonGhost.position.length() - MOON_VIS_R * 0.9), 0.7, 0.32);
    moonVecLabel.position.copy(moonGhost.position).multiplyScalar(0.55)
        .add(new THREE.Vector3(0, 0.55, 0));
    {
        const pts = [];
        for (let i = 0; i <= 120; i++) {
            const jd = geo.jdL - 13.66 + (i / 120) * 27.32;
            const m = moonEq(jd);
            pts.push(eqToScene(V.scale(m.u, m.dist * SCALE)));
        }
        moonOrbitLine.geometry.dispose();
        moonOrbitLine.geometry = new THREE.BufferGeometry().setFromPoints(pts);
    }
    setPlaneQuat(lunarPlaneMesh, moonPlaneNormal(geo.jdL));
    setPlaneQuat(ltoPlaneMesh, geo.n);
    {
        const e = obliquityDeg(jdTT(geo.jdL)) * DEG;
        setPlaneQuat(eclipticPlaneMesh, [0, -Math.sin(e), Math.cos(e)]);  // 黄道法向 = Rx(ε)·(0,0,1)
    }

    /* 轨迹几何 (惯性系):
       u1=发射场方向, v=n×u1 运动方向;
       入轨点 @ θa, TLI @ Θ = θa + ω·coast, LTO: ν 从 0 起, 位置角 = Θ + ν */
    const u1 = eqToScene(geo.u1), n3 = eqToScene(geo.n).normalize();
    const v3 = new THREE.Vector3().crossVectors(n3, u1).normalize();
    const rpVis = EARTH_VIS_R * (params.rp / R_MEAN);
    const raVis = geo.moonDist * SCALE;
    const omPark = Math.sqrt(MU / Math.pow(params.rp, 3)) / DEG;  // deg/s
    const thA = params.thAscent;
    const TH = thA + omPark * sol.coast;      // TLI 位置角
    const dirAt = ang => u1.clone().multiplyScalar(Math.cos(ang * DEG))
                          .add(v3.clone().multiplyScalar(Math.sin(ang * DEG)));

    // 视觉椭圆: r(0)=rpVis, r(Δν)=raVis 两点定形
    const eV = (raVis - rpVis) / (rpVis - raVis * Math.cos(sol.dnu * DEG));
    const pV = rpVis * (1 + eV);
    const ltoR = nu => pV / (1 + eV * Math.cos(nu * DEG));

    visBasis = { u1, v3, n3, rpVis, omPark, thA, TH, eV, pV, geo, sol, params, dirAt, ltoR };

    // 上升段 (绿)
    {
        const pts = [];
        for (let i = 0; i <= 30; i++) {
            const f = i / 30;
            const r = EARTH_VIS_R + (rpVis - EARTH_VIS_R) * Math.pow(f, 0.7);
            pts.push(dirAt(f * thA).multiplyScalar(r));
        }
        ascentGroup.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts),
            new THREE.LineBasicMaterial({ color: 0x22cc66, linewidth: 2 })));
    }
    // 停泊轨道整圈 (暗) + 滑行弧 (青, 亮)
    {
        const ring = [];
        for (let i = 0; i <= 96; i++) ring.push(dirAt(i / 96 * 360).multiplyScalar(rpVis));
        coastGroup.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(ring),
            new THREE.LineBasicMaterial({ color: 0x2ba8a0, transparent: true, opacity: 0.25 })));
        const arc = [];
        const span = TH - thA;
        for (let i = 0; i <= 48; i++) arc.push(dirAt(thA + i / 48 * span).multiplyScalar(rpVis));
        coastGroup.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(arc),
            new THREE.LineBasicMaterial({ color: 0x2ba8a0, linewidth: 2 })));
    }
    // TLI 标记
    {
        const p = dirAt(TH).multiplyScalar(rpVis);
        tliMarker.visible = true;
        tliMarker.position.copy(p);
        tliLabel.position.copy(p.clone().multiplyScalar(1.35));
    }
    // LTO (橙) ν: 0 → Δν
    {
        const pts = [];
        for (let i = 0; i <= 140; i++) {
            const nu = i / 140 * sol.dnu;
            pts.push(dirAt(TH + nu).multiplyScalar(ltoR(nu)));
        }
        ltoGroup.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts),
            new THREE.LineBasicMaterial({ color: 0xf0883e, linewidth: 2 })));
        // 剩余轨道虚线
        const rest = [];
        for (let i = 1; i <= 40; i++) {
            const nu = sol.dnu + i / 40 * Math.min(40, 355 - sol.dnu);
            const r = ltoR(nu);
            if (r > 0 && r < 60) rest.push(dirAt(TH + nu).multiplyScalar(r));
        }
        if (rest.length > 1) {
            const dl = new THREE.Line(new THREE.BufferGeometry().setFromPoints(rest),
                new THREE.LineDashedMaterial({ color: 0xf0883e, dashSize: 0.4, gapSize: 0.3,
                                               transparent: true, opacity: 0.3 }));
            dl.computeLineDistances();
            ltoGroup.add(dl);
        }
    }
    applyToggles();
}

/* ---- 时间轴动画 ---- */
const anim = { playing: false, tSec: 0, lastMs: 0, trail: [] };

function resetAnimation() {
    anim.playing = false; anim.tSec = 0; anim.trail = [];
    document.getElementById('tl-play').textContent = '▶';
    document.getElementById('tl-slider').value = 0;
    applyAnimation();
}
function totalMissionSec() {
    return selected ? selected.transferDays * 86400 : 1;
}
function applyAnimation() {
    if (!selected || !visBasis) return;
    const { geo, sol, params } = selected;
    const { dirAt, ltoR, rpVis, omPark, thA, TH } = visBasis;
    const t = anim.tSec;
    const jdNow = geo.jdL + t / 86400;

    // 地球自转 + 昼夜
    earthGroup.rotation.y = gmstDeg(jdNow) * DEG;
    if (dayNightMat) {
        const sl2 = sunEclipticLon(jdNow) * DEG;
        const sEq = eclToEq([Math.cos(sl2), Math.sin(sl2), 0], jdNow);
        const sS = eqToScene(sEq).normalize();
        sunLight.position.copy(sS.clone().multiplyScalar(60));
        const inv = sS.clone().applyAxisAngle(new THREE.Vector3(0, 1, 0), -earthGroup.rotation.y);
        dayNightMat.uniforms.sunDir.value.copy(inv);
    }
    // 月球实时位置
    const m = moonEq(jdNow);
    moonMesh.position.copy(eqToScene(V.scale(m.u, m.dist * SCALE)));
    moonLabel.position.copy(moonMesh.position).add(new THREE.Vector3(0, 0.75, 0));

    // 航天器
    let pos = null, phase = '—';
    const tTLI = params.tAscent + sol.coast;
    if (t <= 0) { pos = null; }
    else if (t < params.tAscent) {
        const f = t / params.tAscent;
        const r = EARTH_VIS_R + (rpVis - EARTH_VIS_R) * Math.pow(f, 0.7);
        pos = dirAt(f * thA).multiplyScalar(r);
        phase = '上升段';
    } else if (t < tTLI) {
        const ang = thA + omPark * (t - params.tAscent);
        pos = dirAt(ang).multiplyScalar(rpVis);
        phase = '停泊滑行';
    } else {
        const nu = kepNuAt(params.rp, sol.e, t - tTLI);
        const nuC = Math.min(nu, sol.dnu + 5);
        pos = dirAt(TH + nuC).multiplyScalar(ltoR(nuC));
        phase = 'LTO 巡航';
    }
    if (pos) {
        scMarker.visible = true;
        scMarker.position.copy(pos);
        anim.trail.push(pos.clone());
        if (anim.trail.length > 400) anim.trail.shift();
        scTrail.geometry.dispose();
        scTrail.geometry = new THREE.BufferGeometry().setFromPoints(anim.trail);
    } else {
        scMarker.visible = false;
    }
    document.getElementById('tl-phase').textContent = phase;
    document.getElementById('tl-time').textContent =
        fmtJD(jdNow, selected.site.tz) + ' ' + selected.site.tzName;
    document.getElementById('tl-slider').value = Math.round(t / totalMissionSec() * 1000);
}

/* ---- 渲染循环 ---- */
function animate(nowMs) {
    requestAnimationFrame(animate);
    if (anim.playing && selected) {
        const dt = Math.min(0.1, (nowMs - anim.lastMs) / 1000) || 0;
        anim.tSec += dt * parseFloat(document.getElementById('tl-speed').value);
        if (anim.tSec >= totalMissionSec()) { anim.tSec = totalMissionSec(); anim.playing = false;
            document.getElementById('tl-play').textContent = '▶'; }
        applyAnimation();
    }
    anim.lastMs = nowMs;
    controls.update();
    rescaleLabels();
    renderer.render(scene, camera);
}

/* ---- UI 事件 ---- */
function applyToggles() {
    const g = id => document.getElementById(id).checked;
    ascentGroup.visible = g('show-ascent');
    coastGroup.visible = g('show-coast');
    ltoGroup.visible = g('show-lto');
    lunarPlaneMesh.visible = g('show-lunar-plane');
    moonOrbitLine.visible = g('show-lunar-plane');
    ltoPlaneMesh.visible = g('show-lto-plane');
    equatorPlaneMesh.visible = g('show-equator');
    eclipticPlaneMesh.visible = g('show-ecliptic');
    corridorGroup.visible = g('show-corridor');
    [earthLabel, moonLabel, moonArrLabel, siteLabel, tliLabel, sunLabelSpr]
        .forEach(s => s.visible = g('show-labels'));
    const hasSel = !!selected;
    siteVecArrow.visible = hasSel && g('show-site-vec');
    moonVecArrow.visible = hasSel && g('show-moon-vec');
    siteVecLabel.visible = siteVecArrow.visible && g('show-labels');
    moonVecLabel.visible = moonVecArrow.visible && g('show-labels');
    applyEarthMaterial();
}
['show-ascent','show-coast','show-lto',
 'show-lunar-plane','show-lto-plane','show-equator','show-ecliptic',
 'show-corridor','show-labels','show-daynight','show-site-vec','show-moon-vec'].forEach(id =>
    document.getElementById(id).addEventListener('change', applyToggles));

document.getElementById('btn-calc').addEventListener('click', calcDay);
document.getElementById('btn-month').addEventListener('click', scanMonth);

document.getElementById('launch-site').addEventListener('change', () => {
    const s = currentSite();
    document.getElementById('site-lat').textContent =
        `${Math.abs(s.lat).toFixed(3)}°${s.lat >= 0 ? 'N' : 'S'}`;
    document.getElementById('site-lon').textContent =
        `${Math.abs(s.lon).toFixed(3)}°${s.lon >= 0 ? 'E' : 'W'}`;
    document.getElementById('az-min').value = s.azMin;
    document.getElementById('az-max').value = s.azMax;
    document.getElementById('site-corridor').textContent = `${s.azMin}°–${s.azMax}°`;
    calcDay();
});
['az-min','az-max'].forEach(id => document.getElementById(id).addEventListener('change', () => {
    document.getElementById('site-corridor').textContent =
        `${document.getElementById('az-min').value}°–${document.getElementById('az-max').value}°`;
    calcDay();
}));
document.getElementById('transfer-time').addEventListener('input', e => {
    const d = parseFloat(e.target.value);
    document.getElementById('transfer-time-value').textContent =
        `${d.toFixed(2)} 天 (${(d * 24).toFixed(1)} h)`;
});
document.getElementById('transfer-time').addEventListener('change', calcDay);
['launch-date','parking-alt','max-coast','ascent-time','ascent-arc',
 'peri-target','lunar-grav','bias-side','polar-lock'].forEach(id =>
    document.getElementById(id).addEventListener('change', calcDay));

document.getElementById('in-window').addEventListener('input', e => {
    if (!selected || !selected.window) return;
    selectInWindow(selected.window, parseInt(e.target.value) / 100);
});

document.getElementById('tl-play').addEventListener('click', () => {
    if (!selected) return;
    anim.playing = !anim.playing;
    if (anim.playing && anim.tSec >= totalMissionSec()) { anim.tSec = 0; anim.trail = []; }
    document.getElementById('tl-play').textContent = anim.playing ? '⏸' : '▶';
});
document.getElementById('tl-slider').addEventListener('input', e => {
    if (!selected) return;
    anim.playing = false;
    document.getElementById('tl-play').textContent = '▶';
    anim.tSec = parseInt(e.target.value) / 1000 * totalMissionSec();
    anim.trail = [];
    applyAnimation();
});

document.querySelectorAll('.view-btn').forEach(btn => {
    btn.addEventListener('click', () => {
        document.querySelectorAll('.view-btn').forEach(b => b.classList.remove('active'));
        btn.classList.add('active');
        const v = { default: [24, 22, 34], top: [0, 52, 0.1], side: [46, 3, 0], earth: [3.2, 2.2, 3.2] }[btn.dataset.view];
        camera.position.set(v[0], v[1], v[2]);
        controls.target.set(0, 0, 0);
        controls.update();
    });
});

function onResize() {
    const w = viewport.clientWidth, h = viewport.clientHeight;
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
    renderer.setSize(w, h);
}
window.addEventListener('resize', onResize);

/* ---- 初始化 ---- */
document.getElementById('launch-date').value = new Date().toISOString().slice(0, 10);
onResize();
calcDay();
requestAnimationFrame(animate);
