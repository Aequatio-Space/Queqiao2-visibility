/* coast_tool/rockets.js — 运载火箭能力与历史任务验证案例（公开资料量级）
 *
 * 供 coast_tool 前端使用：
 *   - ROCKETS            火箭能力参数（LTO 运力锚点 / 上升段假设 / 滑行与窗口公开量级）
 *   - rocketCapability() C3 → 可用运力与余量（Tsiolkovsky 线性化；移植自
 *                        launch_window/launcher.py 的 c3_curve_linear）
 *   - widthFromMarginSec() Δv 偏航裕度 → 窗口宽度（移植自 launcher.py 的
 *                        width_from_margin）
 *   - VERIFY_CASES       三个历史任务验证案例（数据来自原 lunarlaunchwindow
 *                        仓库 verify_scenario.py / verify_results.txt）
 *
 * 说明：数值为公开资料量级，仅用于工程量级核对，不代表运载方性能承诺：
 *   - CZ-5：嫦娥五号任务论文《无人月球采样返回总体设计》——
 *     "长征五号运载火箭二级滑行时间 1000 s 时，地月转移轨道运载能力
 *      不小于 8200 kg"；窗口 50 min/天，"窄窗口多轨道"（5 条弹道）。
 *   - CZ-8：官方报道——地月转移轨道运力不低于 1.3 t；鹊桥二号 1.2 t，
 *     飞行 24 min 后星箭分离（直接入 200 × 42 万 km LTO）。
 */

/* 标准重力（km/s²），与 launcher.py 的 G0 一致 */
const G0_KM_S2 = 9.80665e-3;
/* 地球自转角速度（rad/s），与 launcher.py 的 OMEGA_EARTH 一致 */
const OMEGA_EARTH = 7.2921159e-5;

const ROCKETS = {
    cz5: {
        name: '长征五号 CZ-5',
        ascentTime: 740,        // s，实测器箭分离时序量级（融合引擎默认）
        ascentArc: 20.0,        // deg
        parkPeri: 200, parkApo: 380,   // km，200 × 380 km 停泊轨道
        mRefKg: 8200, c3Ref: -2.0, isp: 442,
        payloadKg: 8200,        // CE-5/6 探测器质量
        dvMarginMps: 30,        // 射向偏航 Δv 预算（工程量级假设）
        windowNote: '公开：50 min/天 · 5 条弹道（窄窗口多轨道）',
        source: 'CE-5 论文：二级滑行 1000 s 时 LTO 运力 ≥ 8200 kg'
    },
    cz8: {
        name: '长征八号 CZ-8',
        ascentTime: 800,        // s（公开未给精确值，取与验证脚本一致的量级）
        ascentArc: 20.0,
        parkPeri: 200, parkApo: 380,
        mRefKg: 1300, c3Ref: -2.0, isp: 442,
        payloadKg: 1200,        // 鹊桥二号质量
        dvMarginMps: 25,
        windowNote: '公开：3 天 × 2 轨道（鹊桥二号）',
        source: '官方：LTO 运力 > 1.3 t；鹊桥二号 1.2 t，T+24 min 星箭分离'
    }
};

/* C3 → 可用运力（kg）：Tsiolkovsky 线性化（launcher.c3_curve_linear）。
 * Δv(C3) = sqrt(C3 + 2μ/r) − sqrt(μ/r)
 * m_avail = m_ref · exp(−(Δv(C3) − Δv(C3_ref)) / (Isp·g0))
 */
function rocketCapability(rocket, c3, hParkKm, iParkDeg, siteLatDeg) {
    const r = R_EARTH + hParkKm;
    const muOverR = MU / r;
    const dvOf = c => Math.sqrt(Math.max(0, c) + 2 * muOverR) - Math.sqrt(muOverR);
    const dvRef = dvOf(rocket.c3Ref);
    const mAvail = rocket.mRefKg * Math.exp(-(dvOf(c3) - dvRef) / (rocket.isp * G0_KM_S2));
    const margin = mAvail - rocket.payloadKg;
    const vPark = Math.sqrt(MU / r);
    const width = widthFromMarginSec(iParkDeg, siteLatDeg, rocket.dvMarginMps, vPark);
    return {
        mAvail, margin, vPark,
        dvTli: dvOf(c3),
        widthS: width.widthS, widthFeasible: width.feasible
    };
}

/* Δv 偏航裕度 → 窗口半宽/全宽（launcher.width_from_margin）：
 *   Δv_yaw ≈ 2·v_park·sin(δΨ/2)，δΨ ≈ ω_E·cos(φ)·sin(i)·δt
 *   δt_max = 2·asin(Δv/(2 v_park)) / (ω_E·cos(φ)·sin(i))，width = 2·δt_max
 * 仅含射向偏航裕度，不含"变滑行时间"带来的额外窗口增益。
 */
function widthFromMarginSec(iParkDeg, siteLatDeg, dvMarginMps, vParkKms) {
    const i = iParkDeg * DEG, phi = siteLatDeg * DEG;
    const dv = dvMarginMps * 1.0e-3;
    if (dv <= 0) return { widthS: 0, feasible: true };
    if (Math.abs(Math.sin(phi)) > Math.abs(Math.sin(i)) + 1e-12) {
        return { widthS: 0, feasible: false };
    }
    const sinI = Math.sin(i), cosPhi = Math.cos(phi);
    if (sinI * cosPhi <= 0) return { widthS: 0, feasible: false };
    const arg = Math.max(-1, Math.min(1, dv / (2 * vParkKms)));
    const dPsi = 2 * Math.asin(arg);
    const dtMax = dPsi / (OMEGA_EARTH * cosPhi * sinI);
    return { widthS: 2 * dtMax, feasible: true };
}

/* 三个历史任务验证案例（原 lunarlaunchwindow 仓库 verify_scenario.py / verify_results.txt）。
 * 参考解 ref 为该仓库 Model B（LTO 面 = span{发射场@发射, 月球@到达}，切向 TLI 闭合）；
 * 实测 meas 为 Bill Gray Find_Orb 光学定轨；公开 public 为任务公开飞行剖面/窗口。 */
const VERIFY_CASES = {
    ce5: {
        name: '嫦娥五号 Chang\'e 5',
        rocket: 'cz5',
        launchUtc: '2020-11-23T20:30:22Z',
        loiUtc: '2020-11-28T12:58:00Z',
        transferH: 112.5,
        ref: { inc: 21.51, az: 99.3, coastMin: 21.7, coast2Min: 110.1, apoKm: 404496 },
        meas: { inc: 21.298, incSource: 'Find_Orb 实测' },
        public: {
            coastS: 1000,
            coastLabel: '二级滑行设计值 1000 s（≈16.7 min）',
            sep: '器箭分离 ~T+2200 s',
            lto: '200 × 410000 km',
            window: '3 天 × 50 min（5 条弹道）'
        },
        knownError: '计算滑行第一支 21.7 min 比二级滑行设计值 1000 s 偏长 ~5 min'
    },
    qq2: {
        name: '鹊桥二号 QueQiao-2',
        rocket: 'cz8',
        launchUtc: '2024-03-20T00:31:28Z',
        loiUtc: '2024-03-24T16:46:00Z',
        transferH: 112.2,
        ref: { inc: 22.93, az: 102.3, coastMin: 11.4, coast2Min: 99.8, apoKm: 405845 },
        meas: { inc: 21.976, incSource: 'Find_Orb 实测' },
        public: {
            coastS: 640,
            coastLabel: '由器箭分离 T+24 min 推得 ≈10.7 min（推算）',
            sep: '器箭分离 T+24 min',
            lto: '200 × 420000 km',
            window: '3 天 × 2 轨道'
        },
        knownError: '计算滑行 11.4 min 与公开推算 ~10.7 min 接近（+0.7 min）'
    },
    ce6: {
        name: '嫦娥六号 Chang\'e 6',
        rocket: 'cz5',
        launchUtc: '2024-05-03T09:27:29Z',
        loiUtc: '2024-05-08T02:12:00Z',
        transferH: 112.7,
        ref: { inc: 23.54, az: 103.5, coastMin: 23.1, coast2Min: 111.5, apoKm: 360224 },
        meas: null,
        public: {
            coastS: 1000,
            coastLabel: '二级滑行 ≤ ~1000 s（公开剖面量级）',
            sep: '器箭分离 ~T+37 min',
            lto: '200 × ~380000 km',
            window: '2 天 × 50 min（10 条轨道）'
        },
        knownError: '计算滑行 23.1 min 超出公开滑行量级（≤~16.7 min）约 6 min，'
            + '且近远地点到达（到达距离 366,585 km）对 LOI 时刻/月球引力敏感'
    }
};
