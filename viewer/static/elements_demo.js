'use strict';

/* ════════════════════════════════════════════════════════════════════
 * 轨道六根数演示 — UI / 3D 场景 / 文本渲染
 * 全部二体数学来自 ElementsCore（static/elements_core.js），
 * 加载时用 poliastro 基准（static/elements_selftest.js）对拍自检。
 * ════════════════════════════════════════════════════════════════════ */

const C = window.ElementsCore;
const X3 = window.ElementsCr3bp;   // 月球/CR3BP 模块
const $ = (id) => document.getElementById(id);

/* 全局错误钉在顶栏（含资源加载失败等 Event），便于无 console 环境定位 */
window.__errs = [];
function __pushErr(msg) {
  window.__errs.push(msg);
  const note = $('status-note');
  if (note) note.textContent = ' · 脚本错误：' + window.__errs.slice(-2).join(' | ');
  const badge = $('selftest');
  if (badge) {
    badge.textContent = '✗ 脚本错误（见顶栏）';
    badge.className = 'bad';
  }
}
window.addEventListener('error', (e) => {
  __pushErr(`${e.message} @${(e.filename || '').split('/').pop()}:${e.lineno}`);
}, true);
window.addEventListener('unhandledrejection', (e) => {
  __pushErr('Promise: ' + (e.reason && e.reason.message ? e.reason.message : e.reason));
});

if (!window.THREE) {
  __pushErr('three.js 未加载');
  $('status').innerHTML = '<h1 style="color:var(--bad)">three.js 加载失败</h1>' +
    '<p>请通过 Flask 入口访问（server.py → /elements_demo），或检查 static/vendor/three.min.js。</p>';
  throw new Error('THREE missing');
}

/* ── 常量与状态 ──────────────────────────────────────────────────── */
const KM = 1 / 1000;            // 场景单位：1 unit = 1000 km
const EARTH_R = C.R_EARTH * KM; // 6.378137 units
const EARTH_OMEGA = 7.2921159e-5; // rad/s，地球恒星自转角速度（与 elements_perturb.js 一致；core 未导出，本地定义）
/* 会合系旋转轴（场景坐标）：轨道面法向 e3 = Mᵀ·ẑ_orbit，即显示系轨道法向。
   L2 视角下会合系静态环（θ=0 锚定）按此轴实时刚性旋转，环心恒压在实时 L2 上；
   等价于 rotStateToDisplay(y,t) = R(轴, ω·t)·rotStateToDisplay(y,0)（数值验证残差 <0.1 km）。 */
const FRAME_AXIS = eciToScene(X3.orbitToDisp([0, 0, 1])).normalize();

/* 范围表见 BODY_RANGES（下方定义；RANGES 引用其 earth 条目） */

/* 真实 TLE 根数（2026-10-11 历元，NORAD 25544 / 48274）。
   a 由平均运动 n 反算；nu 由 M→E→ν 转换。
   TLE 验证见 pert-ctl 注释：ISS 30d Ω 漂移 模型 -4.98°/d vs 实测 -4.95°/d（误差<1%）。 */
const PRESETS = {
  iss:     { a: 6798.281, e: 0.0007, i: 51.6314, Om: 82.9019, w: 250.7442, nu: 109.3546 },
  css:     { a: 6764.043, e: 0.0001, i: 41.4699, Om: 333.4359, w: 31.4962, nu: 328.5877 },
  leo:     { a: 6778.137, e: 0.0006, i: 51.6,  Om: 120,  w: 30,  nu: 90 },
  sso:     { a: 7078.137, e: 0.0012, i: 98.2,  Om: 210,  w: 15,  nu: 33.7 },
  molniya: { a: 26562,    e: 0.74,   i: 63.4,  Om: 60,   w: 270, nu: 180 },
  geo:     { a: 42164,    e: 0.0004, i: 0.05,  Om: 0,    w: 10,  nu: 270 },
  /* 北斗 MEO：T=12 h 反解 a ≈ 26 610 km（高度 ≈ 20 230 km，约 2 万公里），i=55° 为北斗 MEO 星座标称倾角 */
  bds:     { a: Math.cbrt(C.MU * (12 * 3600 / (2 * Math.PI)) ** 2), e: 0.002, i: 55, Om: 0, w: 0, nu: 0 },
  tli:     { a: 195489,   e: 0.9664, i: 28.5,  Om: 10,   w: 180, nu: 0 },
};

/* ── 月球预设（二体根数型；CR3BP 型见 MOON_PRESETS_CR3BP） ──────
 * CLO：B10 §2.1.1，12 h 环月圆轨道（a 由 T=12h 反解，高度 ≈ 4405 km）。
 * ELFO：鹊桥二号使命轨道（S13 表 1，太阳同步冻结测控回归轨道，T 取 24 h）：
 *   冻结条件 e=√(1−5cos²i/3)（B10 式(3)）与 i=119.203°（S13 收敛值，逆行）
 *   → e=0.7767，近月点高度 ≈ 440 km，ω 冻结在 90°。
 * ELL：15×200 km 环月椭圆轨道（hp=15 km、ha=200 km 反解 a/e：
 *   a = (rp+ra)/2 ≈ 1844.9 km，e = (ra−rp)/(ra+rp) ≈ 0.0501）。 */
const MOON_PRESETS = {
  clo:  { a: Math.cbrt(X3.MU_MOON * (12 * 3600 / (2 * Math.PI)) ** 2), e: 0.001, i: 90, Om: 0, w: 0, nu: 0 },
  ell:  (() => { const rp = X3.R_MOON + 15, ra = X3.R_MOON + 200;
    return { a: (rp + ra) / 2, e: (ra - rp) / (ra + rp), i: 90, Om: 0, w: 0, nu: 0 }; })(),
  elfo: { a: Math.cbrt(X3.MU_MOON * (24 * 3600 / (2 * Math.PI)) ** 2), e: 0.7767, i: 119.203, Om: 0, w: 90, nu: 0 },
};

/* CR3BP 型预设（会合系无量纲初值在 X3.PRESET_IC，此处只存 UI 元数据） */
const MOON_PRESETS_CR3BP = {
  halo: { label: '鹊桥一号 Halo', preset: 'halo' },
  nrho: { label: 'NRHO 9:2', preset: 'nrho' },
  dro:  { label: 'DRO 共振轨道', preset: 'dro' },
};

/* ── 中心天体状态（选项卡） ── */
let body = 'earth';        // 'earth' | 'moon'
const BODY_RANGES = {
  earth: {
    a:  [6600, 200000, 1],
    e:  [0, 0.98, 0.001],
    i:  [0, 180, 0.1],
    Om: [0, 360, 0.5],
    w:  [0, 360, 0.5],
    nu: [0, 360, 0.1],
  },
  moon: {
    a:  [X3.R_MOON + 15, X3.R_MOON + 20000, 1],
    e:  [0, 0.98, 0.001],
    i:  [0, 180, 0.1],
    Om: [0, 360, 0.5],
    w:  [0, 360, 0.5],
    nu: [0, 360, 0.1],
  },
};
let cr3bpState = { active: false, preset: 'nrho', azKm: 13000, phase: 0, stationKeeping: false };
/* 当前月球预设键：applyMoonPreset 写入；手动拖根数滑块（非 ν）、点地球预设时清空。
   用于「ELFO 在摄动二体/CR3BP 下不再弹模型推荐警示」这类预设专属判定 */
let moonPresetKey = null;
/* CR3BP 三维周期轨道无法由开普勒六根数描述：禁用六根数滑块并标注，切回根数型预设即恢复 */
function setElSlidersEnabled(on) {
  document.querySelectorAll('.el-row').forEach((row) => {
    row.querySelector('input[type="range"]').disabled = !on;
    row.querySelector('input[type="number"]').disabled = !on;
  });
  /* 用 .show 类切换：.hintbox 基类自带 display:none，清内联样式不够 */
  $('hint-cr3bp').classList.toggle('show', !on);
}

/* RANGES 必须在 BODY_RANGES 之后（选项卡切换时经 applyBodyRanges 整体重赋值并刷到滑块） */
let RANGES = BODY_RANGES.earth;

let els = { ...PRESETS.molniya };   // 角度一律用「度」存；喂给 core 前转弧度
let prop = null;                     // 当前传播器实例
let playing = false;
let speedX = 600;                  // 仿真秒/真实秒（对应下拉框 10 min/s）
let tSincePeri = 0;                  // 自最近过近地点起算的秒数（播放相位）
/* UTC 历元：tSincePeri=0（近地点）对应的真实 UTC 毫秒。仿真 UTC = epochUTCms + tSincePeri×1000。
   播放时 tSincePeri 按 dt×speedX 推进 ⇒ 仿真时钟与真实墙钟同步推进（倍率 ×speedX）；
   倍率 ×1（实时）时仿真 UTC 与实际 UTC 完全一致。null = 尚未锚定（首次 rebuild 时对齐到现在）。 */
let epochUTCms = null;
function syncEpochToNow() { epochUTCms = Date.now() - tSincePeri * 1000; }
function fmtUTC(ms) {
  const d = new Date(ms);
  const p = (n) => String(n).padStart(2, '0');
  return `${d.getUTCFullYear()}-${p(d.getUTCMonth() + 1)}-${p(d.getUTCDate())} ` +
         `${p(d.getUTCHours())}:${p(d.getUTCMinutes())}:${p(d.getUTCSeconds())}`;
}
function updateUtcClock() {
  if (epochUTCms === null) return;
  $('utc-sim').textContent = `UTC ${fmtUTC(epochUTCms + tSincePeri * 1000)}`;
  $('utc-real').textContent = `实际 ${fmtUTC(Date.now())} · ×${speedX}`;
}
let texFailed = false;

const elsRad = () => ({
  a: els.a, e: els.e,
  i: els.i * C.D2R, Om: els.Om * C.D2R, w: els.w * C.D2R,
});

/* ── 数值格式化 ──────────────────────────────────────────────────── */
function fmt(x, dp) {
  if (x === null || !isFinite(x)) return '—';
  if (Math.abs(x) < 1e-9) return '0';
  const a = Math.abs(x);
  if (dp === undefined) {
    dp = a < 1 ? 6 : a < 10 ? 4 : a < 100 ? 3 : a < 1000 ? 2 : a < 1e6 ? 1 : 0;
  }
  return x.toLocaleString('zh-CN', { minimumFractionDigits: dp, maximumFractionDigits: dp });
}
function fmtVec(v, unit) {
  return `(${fmt(v[0])}, ${fmt(v[1])}, ${fmt(v[2])})${unit ? ' ' + unit : ''}`;
}
function fmtDeg(x, dp) { return `${fmt(x, dp === undefined ? 2 : dp)}°`; }
function fmtT(sec) {
  if (!isFinite(sec)) return '—';
  if (sec < 120) return `${fmt(sec, 1)} s`;
  if (sec < 7200) return `${fmt(sec / 60, 1)} min`;
  return `${fmt(sec / 3600, 3)} h`;
}

/* ── 3D 场景 ─────────────────────────────────────────────────────── */
const stage = $('stage');
const renderer = new THREE.WebGLRenderer({ antialias: true, logarithmicDepthBuffer: true });
/* 探测 GL 渲染器：纯软件渲染（SwiftShader / llvmpipe）下降低像素比与帧率，
   防止无 GPU 环境（含 headless 自动化）下软件光栅化吃满 CPU */
let _glRenderer = '';
try {
  const _gl = renderer.getContext();
  const _ext = _gl.getExtension('WEBGL_debug_renderer_info');
  _glRenderer = _ext ? _gl.getParameter(_ext.UNMASKED_RENDERER_WEBGL) : _gl.getParameter(_gl.VERSION);
} catch (_) {}
const SOFT_GL = /swiftshader|llvmpipe|software|Microsoft Basic/i.test(_glRenderer || '');
renderer.setPixelRatio(SOFT_GL ? 1 : Math.min(window.devicePixelRatio || 1, 2));
renderer.setSize(window.innerWidth, window.innerHeight);
/* FatLine 粗线：uResolution 用绘图缓冲（物理像素），宽按 CSS 像素 ×pixelRatio */
THREE.setFatLinePixelRatio(renderer.getPixelRatio());
THREE.setFatLineResolution(renderer.domElement.width, renderer.domElement.height);
stage.appendChild(renderer.domElement);

const scene = new THREE.Scene();
scene.background = new THREE.Color(0x070b14);
const camera = new THREE.PerspectiveCamera(50, window.innerWidth / window.innerHeight, 0.05, 4000);

/* ECI (x→♈, y, z→北极) → three.js 场景 (x, y↑, z)，保持右手系 */
function eciToScene(v) { return new THREE.Vector3(v[0] * KM, v[2] * KM, -v[1] * KM); }

/* 贴图加载器必须放在首次调用之前（月球贴图在下方立即发起加载）。
   曾把 loadTex 定义放在文件后部：调用落在 const texLoader 暂时性死区，
   ReferenceError 被 Promise 构造器转成 rejection，.catch 误当成"网络失败"，
   于是程序化替代 + 提示请走 http，而真实图片请求从未发出。 */
const texLoader = new THREE.TextureLoader();
texLoader.setCrossOrigin('anonymous');
function loadTex(url) {
  return new Promise((resolve, reject) => {
    texLoader.load(url, resolve, undefined, () => reject(new Error('load fail: ' + url)));
  });
}

/* — 地球：昼/夜贴图 shader（贴图失败回退程序化） — */
const SUN_DIR = new THREE.Vector3(1, 0.3, 0.45).normalize();
const earthUniforms = {
  dayTex:   { value: null },
  nightTex: { value: null },
  sunDir:   { value: SUN_DIR },
};
const earthMat = new THREE.ShaderMaterial({
  uniforms: earthUniforms,
  vertexShader: `
    varying vec3 vNormalW; varying vec2 vUv;
    void main() {
      vNormalW = normalize((modelMatrix * vec4(normal, 0.0)).xyz);
      vUv = uv;
      gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
    }`,
  fragmentShader: `
    uniform sampler2D dayTex; uniform sampler2D nightTex; uniform vec3 sunDir;
    varying vec3 vNormalW; varying vec2 vUv;
    void main() {
      float d = dot(normalize(vNormalW), normalize(sunDir));
      float dayMix = smoothstep(-0.15, 0.15, d);
      vec3 day = texture2D(dayTex, vUv).rgb * (0.22 + 0.95 * max(d, 0.0));
      vec3 night = texture2D(nightTex, vUv).rgb * 1.45 + vec3(0.008, 0.012, 0.022);
      gl_FragColor = vec4(mix(night, day, dayMix), 1.0);
    }`,
});
const earth = new THREE.Mesh(new THREE.SphereGeometry(EARTH_R, 96, 64), earthMat);
scene.add(earth);

/* — 月球（地球选项卡：绕地球圆化轨道，SPICE 历元锁定；月球选项卡：中心天体） — */
const MOON_R = X3.R_MOON * KM;               // 1.7374 units
const moonUniforms = {
  dayTex: { value: null },
  sunDir: { value: SUN_DIR },
};
const moonMat = new THREE.ShaderMaterial({
  uniforms: moonUniforms,
  vertexShader: earthMat.vertexShader,
  /* 月球专用：无夜面灯光，背阳面纯黑。
     不能复用地球 fragmentShader——其 nightTex 项在 moonUniforms 中未声明，
     采样器默认绑到 0 号纹理单元（昼图），暗面被错误点亮。 */
  fragmentShader: `
    uniform sampler2D dayTex; uniform vec3 sunDir;
    varying vec3 vNormalW; varying vec2 vUv;
    void main() {
      float d = dot(normalize(vNormalW), normalize(sunDir));
      float dayMix = smoothstep(-0.12, 0.12, d);
      vec3 day = texture2D(dayTex, vUv).rgb * (0.25 + 0.95 * max(d, 0.0));
      gl_FragColor = vec4(day * dayMix, 1.0);
    }`,
});
const moonMesh = new THREE.Mesh(new THREE.SphereGeometry(MOON_R, 64, 48), moonMat);
moonMesh.name = 'moon';
moonMesh.visible = false;
scene.add(moonMesh);
/* 程序化月球回退贴图（灰色陨石坑风格）。顶层定义：加载失败与渲染拦截两条
   回退路径（loadTex .catch / tick 的 WebGL 抛错）都要用它。 */
let moonTexFailed = false;
function proceduralMoonTexture() {
  const cv = document.createElement('canvas');
  cv.width = 1024; cv.height = 512;
  const ctx = cv.getContext('2d');
  ctx.fillStyle = '#8a8a8a'; ctx.fillRect(0, 0, cv.width, cv.height);
  for (let k = 0; k < 260; k++) {
    const r = 2 + Math.random() * 26;
    const g = ctx.createRadialGradient(0, 0, 0, 0, 0, r);
    const shade = 100 + Math.floor(Math.random() * 60);
    g.addColorStop(0, `rgba(${shade},${shade},${shade},0.85)`);
    g.addColorStop(0.75, `rgba(${shade - 30},${shade - 30},${shade - 30},0.5)`);
    g.addColorStop(1, 'rgba(200,200,200,0.35)');
    ctx.save();
    ctx.translate(Math.random() * cv.width, cv.height * (0.06 + 0.88 * Math.random()));
    ctx.fillStyle = g;
    ctx.beginPath(); ctx.arc(0, 0, r, 0, Math.PI * 2); ctx.fill();
    ctx.restore();
  }
  return new THREE.CanvasTexture(cv);
}
function useProceduralMoon(reason) {
  moonUniforms.dayTex.value = proceduralMoonTexture();
  moonTexFailed = true;
  const n = $('status-note');
  if (n) n.textContent = `（月球贴图不可用，已用程序化替代：${reason}）`;
}
loadTex('../../assets/moon_surface.jpg')
  .then((tex) => { moonUniforms.dayTex.value = tex; })
  .catch((err) => {
    /* 失败必须可见：静默回退会让用户误以为贴图本身错了（此前"月球没有正确贴图"即此静默回退所致） */
    useProceduralMoon((err && err.message) || '加载失败');
  });

/* 地球选项卡：月球轨道圈（圆化模型，仅视觉参考） */
let moonOrbitRing = null;
function buildMoonOrbitRing() {
  const pts = [];
  const O = window.ElementsPerturb.MOON_ORBIT;
  for (let k = 0; k <= 128; k++) {
    const th = O.theta0 + (k / 128) * C.TWO_PI;
    const c = Math.cos(th), s = Math.sin(th);
    pts.push(new THREE.Vector3(
      O.a * (c * O.node[0] + s * O.eperp[0]) * KM,
      O.a * (c * O.node[2] + s * O.eperp[2]) * KM,
      -O.a * (c * O.node[1] + s * O.eperp[1]) * KM));
  }
  const ring = new THREE.FatLine(
    new THREE.BufferGeometry().setFromPoints(pts),
    { color: 0x9aa7b8, transparent: true, opacity: 0.35 }
  );
  return ring;
}

/* 月球选项卡：地球迷你球（贴图复用 earthUniforms 的纹理） */
const EARTH_MINI_R = EARTH_R;   // 同比例（真实比例 6378/1737≈3.67）
const earthMini = new THREE.Mesh(new THREE.SphereGeometry(EARTH_MINI_R, 48, 32), earthMat);
earthMini.name = 'earth-mini';
earthMini.visible = false;
scene.add(earthMini);

/* L1/L2 标记与连线（月球选项卡，CR3BP 模式显示） */
const lpointGroup = new THREE.Group();
scene.add(lpointGroup);
function buildLpointMarks() {
  lpointGroup.clear();
  for (const [name, color] of [['L1', 0x7dffa8], ['L2', 0xffa552]]) {
    const mk = new THREE.Mesh(new THREE.SphereGeometry(0.55, 16, 12),
      new THREE.MeshBasicMaterial({ color }));
    mk.name = name;
    lpointGroup.add(mk);
    const lb = makeLabel(name, name === 'L1' ? '#7dffa8' : '#ffa552');
    lb.name = name + '-label';
    lpointGroup.add(lb);
  }
}
buildLpointMarks();
lpointGroup.visible = false;

/* 大气辉光（BackSide 加色） */
const glow = new THREE.Mesh(
  new THREE.SphereGeometry(EARTH_R * 1.05, 64, 48),
  new THREE.ShaderMaterial({
    vertexShader: `
      varying vec3 vNormal;
      void main() {
        vNormal = normalize(normalMatrix * normal);
        gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
      }`,
    fragmentShader: `
      varying vec3 vNormal;
      void main() {
        float i = pow(max(0.0, 0.68 - dot(vNormal, vec3(0.0, 0.0, 1.0))), 4.0);
        gl_FragColor = vec4(0.30, 0.55, 1.0, 1.0) * i;
      }`,
    side: THREE.BackSide, blending: THREE.AdditiveBlending,
    transparent: true, depthWrite: false,
  })
);
scene.add(glow);

/* 程序化贴图回退（file:// 直开时 WebGL 取本地贴图会被浏览器安全策略拦截） */
function proceduralEarthTextures() {
  function canvasTex(draw) {
    const cv = document.createElement('canvas');
    cv.width = 1024; cv.height = 512;
    draw(cv.getContext('2d'), cv.width, cv.height);
    return new THREE.CanvasTexture(cv);
  }
  const day = canvasTex((ctx, w, h) => {
    const g = ctx.createLinearGradient(0, 0, 0, h);
    g.addColorStop(0, '#2b4d6e'); g.addColorStop(0.5, '#17548e'); g.addColorStop(1, '#2b4d6e');
    ctx.fillStyle = g; ctx.fillRect(0, 0, w, h);
    const land = ['#3f7a4f', '#5d8a4a', '#8a9a5b', '#b5a642', '#7a6d4f'];
    for (let k = 0; k < 16; k++) {
      ctx.fillStyle = land[k % land.length];
      ctx.globalAlpha = 0.85;
      const cx = Math.random() * w, cy = h * (0.18 + 0.64 * Math.random());
      ctx.beginPath();
      ctx.ellipse(cx, cy, 40 + Math.random() * 130, 25 + Math.random() * 70,
                  Math.random() * Math.PI, 0, Math.PI * 2);
      ctx.fill();
    }
    ctx.globalAlpha = 0.9; ctx.fillStyle = '#e8eef5';
    ctx.fillRect(0, 0, w, h * 0.055); ctx.fillRect(0, h * 0.945, w, h * 0.055);
    ctx.globalAlpha = 1;
  });
  const night = canvasTex((ctx, w, h) => {
    ctx.fillStyle = '#04060d'; ctx.fillRect(0, 0, w, h);
    for (let k = 0; k < 900; k++) {
      ctx.fillStyle = Math.random() < 0.7 ? '#ffd28a' : '#ffab6b';
      ctx.globalAlpha = 0.25 + 0.65 * Math.random();
      const r = 0.4 + 1.4 * Math.random();
      ctx.beginPath();
      ctx.arc(Math.random() * w, h * (0.1 + 0.8 * Math.random()), r, 0, Math.PI * 2);
      ctx.fill();
    }
    ctx.globalAlpha = 1;
  });
  return { day, night };
}
function useProceduralEarth(reason) {
  const t = proceduralEarthTextures();
  earthUniforms.dayTex.value = t.day;
  earthUniforms.nightTex.value = t.night;
  texFailed = true;
  $('status-note').textContent = `（地球贴图加载失败，已用程序化替代：${reason}）`;
}

Promise.all([loadTex('../../assets/earth_day.jpg'), loadTex('../../assets/earth_night.jpg')])
  .then(([day, night]) => {
    if (texFailed) return; // 已回退则不再切换
    earthUniforms.dayTex.value = day;
    earthUniforms.nightTex.value = night;
  })
  .catch((err) => useProceduralEarth('网络/跨源'));

/* — 星空背景 — */
{
  const n = 1300, pos = new Float32Array(n * 3);
  for (let k = 0; k < n; k++) {
    const v = new THREE.Vector3().randomDirection
      ? new THREE.Vector3().randomDirection()
      : new THREE.Vector3(Math.random() * 2 - 1, Math.random() * 2 - 1, Math.random() * 2 - 1).normalize();
    v.multiplyScalar(1200);
    pos[k * 3] = v.x; pos[k * 3 + 1] = v.y; pos[k * 3 + 2] = v.z;
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
  scene.add(new THREE.Points(g, new THREE.PointsMaterial({
    size: 1.6, sizeAttenuation: false, color: 0xaebdd0, transparent: true, opacity: 0.75,
  })));
}

/* — 可重建对象（轨道线 / 网格 / 节线 / 标记 / 标签 / 轴） — */
let orbitLine = null, grid = null, nodeLine = null;
let periMark = null, apoMark = null;
/* 当前模式是否有拱点标记（近/远月点）：由 rebuild 设置，applyDisplay 读取——
   无拱点时标签开关不得把 peri/apo 标签强制点亮（会停在默认原点 = 月心） */
let apsisAvailable = false;
const labels = {};   // name → sprite
let gridR = 30;

function makeLabel(text, color) {
  /* 4× 超采样栅格化：sprite 世界尺寸不变、纹理分辨率 ×4，
     消除标签放大后的马赛克（点标记为几何体，天生矢量清晰） */
  const SS = 4;
  const cv = document.createElement('canvas');
  const ctx = cv.getContext('2d');
  const font = `600 ${42 * SS}px -apple-system, "PingFang SC", sans-serif`;
  ctx.font = font;
  cv.width = Math.ceil(ctx.measureText(text).width) + 26 * SS;
  cv.height = 64 * SS;
  ctx.font = font;
  ctx.textBaseline = 'middle'; ctx.textAlign = 'center';
  ctx.lineWidth = 8 * SS; ctx.strokeStyle = 'rgba(5,8,15,0.9)';
  ctx.strokeText(text, cv.width / 2, cv.height / 2);
  ctx.fillStyle = color;
  ctx.fillText(text, cv.width / 2, cv.height / 2);
  const tex = new THREE.CanvasTexture(cv);
  const sp = new THREE.Sprite(new THREE.SpriteMaterial({
    map: tex, transparent: true, depthTest: false,
  }));
  sp.renderOrder = 20;
  sp.userData.text = text;   // 供诊断/测试按文字识别标签
  sp.userData.aspect = cv.width / cv.height;
  scene.add(sp);
  return sp;
}
/* 标签与点标记/航天器同一缩放律：世界尺寸 ∝ 相机视距，屏占比跨轨道尺度恒定、
   相互比例固定（「绑定」效果）；1.7 使文字视大小保持 DRO/NRHO 下已认可的水平 */
function labelScale() { return zoomK() * 1.7; }
function placeLabel(sp, scenePos, dy) {
  const s = labelScale();
  sp.scale.set(s * sp.userData.aspect, s, 1);
  sp.position.copy(scenePos);
  sp.position.y += dy !== undefined ? dy * s : s * 0.7;
}

/* ECI 坐标轴（春分点方向是 Ω 的量基准） */
const axesGroup = new THREE.Group();
scene.add(axesGroup);
{
  const mk = (dir, len, color) => {
    const a = new THREE.ArrowHelper(dir, new THREE.Vector3(0, 0, 0), len, color, len * 0.12, len * 0.06);
    axesGroup.add(a);
    return a;
  };
  mk(new THREE.Vector3(1, 0, 0), 9, 0x66ddaa);            // X → ♈
  mk(new THREE.Vector3(0, 0, -1), 7, 0x55637a);           // Y
  mk(new THREE.Vector3(0, 1, 0), 9, 0x7dffa8);            // Z → 北极
}

/* 航天器 + 速度/地心矢量 + 当地水平盘 + 位置虚线 */
const scGroup = new THREE.Group();
scene.add(scGroup);
const scMesh = new THREE.Mesh(
  new THREE.SphereGeometry(0.42, 24, 18),
  new THREE.MeshBasicMaterial({ color: 0x4ec9ff })
);
scGroup.add(scMesh);
const scHalo = new THREE.Mesh(
  new THREE.SphereGeometry(0.75, 16, 12),
  new THREE.MeshBasicMaterial({ color: 0x4ec9ff, transparent: true, opacity: 0.22 })
);
scGroup.add(scHalo);

const velArrow = new THREE.ArrowHelper(new THREE.Vector3(1, 0, 0), new THREE.Vector3(), 1, 0xffd166, 0.3, 0.16);
scene.add(velArrow);
const rArrow = new THREE.ArrowHelper(new THREE.Vector3(1, 0, 0), new THREE.Vector3(), 1, 0x9aa7b8, 0.3, 0.16);
scene.add(rArrow);
const disc = new THREE.Mesh(
  new THREE.RingGeometry(0.82, 1.0, 48),
  new THREE.MeshBasicMaterial({ color: 0x4ec9ff, transparent: true, opacity: 0.32, side: THREE.DoubleSide })
);
scene.add(disc);

const posLine = new THREE.FatLine(
  new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(), new THREE.Vector3()]),
  { color: 0x9aa7b8, dashed: true, dashSize: 0.55, gapSize: 0.38, transparent: true, opacity: 0.7 }
);
scene.add(posLine);

/* ── CR3BP 受力分析：月球引力（白）/ 地球引力（蓝）/ 合力（红）三箭头 + 标签 ──
   仅三体模式显示（row-forces 由 rebuild 按 isCr3bp 显隐）。
   物理量：a_月 = −μ☾·r/|r|³、a_⊕ = μ⊕·(r⊕−r)/|r⊕−r|³（牛顿直接项，月心显示系），
   a_合 = a_月 + a_⊕。方向严格真实；长度动态范围跨轨道达数百倍（NRHO 近/远月点
   月引力 ~480×），故按当帧最大力归一并取 √ 压缩：L = L0·√(a/a_max)。 */
const forceGroup = new THREE.Group();
scene.add(forceGroup);
const fMoonArr = new THREE.ArrowHelper(new THREE.Vector3(1, 0, 0), new THREE.Vector3(), 1, 0xe8e8e8, 0.3, 0.16);
const fEarthArr = new THREE.ArrowHelper(new THREE.Vector3(1, 0, 0), new THREE.Vector3(), 1, 0x6fa8ff, 0.3, 0.16);
const fTotArr = new THREE.ArrowHelper(new THREE.Vector3(1, 0, 0), new THREE.Vector3(), 1, 0xff6b6b, 0.3, 0.16);
forceGroup.add(fMoonArr, fEarthArr, fTotArr);
/* makeLabel 直接挂 scene，再 add 进 forceGroup 即完成换父（组无变换，场景坐标一致） */
const fLabels = {
  moon: makeLabel('F☾', '#e8e8e8'),
  earth: makeLabel('F⊕', '#6fa8ff'),
  tot: makeLabel('F合', '#ff6b6b'),
};
for (const k in fLabels) forceGroup.add(fLabels[k]);
forceGroup.visible = false;
window.__forceGroup = forceGroup;   // 供诊断/测试读取（无副作用）

/* 每帧更新受力箭头（仅 CR3BP；p = 航天器场景位置，zk = zoomK） */
function updateForces(p, st, zk) {
  const r = st.rE;                                  /* 航天器月心位置 km（显示系） */
  const rE = X3.earthPosMoonTab(tSincePeri);        /* 地球月心位置 km（同系） */
  const rM = Math.hypot(r[0], r[1], r[2]);
  const dE = [rE[0] - r[0], rE[1] - r[1], rE[2] - r[2]];
  const rD = Math.hypot(dE[0], dE[1], dE[2]);
  const aM = X3.MU_MOON / (rM * rM);                /* 月球引力加速度 km/s² */
  const aE = C.MU / (rD * rD);                      /* 地球引力加速度 km/s² */
  const vM = [-r[0] / rM, -r[1] / rM, -r[2] / rM];  /* 指向月球 */
  const vE = [dE[0] / rD, dE[1] / rD, dE[2] / rD];  /* 指向地球 */
  const aS = [aM * vM[0] + aE * vE[0], aM * vM[1] + aE * vE[1], aM * vM[2] + aE * vE[2]];
  const mS = Math.hypot(aS[0], aS[1], aS[2]);
  const aMax = Math.max(aM, aE, mS, 1e-30);
  const L0 = 3.2 * zk;                              /* 最大箭头长度（与速度箭头同一 zk 缩放律） */
  const setF = (arr, lb, vec, mag) => {
    const dir = eciToScene(vec).normalize();
    const len = Math.max(L0 * Math.sqrt(mag / aMax), 0.18 * zk);   /* √ 压缩 + 可见下限 */
    arr.position.copy(p);
    arr.setDirection(dir);
    arr.setLength(len, Math.max(0.12 * zk, len * 0.22), Math.max(0.07 * zk, len * 0.11));
    placeLabel(lb, p.clone().addScaledVector(dir, len), 0.25);
  };
  setF(fMoonArr, fLabels.moon, vM, aM);
  setF(fEarthArr, fLabels.earth, vE, aE);
  setF(fTotArr, fLabels.tot, aS, mS);
}

/* — 相机控制：绑定在 window 级 + 画布命中检测，规避任何子元素/层叠吞事件 — */
const camCtl = { az: -1.05, el: 0.42, dist: 60, target: new THREE.Vector3() };
window.__camCtl = camCtl;   // 供诊断/测试读取相机状态（无副作用）
window.__cam = camera;      // 供诊断/测试读取相机对象（无副作用）
window.__scene = scene;     // 供诊断/测试读取场景对象（无副作用）
window.__scGroup = scGroup; // 供诊断/测试读取航天器组位置（无副作用）
window.__cr3bpState = cr3bpState; // 供诊断/测试读写 CR3BP 预设状态（同引用，改后需触发 rebuild）
window.__applyMoonPreset = applyMoonPreset; // 供诊断/测试套用月球预设（等同点击预设按钮）
let fpv = false;                     // 第一视角（航天器视角）模式
const fpvCam = { az: 0, el: 0 };     // 第一视角下的视线偏移（拖拽微调）
let scState = null;                  // 每帧由 updateScene 刷新：{ p, vScene, rScene }
function setFov(f) { camera.fov = Math.min(95, Math.max(25, f)); camera.updateProjectionMatrix(); }
function applyCamera() {
  if (fpv && scState) {
    /* 第一视角：相机贴在航天器上方一点，沿速度方向看向前方 */
    const { p, vScene, rScene } = scState;
    const fwd = vScene.clone().applyAxisAngle(rScene, fpvCam.az);
    const right = new THREE.Vector3().crossVectors(fwd, rScene).normalize();
    fwd.applyAxisAngle(right, fpvCam.el);
    camera.up.copy(rScene);
    camera.position.copy(p).addScaledVector(rScene, 0.6);
    camera.lookAt(camera.position.clone().add(fwd));
    return;
  }
  camera.up.set(0, 1, 0);
  const { az, el, dist, target } = camCtl;
  /* 俯视·地球居中：目标每帧锁实时地球位置（地球绕原点公转），惯性系不共转——
     地球静止于画面中心，月球与航天器环绕地球运动 */
  if (topView) target.copy(earthMini.position);
  camera.position.set(
    target.x + dist * Math.cos(el) * Math.cos(az),
    target.y + dist * Math.sin(el),
    target.z + dist * Math.cos(el) * Math.sin(az)
  );
  /* L2 会合系视角：相机绕 FRAME_AXIS 随会合系共转（与会合系环
     同一 ω_EM·t 角）。位置与 up 必须一起转——只转位置不转 up 时，绕倾斜轴（轨道面
     法向）的共转会产生滚转角变化，使 L2 等会合系固定点在屏幕上漂移（2026-10-14
     验证发现）。拖动期间用冻结值，松手后无跳变续转。 */
  if (l2View) {
    const spin = l2SpinFrozen !== null ? l2SpinFrozen : l2SpinBase + X3.OMEGA_EM * tSincePeri;
    const qL2 = new THREE.Quaternion().setFromAxisAngle(FRAME_AXIS, spin);
    camera.position.sub(target).applyQuaternion(qL2).add(target);
    camera.up.applyQuaternion(qL2);
  }
  camera.lookAt(target);
}
/* 命中检测：该坐标最顶层的元素是否就是 3D 画布本身。
   canvas 全屏铺底，控制面板 / 输出卡片只是浮在它上面，
   绝不能只用 canvas 的几何矩形判断——否则在滑块、滚动条上的
   pointerdown / wheel 会被误当成画布事件吞掉（preventDefault +
   setPointerCapture 会让滑块无法拖动、面板无法滚动）。 */
function overCanvas(clientX, clientY) {
  const top = document.elementFromPoint(clientX, clientY);
  return top === renderer.domElement;
}
{
  const cvs = renderer.domElement;
  const pointers = new Map();
  let pinchDist = 0;
  cvs.style.touchAction = 'none';
  /* 记录活动指针是否在画布内开始（pointerdown 必须在画布上才算数） */
  window.addEventListener('pointerdown', (ev) => {
    if (!overCanvas(ev.clientX, ev.clientY)) return;
    pointers.set(ev.pointerId, { x: ev.clientX, y: ev.clientY });
    try { cvs.setPointerCapture(ev.pointerId); } catch (_) { /* 某些 Chrome 版本对鼠标会抛，忽略 */ }
    /* L2 共转视角：拖动开始即冻结共转角，避免视角与用户抢鼠标（俯视·地球居中
       为惯性系视角，无共转，无需冻结） */
    if (pointers.size === 1 && l2View) l2SpinFrozen = l2SpinBase + X3.OMEGA_EM * tSincePeri;
    if (pointers.size === 2) {
      const p = [...pointers.values()];
      pinchDist = Math.hypot(p[0].x - p[1].x, p[0].y - p[1].y);
    }
    ev.preventDefault();
  });
  /* pointermove 绑在 window：即使 capture 失败/事件被中间层截获，拖拽仍连续 */
  window.addEventListener('pointermove', (ev) => {
    if (!pointers.has(ev.pointerId)) return;
    const prev = pointers.get(ev.pointerId);
    pointers.set(ev.pointerId, { x: ev.clientX, y: ev.clientY });
    if (pointers.size === 1) {
      /* 概览模式拖动相机球坐标；第一视角模式拖动视线偏移（转头） */
      const t = fpv ? fpvCam : camCtl;
      t.az += (ev.clientX - prev.x) * 0.0052;
      t.el = Math.min(1.45, Math.max(-1.45, t.el + (ev.clientY - prev.y) * 0.0052));
    } else if (pointers.size === 2) {
      const p = [...pointers.values()];
      const d = Math.hypot(p[0].x - p[1].x, p[0].y - p[1].y);
      if (pinchDist > 0) {
        if (fpv) setFov(camera.fov * pinchDist / d);
        else camCtl.dist = clampDist(camCtl.dist * pinchDist / d);
      }
      pinchDist = d;
    }
  });
  const release = (ev) => {
    pointers.delete(ev.pointerId);
    if (pointers.size < 2) pinchDist = 0;
    /* L2 共转视角：松手时把拖动期累计的自转角折回 base，从当前朝向无跳变续转 */
    if (pointers.size === 0 && l2SpinFrozen !== null) {
      l2SpinBase = l2SpinFrozen - X3.OMEGA_EM * tSincePeri;
      l2SpinFrozen = null;
    }
  };
  window.addEventListener('pointerup', release);
  window.addEventListener('pointercancel', release);
  /* 滚轮缩放：仅在画布区域内响应 */
  window.addEventListener('wheel', (ev) => {
    if (!overCanvas(ev.clientX, ev.clientY)) return;
    ev.preventDefault();
    if (fpv) { setFov(camera.fov * Math.exp(ev.deltaY * 0.0012)); return; }  // 第一视角：滚轮变焦
    camCtl.dist = clampDist(camCtl.dist * Math.exp(ev.deltaY * 0.0012));
  }, { passive: false });
  window.addEventListener('contextmenu', (ev) => { if (overCanvas(ev.clientX, ev.clientY)) ev.preventDefault(); });

  /* 诊断：?diag=1 时把相机控制事件触发情况写进顶栏（供 dump-dom 读取） */
  if (new URLSearchParams(location.search).has('diag')) {
    let pd = 0, pm = 0, wh = 0;
    window.addEventListener('pointerdown', () => { pd++; __diag(); });
    window.addEventListener('pointermove', () => { pm++; __diag(); });
    window.addEventListener('wheel', () => { wh++; __diag(); }, { passive: true });
    window.__diag = () => {
      const n = $('status-note');
      if (n) n.textContent = ` [diag] pd=${pd} pm=${pm} wh=${wh} ` +
        `az=${camCtl.az.toFixed(3)} el=${camCtl.el.toFixed(3)} dist=${camCtl.dist.toFixed(2)} ` +
        `three=${(window.THREE && THREE.REVISION) || '?'}`;
    };
    window.__diag();
  }
}
function clampDist(d) { return Math.min(800, Math.max(7.5, d)); }
/* 缩放自适应：标注对象世界尺寸 ∝ 相机视距（基准 60），透视投影下屏占比恒定——
   与轨道尺度无关，环月小轨道与 DRO/NRHO 大轨道上文字/点/航天器视觉大小一致。 */
function zoomK() { return camCtl.dist / 60; }
/* L2 会合系视角状态与轨迹包围球（fitView 取景用） */
let l2View = false;          // L2 会合系（旋转坐标系）视角，仅 CR3BP 可用
let topView = false;         // 俯视·地球居中视角（全部月球轨道）：地心惯性俯视，地球锁画面中心
let geoLine = null;          // 地心惯性轨迹线（CR3BP，绿；几何 = 月心显示系位置 − 地球位置）
let moonRing = null;         // 月球地心轨道环（浅蓝；几何 = −earthPosMoonTab，即月心地心路径）
/* L2 视角相机随会合系共转：spin = l2SpinBase + ω_EM·tSincePeri（与会合系环同一
   转角，故环/地月连线在画面中静止）。拖动期间冻结为 l2SpinFrozen 避免与用户
   抢视角；松手时把拖动期累计的自转角折回 l2SpinBase，无跳变地恢复共转。
   所有旋转绕同一固定轴 FRAME_AXIS，标量角相加即精确合成。 */
let l2SpinBase = 0;
let l2SpinFrozen = null;     // 非 null 表示用户正在拖动
let orbitLineSyn = null;     // 会合系静态轨道环（θ=0 锚定变换，地月连线固定）
let inertialInfo = null;     // 惯性系轨迹包围球 { center, radius }
let synInfo = null;          // 会合系轨迹包围球 { center, radius }
function bboxInfo(pts) {
  if (!pts || !pts.length) return null;
  const min = new THREE.Vector3(Infinity, Infinity, Infinity);
  const max = new THREE.Vector3(-Infinity, -Infinity, -Infinity);
  for (const p of pts) { min.min(p); max.max(p); }
  const center = min.clone().add(max).multiplyScalar(0.5);
  let radius = 0;
  for (const p of pts) radius = Math.max(radius, p.distanceTo(center));
  return { center, radius: Math.max(radius, 1e-6) };
}
function fitView() {
  if (fpv) {           // 第一视角下：复位视线偏移与视场
    fpvCam.az = 0; fpvCam.el = 0;
    camera.fov = 50; camera.updateProjectionMatrix();
    return;
  }
  /* 俯视·地球居中（全部月球轨道）：取景中心 = 实时地球位置（applyCamera 每帧
     把 target 压到 earthMini），半径覆盖地月距 + 航天器月心活动包络（CR3BP 取
     惯性轨迹包围球，二体/摄动取远月点距离）——地心轨迹环半径可达 384+90 单位 */
  if (topView && body === 'moon') {
    camCtl.target.copy(earthMini.position);
    const env = (prop && prop.isCr3bp && inertialInfo)
      ? inertialInfo.center.length() + inertialInfo.radius
      : els.a * (1 + els.e) * KM * 1.2;
    camCtl.dist = clampDist((earthMini.position.length() + env) /
      (0.87 * Math.tan(camera.fov * 0.5 * Math.PI / 180)));
    return;
  }
  /* CR3BP：相机对准轨道包围球中心，距离取轨道直径 ≈ 屏高 50%
     （r / (d·tan(fov/2)) = 0.5 → d = 2r/tan(25°)；Halo/NRHO/DRO 通用） */
  if (prop && prop.isCr3bp) {
    if (l2View && synInfo) {
      /* L2 视角：相机锁定月心（原点），使月球居中便于定位。会合系环实时绕
         月心旋转（环心偏向 L2），故取景半径须覆盖 |环心| + 环半径 的整包络。 */
      camCtl.target.set(0, 0, 0);
      const need = synInfo.center.length() + synInfo.radius;
      camCtl.dist = clampDist(need / (0.87 * Math.tan(camera.fov * 0.5 * Math.PI / 180)));
      return;
    }
    if (inertialInfo) {
      camCtl.target.copy(inertialInfo.center);
      camCtl.dist = clampDist(inertialInfo.radius / (0.5 * Math.tan(camera.fov * 0.5 * Math.PI / 180)));
      return;
    }
  }
  camCtl.target.set(0, 0, 0);
  const apoUnits = els.a * (1 + els.e) * KM;
  const baseR = body === 'moon' ? MOON_R : EARTH_R;
  camCtl.dist = clampDist(Math.max(14, apoUnits * 1.7 + baseR * 1.2));
}

/* ── 角度标注弧线 ───────────────────────────────────────────────────
 * ω（近地点幅角）：轨道面内 升交点方向 → 近地点方向，紫色
 * ν（真近点角）  ：轨道面内 近地点方向   → 航天器方向，黄色
 * 两弧都以地心（焦点）为圆心，沿角动量 +h 方向（即运动方向）扫掠，
 * 角度一律归一到 [0, 2π)，保证画出来的是「实际走过」的那段弧。 */
let wArc = null, nuArc = null;
let angleInfo = null;   // { h, periDir, v, rV }：ν 弧每帧更新所需
const ARC_SEG = 80;

function wrap01(x) { x %= C.TWO_PI; if (x < 0) x += C.TWO_PI; return x; }
/* 从 from 绕 normal 到 to 的有向角，归一到 [0,2π) */
function angleAbout(from, to, normal) {
  const c = THREE.MathUtils.clamp(from.dot(to), -1, 1);
  const s = new THREE.Vector3().crossVectors(from, to).dot(normal);
  return wrap01(Math.atan2(s, c));
}
/* 轨道面内半径 r、从 from 沿 +normal 扫 θ 的弧线点列 */
function arcPositions(r, from, vPerp, theta) {
  const pts = [];
  for (let k = 0; k <= ARC_SEG; k++) {
    const t = theta * k / ARC_SEG, ct = Math.cos(t), st = Math.sin(t);
    pts.push(new THREE.Vector3(
      r * (from.x * ct + vPerp.x * st),
      r * (from.y * ct + vPerp.y * st),
      r * (from.z * ct + vPerp.z * st)));
  }
  return pts;
}
/* 弧线中点方向（放角度文字标签用） */
function arcMidDir(from, vPerp, theta) {
  const t = theta / 2;
  return from.clone().multiplyScalar(Math.cos(t)).addScaledVector(vPerp, Math.sin(t)).normalize();
}
/* 角度标签：沿弧线中点径向放在弧外侧 */
function placeArcLabel(sp, midDir, r) {
  const s = labelScale();
  sp.scale.set(s * sp.userData.aspect, s, 1);
  sp.position.copy(midDir).multiplyScalar(r + s * 0.9);
}

/* ── 重建（根数变化时） ──────────────────────────────────────────── */
let trailObj = null;       // 摄动拖尾（含 vertexColors 渐隐）
let trailPts = null;       // {x,y,z} ECI km 数组（ecilToScene 前）
let trailN = 0;            // 拖尾总点数（setDrawRange 上限）
let oscSeries = null;      // 90 d 密切根数等间隔序列（图表用）

/* 演化色（青→橙→红）共享：拖尾顶点色、进度条色条、图表指针色 三处同源 */
const EVOL_C = { a: new THREE.Color(0x4ec9ff), b: new THREE.Color(0xffa552), c: new THREE.Color(0xff6b6b), t: new THREE.Color() };
function evolColorHex(frac) {
  const f = Math.min(1, Math.max(0, frac));
  if (f < 0.5) EVOL_C.t.copy(EVOL_C.a).lerp(EVOL_C.b, f * 2);
  else EVOL_C.t.copy(EVOL_C.b).lerp(EVOL_C.c, (f - 0.5) * 2);
  return `rgb(${(EVOL_C.t.r * 255) | 0},${(EVOL_C.t.g * 255) | 0},${(EVOL_C.t.b * 255) | 0})`;
}
let lastChartTNow = -1;    // 上次重绘图表时的 tNow（避免暂停时重复重绘）

function sceneFromEciArr(x, y, z) {
  return new THREE.Vector3(x * KM, z * KM, -y * KM);
}

/* 拖尾滑动窗口：只显示 tNow 之前一段的轨迹段（drawRange[start, count]），
   不画完整 90 d 拖尾（LEO 千余圈糊成密带看不清），也不从历元累积。
   地球：固定 7 d（LEO/中轨多圈拖尾，用户既有习惯）。
   月球：1 个轨道周期（2026-10-14 ELFO 反馈：月球摄动轨道短周期（ELFO ≈2 h），
   7 d = 84 圈糊带 + 渐进填充期间轨道大半未画出、近/远月点标记悬空在
   无轨迹支撑的一侧，观感「标签跑到月球质心」。改为「最近一圈」——
   航天器渐进画完第一圈后整环恒保留（二体主导，圈圈近似重合），随 J2/地球
   三体摄动缓慢变形，与 CR3BP 三阶段尾迹同一「画完保留」语义）。 */
const TRAIL_WINDOW_S = 7 * 86400;   // 地球拖尾窗口长度：7 天
function trailWindowS() {           // 当前模式拖尾窗口（秒）
  if (body === 'moon' && prop && prop.isPerturbed) {
    const T = prop.periodSeconds();
    if (T > 0) { return T; }
  }
  return TRAIL_WINDOW_S;
}
function trailRevealFrac() {       // 进度条用：tNow 在整条 90 d 窗口中的进度
  if (!prop || !prop.isPerturbed) { return 0; }
  const win = prop.windowSeconds();
  return win > 0 ? Math.min(1, Math.max(0, prop.tNow / win)) : 0;
}
function updateTrailReveal() {
  if (!trailObj || !trailN) { return; }
  const win = prop.windowSeconds();
  if (win <= 0) { trailObj.geometry.setDrawRange(0, 0); return; }
  const end = Math.floor((prop.tNow / win) * trailN);
  const start = Math.floor((Math.max(0, prop.tNow - trailWindowS()) / win) * trailN);
  let count = end - start;
  if (count < 1) count = 1;        /* ≥1 顶点，避免空 range；tNow=0 时只占首点不可见 */
  trailObj.geometry.setDrawRange(start, count);
}

/* CR3BP 轨道渐进尾迹：惯性线（orbitLine，DRO 完整图案 [0,nT]、其余 [0,T] 单份弧）
   画完后恒保留；会合系线（orbitLineSyn，周期 T 真闭合、双圈拼接）每帧按 tSincePeri
   四阶段揭示/消退（画完一圈 → 完整保留一个周期 → 消退 → 半圈尾迹，见函数内注释）；
   反向拖动进度条/ν 相位滑条时同步收缩。cnt=1 时单顶点不画线，起点处不可见。 */
function updateCr3bpReveal() {
  if (!orbitLine || !(prop && prop.isCr3bp)) { return; }
  const T = prop.periodSeconds();
  if (!(T > 0)) { return; }
  const tC = tSincePeri;
  /* 俯视勾选闸门（每帧权威）：航天器蓝线（orbitLine，月心惯性轨迹短弧）与
     绿线（geoLine，地心惯性轨迹）的可见性在俯视 CR3BP 下由 chk-trail-moon /
     chk-trail-earth 决定。必须在每帧揭示前强制同步——否则拖动进度条 / 改根数 /
     切预设会经 setL2View·setTopView·rebuildCommonAfterProp 的
     「orbitLine.visible = !l2View」把勾选的 false 覆盖回 true，导致用户取消
     勾选后短弧仍残留（2026-10-16 用户反馈）。非俯视或非 CR3BP 不动 visible，
     由 setL2View/setTopView/applyDisplay 管理。moonRing 恒显，不在此列。 */
  if (topView) {
    const cm = $('chk-trail-moon'), ce = $('chk-trail-earth');
    if (cm) {
      orbitLine.visible = cm.checked;
      /* trailObj（vertexColors 渐隐拖尾，月心显示系坐标）同样是「航天器相对
         月球的轨迹」，必须与蓝线同一勾选管辖——此前漏管导致取消勾选后
         蓝线消失但渐隐拖尾仍残留（2026-10-16 用户二次反馈） */
      if (trailObj) trailObj.visible = cm.checked;
    }
    if (ce && geoLine) geoLine.visible = ce.checked;
  }
  /* 惯性线（DRO 为完整惯性图案 [0,nT]，其余为单份 [0,T] 第一圈弧）：渐进画出后
     恒保留——Halo 惯性轨迹缓慢进动不闭合，画完即止（2026-10-14 DRO 横穿弦修复，
     语义 = 「保留第一圈轨迹」）；DRO n:1 共振惯性图案 nT 才闭合，2026-10-15 起
     采样跨度即为 nT（inertialSpanSeconds），揭示同步按 nT——修复「固定残存半个
     周期」（原 [0,T] 弧只是惯性图案的一半，航天器后半程飞出已画轨迹） */
  const TI = prop.inertialSpanSeconds ? prop.inertialSpanSeconds() : T;
  const nI = orbitLine.geometry.attributes.position.count - 1;
  const i1I = Math.round(Math.min(Math.max(tC, 0), TI) / TI * nI);
  orbitLine.geometry.setDrawRange(0, i1I + 1);
  /* 地心轨迹（俯视·地球居中）：跨度 = 恒星月（userData.span，几何构建时写入），
     与蓝线节奏解耦——halo/nrho 蓝线按 T 画完即止，绿线继续绕地球画满一整圈
     （2026-10-15 修复「绿线画半圈即停」：原与蓝线共用 inertialSpanSeconds） */
  if (geoLine) {
    const TG = geoLine.userData.span || C.TWO_PI / X3.OMEGA_EM;
    const nG = geoLine.geometry.attributes.position.count - 1;
    geoLine.geometry.setDrawRange(0, Math.round(Math.min(Math.max(tC, 0), TG) / TG * nG) + 1);
  }
  /* 会合系线（周期 T 真闭合，双圈拼接 2n+1 顶点，顶点 i 时刻 = i/n·T）四阶段尾迹
     （2026-10-14 用户反馈：画完一圈立即消退太早，应完整保留一个周期再消失）：
     ① t ≤ T          渐进画第一圈 [0, t]，无消退
     ② T < t ≤ 2T     整圈完整保留（画完一圈后再留一个周期，消失从 2T 才开始）
     ③ 2T < t ≤ 3.5T  起点以 2× 消退（sAbs = 2(t−2T)），窗口 2T → 0.5T
     ④ t > 3.5T       半圈尾迹窗口 [t−T/2, t] 恒定跟随
     各阶段 sAbs 连续（T 处 0，2T 处 0，3.5T 处 3T≡0.5T 与稳态式一致）。
     窗口表示：窗口长 > T 或起点相位在后（s > e）→ 跨圈段 [sIdx, n+eIdx]
     （第二圈顶点与第一圈重合，外观不变）；否则同相直段 [sIdx, eIdx]。
     注意跨圈判据必须看窗口长：t∈(2T,3T) 时 s<e 但窗口超一圈，
     直段公式会漏掉第一圈尾段。 */
  if (orbitLineSyn) {
    const n = (orbitLineSyn.geometry.attributes.position.count - 1) / 2;
    const eIdx = Math.round((tC % T) / T * n);
    let i0, i1;
    if (tC <= T) {
      i0 = 0;
      i1 = Math.round(Math.min(Math.max(tC, 0), T) / T * n);
    } else if (tC <= 2 * T) {
      i0 = 0;
      i1 = n + eIdx;
    } else {
      const sAbs = tC <= 3.5 * T ? (tC - 2 * T) * 2 : tC - T / 2;
      const s = sAbs % T, e = tC % T;
      if (tC - sAbs >= T) {
        /* 窗口 ≥ 一圈：轨迹经过每个相位，视觉恒整环——画圈1整圈 + 圈2头。
           不能用 [sIdx, n+eIdx] 首尾段表示：窗口中间的整圈会丢失
           （t∈(2T,3T) 起点推进后首段缩短，环上并集 < 整环 → 视觉跳变） */
        i0 = 0;
        i1 = n + eIdx;
      } else {
        i0 = Math.round(s / T * n);
        i1 = (s > e ? n : 0) + eIdx;
      }
    }
    orbitLineSyn.geometry.setDrawRange(i0, i1 - i0 + 1);
  }
}
/* 进度条 #ph-prog：--prog = tNow/window，色条与拖尾同源（青→橙→红） */
function updateProgBar() {
  const el = $('ph-prog');
  if (el) el.style.setProperty('--prog', trailRevealFrac().toFixed(4));
}

function rebuildTrail() {
  if (trailObj) {
    scene.remove(trailObj);
    trailObj.geometry.dispose();
    trailObj.material.dispose();
    trailObj = null;
  }
  trailPts = null;
  trailN = 0;
  if (!(prop && prop.isPerturbed)) { return; }
  /* 50000 点：LEO（T≈92 min）下 stride≈4 步 → 每圈 ~30 点，避免稀疏采样
     导致 LINE_STRIP 相邻点跨越半个轨道、产生横穿椭圆的弦状伪线 */
  const raw = prop.trailPoints(50000);
  const n = raw.x.length;
  const posArr = new Float32Array(n * 3);
  const colArr = new Float32Array(n * 3);
  const cStart = new THREE.Color(0x4ec9ff);   // 青：起点
  const cMid = new THREE.Color(0xffa552);     // 橙：中段
  const cEnd = new THREE.Color(0xff6b6b);     // 红：终点
  const tmp = new THREE.Color();
  for (let i = 0; i < n; i++) {
    const v = sceneFromEciArr(raw.x[i], raw.y[i], raw.z[i]);
    posArr[i * 3] = v.x; posArr[i * 3 + 1] = v.y; posArr[i * 3 + 2] = v.z;
    const t = i / (n - 1);
    if (t < 0.5) { tmp.copy(cStart).lerp(cMid, t * 2); }
    else { tmp.copy(cMid).lerp(cEnd, (t - 0.5) * 2); }
    colArr[i * 3] = tmp.r; colArr[i * 3 + 1] = tmp.g; colArr[i * 3 + 2] = tmp.b;
  }
  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.BufferAttribute(posArr, 3));
  geo.setAttribute('color', new THREE.BufferAttribute(colArr, 3));
  trailObj = new THREE.FatLine(geo, {
    vertexColors: true, transparent: true, opacity: 0.55,
  });
  scene.add(trailObj);
  trailPts = raw;
  trailN = n;
  updateTrailReveal();   /* 重建后 tNow=0 → 窗口空；播放/拖动时只显 tNow 前 7 d 尾迹 */
}

/* 90 d 密切根数序列：每 prop.h 步取一个 osc（降采样到 ≤ 600 点） */
function buildOscSeries() {
  oscSeries = null;
  if (!(prop && prop.isPerturbed)) { return; }
  const n = prop.nSteps;
  const stride = Math.max(1, Math.ceil(n / 600));
  const s = { t: [], a: [], e: [], i: [], Om: [], w: [] };
  for (let k = 0; k <= n; k += stride) {
    const y = prop.yArr.slice(k * 6, k * 6 + 6);
    const o = window.ElementsPerturb.cartToClassical(y.slice(0, 3), y.slice(3, 6));
    s.t.push(k * prop.h / 86400);          // 天
    s.a.push(o.a);
    s.e.push(o.e);
    s.i.push(o.i * C.R2D);
    s.Om.push(o.Om * C.R2D);
    s.w.push(o.w * C.R2D);
  }
  oscSeries = s;
}

function drawChart(canvas, data, color, label) {
  const ctx = canvas.getContext('2d');
  const W = canvas.width = canvas.clientWidth * (window.devicePixelRatio || 1);
  const H = canvas.height = canvas.clientHeight * (window.devicePixelRatio || 1);
  ctx.clearRect(0, 0, W, H);
  if (!data || data.length < 2) return;
  const pad = 4 * (window.devicePixelRatio || 1);
  let mn = Infinity, mx = -Infinity;
  for (const v of data) { if (v < mn) mn = v; if (v > mx) mx = v; }
  const span = mx - mn || 1;
  const lo = mn - span * 0.08, hi = mx + span * 0.08;
  ctx.beginPath();
  ctx.moveTo(pad, H - pad);
  for (let i = 0; i < data.length; i++) {
    const x = pad + (W - 2 * pad) * i / (data.length - 1);
    const y = H - pad - (H - 2 * pad) * (data[i] - lo) / (hi - lo);
    ctx.lineTo(x, y);
  }
  ctx.strokeStyle = color;
  ctx.lineWidth = 3 * (window.devicePixelRatio || 1);
  ctx.stroke();
  // 当前时刻指针：竖虚线 + 顶部三角，色随演化进度（青→橙→红，与进度条同源）
  if (prop && prop.isPerturbed && oscSeries) {
    const frac = Math.min(1, Math.max(0, prop.tNow / prop.windowSeconds()));
    const fx = pad + (W - 2 * pad) * frac;
    const dpr = (window.devicePixelRatio || 1);
    ctx.beginPath();
    ctx.moveTo(fx, 0); ctx.lineTo(fx, H);
    ctx.strokeStyle = 'rgba(255,255,255,0.4)';
    ctx.lineWidth = 1 * dpr;
    ctx.setLineDash([3 * dpr, 3 * dpr]);
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = evolColorHex(frac);
    const ts = 5 * dpr;
    ctx.beginPath();
    ctx.moveTo(fx - ts, 0); ctx.lineTo(fx + ts, 0); ctx.lineTo(fx, ts * 1.5);
    ctx.closePath(); ctx.fill();
  }
  // 终点值
  ctx.font = `${10 * (window.devicePixelRatio || 1)}px ui-monospace, Menlo, monospace`;
  ctx.fillStyle = color;
  ctx.textAlign = 'right';
  ctx.fillText(data[data.length - 1].toFixed(data[data.length - 1] < 10 ? 5 : 2), W - 4, 12 * (window.devicePixelRatio || 1));
}

function drawCharts() {
  if (!oscSeries) { return; }
  drawChart($('cv-a'), oscSeries.a, '#4ec9ff', 'a');
  drawChart($('cv-e'), oscSeries.e, '#ffd166', 'e');
  drawChart($('cv-i'), oscSeries.i, '#7dffa8', 'i');
  drawChart($('cv-om'), oscSeries.Om, '#b07cff', 'Ω');
}

function updateOscTable(st) {
  if (!(prop && prop.isPerturbed)) { return; }
  const ep = window.ElementsPerturb;
  const now = st ? st.osc : prop.osculatingElements();
  $('osc-a0').textContent = fmt(els.a, 3);
  $('osc-e0').textContent = fmt(els.e, 5);
  $('osc-i0').textContent = fmtDeg(els.i, 3);
  $('osc-om0').textContent = fmtDeg(els.Om, 3);
  $('osc-w0').textContent = fmtDeg(els.w, 3);
  $('osc-anow').textContent = fmt(now.a, 3);
  $('osc-enow').textContent = fmt(now.e, 5);
  $('osc-inow').textContent = fmtDeg(now.i * C.R2D, 4);
  $('osc-omnow').textContent = fmtDeg(now.Om * C.R2D, 3);
  $('osc-wnow').textContent = fmtDeg(now.w * C.R2D, 3);
  const tend = prop.stateAtTime(prop.windowSeconds()).osc;
  $('osc-aend').textContent = fmt(tend.a, 3);
  $('osc-eend').textContent = fmt(tend.e, 5);
  $('osc-iend').textContent = fmtDeg(tend.i * C.R2D, 4);
  $('osc-omend').textContent = fmtDeg(tend.Om * C.R2D, 3);
  $('osc-wend').textContent = fmtDeg(tend.w * C.R2D, 3);
  const f = prop.forceMagnitudes();
  const m2s = (v) => v < 1e-12 ? '<1e-12' : v.toExponential(2);
  const bits = [];
  if ($('pp-grav').checked) bits.push(`<span class="fz">zonal ${m2s(f.zonal * 1000)}</span>`);
  if ($('pp-drag').checked) bits.push(`<span class="fd">drag ${m2s(f.drag * 1000)}</span>`);
  if ($('pp-srp').checked) bits.push(`<span class="fs">srp ${m2s(f.srp * 1000)}</span>`);
  $('force-line').innerHTML = `摄动加速度 (m/s²)<br>${bits.join(' · ') || '（全部关闭）'}`;
}

/* ── 重建（根数变化时） ──────────────────────────────────────────── */
function rebuild() {
  const sel = $('sel-prop');
  const cls = C.propagators[sel ? sel.value : 'kepler'] || C.KeplerPropagator;
  const nu0 = els.nu * C.D2R;
  if (cls === X3.Cr3bpPropagator) {
    /* CR3BP：会合系三体传播器（月球选项卡专属） */
    $('pert-busy-moon').textContent = '三体传播计算中…';
    prop = new X3.Cr3bpPropagator({
      /* 周期轨道族用预解初值；其余（ELFO 等根数预设）由当前根数构建初值 */
      preset: cr3bpState.active ? cr3bpState.preset : 'elements',
      elements: { ...elsRad(), nu: nu0 },
      azKm: cr3bpState.azKm,
      phase: cr3bpState.phase,
      stationKeeping: cr3bpState.stationKeeping && cr3bpState.preset === 'halo',
    });
    $('pert-busy-moon').textContent = '';
    rebuildCommonAfterProp();
    return;
  }
  if (cls === X3.MoonCowellPropagator) {
    $('pert-busy-moon').textContent = '传播计算中…';
    prop = new X3.MoonCowellPropagator(elsRad(), {
      grav: $('mp-grav').checked,
      earth3B: $('mp-earth3b').checked,
      srp: $('mp-srp').checked,
      CR: parseFloat($('pp-cr').value) || 1.5,
      SM: parseFloat($('pp-sm').value) || 0.02,
    });
    $('pert-busy-moon').textContent = '';
  } else if (cls.hasPerturbation) {
    /* 类级摄动标记用契约里的 static hasPerturbation（isPerturbed 仅为实例属性） */
    $('pert-busy').textContent = '传播计算中…';
    const opts = {
      grav: $('pp-grav').checked,
      drag: $('pp-drag').checked,
      srp: $('pp-srp').checked,
      moon3B: $('pp-moon3b').checked,
      Cd: parseFloat($('pp-cd').value) || 2.2,
      CR: parseFloat($('pp-cr').value) || 1.5,
      SM: parseFloat($('pp-sm').value) || 0.02,
      lam0: (parseFloat($('pp-lam').value) || 0) * C.D2R,
    };
    prop = new cls(elsRad(), opts);
    $('pert-busy').textContent = '';
  } else {
    /* 开普勒：按中心天体选 μ/R */
    prop = (body === 'moon')
      ? new C.KeplerPropagator(elsRad(), X3.MU_MOON, X3.R_MOON)
      : new cls(elsRad());
  }
  rebuildCommonAfterProp();
}

/* rebuild 的公共后段（轨道线/标记/网格/弧线/派生量） */
function rebuildCommonAfterProp() {
  const nu0 = els.nu * C.D2R;
  const isCr3bp = !!(prop && prop.isCr3bp);

  /* 轨道椭圆（CR3BP：整周期采样实际轨迹，非开普勒椭圆） */
  const pts = [];
  let cr3bpPts = null;      /* CR3BP 惯性轨迹场景点：拱点极值扫描用（见下） */
  if (isCr3bp) {
    /* DRO 惯性图案跨度 nT（2:1→2T、3:1→3T），采样点数同比放大保持弧顺滑 */
    const spanK = prop.inertialSpanSeconds ? prop.inertialSpanSeconds() / prop.periodSeconds() : 1;
    const path = prop.orbitPath(Math.round(300 * spanK));
    for (const rE of path) pts.push(eciToScene(rE));
    cr3bpPts = pts;
    /* 惯性线保持单份完整图案（2026-10-14 修复 DRO「中间辅助线」；
       2026-10-15 修复 DRO「固定残存半个周期」：采样跨度 T → inertialSpanSeconds
       = nT，惯性系 nT 闭合，见 elements_cr3bp.js 该函数注释）：
       Halo/NRHO 惯性轨迹不以 T 闭合（缓慢进动），重复拼接会在 t=T 处跳回
       t≈0 位置，产生横穿轨道的弦（实测 144 场景单位），故恒为单份 [0,T] 弧。
       渐进揭示语义见 updateCr3bpReveal：图案画完后恒保留。 */
  } else {
    const oscOrb = (prop && prop.isPerturbed) ? prop.osculatingElements() : null;
    for (let k = 0; k <= 256; k++) {
      const nu = (k / 256) * C.TWO_PI;
      const rE = oscOrb ? C.classicalToState(oscOrb, nu).rE : prop.stateAtTrueAnomaly(nu).rE;
      pts.push(eciToScene(rE));
    }
  }
  if (orbitLine) { scene.remove(orbitLine); orbitLine.geometry.dispose(); }
  orbitLine = new THREE.FatLine(
    new THREE.BufferGeometry().setFromPoints(pts),
    { color: 0x4ec9ff, transparent: true, opacity: 0.85 }
  );
  scene.add(orbitLine);

  /* CR3BP：额外构建会合系（旋转坐标系）静态轨道环，L2 视角下显示；
     惯性/会合系轨迹包围球供 fitView 取景（默认视图放大到轨道充满 ~50% 屏高） */
  if (orbitLineSyn) { scene.remove(orbitLineSyn); orbitLineSyn.geometry.dispose(); orbitLineSyn = null; }
  inertialInfo = bboxInfo(pts);
  synInfo = null;
  if (isCr3bp) {
    const synPath = prop.orbitPath(300, true);
    const synPts = synPath.map(eciToScene);
    for (let k = 1; k < synPath.length; k++) synPts.push(eciToScene(synPath[k]));
    orbitLineSyn = new THREE.FatLine(
      new THREE.BufferGeometry().setFromPoints(synPts),
      { color: 0xffa552, transparent: true, opacity: 0.9 }
    );
    orbitLineSyn.visible = l2View;
    scene.add(orbitLineSyn);
    synInfo = bboxInfo(synPts);
    orbitLine.visible = !l2View;
    /* 不预渲染整周期：初始 drawRange=0 隐藏，由航天器随传播进度逐段画出
       （updateCr3bpReveal 每帧按 tSincePeri 揭示，见拖尾函数旁） */
    orbitLine.geometry.setDrawRange(0, 0);
    orbitLineSyn.geometry.setDrawRange(0, 0);
  }

  /* 俯视·地球居中视角资源：地心惯性轨迹（绿，仅 CR3BP）+ 月球地心轨道环（浅蓝）。
     几何一律存「地心坐标」：航天器地心位置 = 月心显示系位置 − 地球位置；
     月球地心路径 = −earthPosMoonTab。线对象 position 每帧压到实时地球位置
     （updateScene），屏幕上即环绕地球的路径。采样跨度 = 恒星月 tSid
     （地心视角绕地球一圈的自然周期，全部 CR3BP 预设统一）：DRO n:1 的 nT
     恰为 1 个恒星月——地心轨迹绕地球整圈准闭合；halo/nrho 原沿用
     inertialSpanSeconds（=T，仅 54%/24% 恒星月）导致绿线画到一半就停
     （2026-10-15 用户反馈，实测 halo 扫掠 187.9°）。halo 的会合系相位
     tSid/T ≈ 1.84 非整数倍 → 一圈末端与起点不重合，属物理真实而非 bug。 */
  if (geoLine) { scene.remove(geoLine); geoLine.geometry.dispose(); geoLine = null; }
  if (moonRing) { scene.remove(moonRing); moonRing.geometry.dispose(); moonRing = null; }
  if (body === 'moon') {
    const nR = 128, tSid = C.TWO_PI / X3.OMEGA_EM;
    const rPts = [];
    for (let k = 0; k <= nR; k++) rPts.push(eciToScene(X3.earthPosMoonTab(k / nR * tSid)).negate());
    moonRing = new THREE.FatLine(
      new THREE.BufferGeometry().setFromPoints(rPts),
      { color: 0x6fa8ff, transparent: true, opacity: 0.35 }
    );
    moonRing.visible = false;   /* 可见性由 applyDisplay 按 topView 集中管理 */
    scene.add(moonRing);
  }
  if (isCr3bp) {
    const spanGeo = C.TWO_PI / X3.OMEGA_EM;   /* 恒星月：地心一圈自然跨度 */
    const geoPath = prop.orbitPath(360, false, spanGeo);
    const nG = geoPath.length - 1;   /* 非 DRO 分支返回 max(240, 2·nPts)+1 点——
        必须从实际长度推导，地球位置时间轴 k/nG·spanGeo 才与采样点时刻对齐 */
    const gPts = geoPath.map((rE, k) =>
      eciToScene(rE).sub(eciToScene(X3.earthPosMoonTab(k / nG * spanGeo))));
    geoLine = new THREE.FatLine(
      new THREE.BufferGeometry().setFromPoints(gPts),
      { color: 0x7dffa8, transparent: true, opacity: 0.9 }
    );
    geoLine.userData.span = spanGeo;   /* updateCr3bpReveal 按 tSid 揭示 */
    geoLine.geometry.setDrawRange(0, 0);   /* 由 updateCr3bpReveal 按 tSincePeri 揭示 */
    geoLine.visible = false;
    scene.add(geoLine);
  }

  /* 近地点 / 远地点标记（摄动模式下用瞬时密切根数；CR3BP 周期轨道族隐藏；
     CR3BP 根数初值模式（ELFO）显示——其 ic 由开普勒根数构建、二体主导，
     与二体同待遇给出拱点。2026-10-15 用户反馈「远月点仍在月心」） */
  const oscNow = (prop && prop.isPerturbed && !isCr3bp) ? prop.osculatingElements() : null;
  const mkEl = oscNow || elsRad();
  const mkElsDeg = oscNow
    ? { a: oscNow.a, e: oscNow.e, i: oscNow.i * C.R2D, Om: oscNow.Om * C.R2D, w: oscNow.w * C.R2D }
    : els;
  const bodyR = body === 'moon' ? X3.R_MOON : C.R_EARTH;
  const periSt = isCr3bp ? null : C.classicalToState(mkEl, 0, body === 'moon' ? X3.MU_MOON : undefined, body === 'moon' ? X3.R_MOON : undefined);
  const periE = isCr3bp ? null : periSt.rE;
  const apoE = isCr3bp ? null : C.classicalToState(mkEl, Math.PI, body === 'moon' ? X3.MU_MOON : undefined, body === 'moon' ? X3.R_MOON : undefined).rE;
  const hp = mkElsDeg.a * (1 - mkElsDeg.e) - bodyR;
  /* 拱点标记是否显示：非三体恒显示；三体中仅根数初值模式（preset 'elements'，
     即 ELFO）显示。周期轨道族（halo/nrho/dro*）的拱点在会合系中不落在惯性
     轨迹上，故隐藏。apsisAvailable 亦被 applyDisplay 读取——否则标签开关会把
     隐藏的拱点标签强制点亮在默认原点（= 月心），即用户所见 bug。 */
  apsisAvailable = !isCr3bp || (prop && prop.preset === 'elements');
  /* 拱点位置（场景坐标）：根数初值模式取 CR3BP 实际轨迹的月心距极值，
     保证标记落在画出的轨迹上（实测远月点较开普勒椭圆高 ~120 km） */
  let periScene = null, apoScene = null;
  if (apsisAvailable) {
    if (isCr3bp) {
      let rMin = Infinity, rMax = -Infinity;
      for (const pt of (cr3bpPts || [])) {
        const r = pt.length();
        if (r < rMin) { rMin = r; periScene = pt.clone(); }
        if (r > rMax) { rMax = r; apoScene = pt.clone(); }
      }
    } else {
      periScene = eciToScene(periE);
      apoScene = eciToScene(apoE);
    }
  }
  if (!periMark) {
    periMark = new THREE.Mesh(new THREE.SphereGeometry(0.38, 16, 12),
      new THREE.MeshBasicMaterial({ color: 0x7dffa8 }));
    apoMark = new THREE.Mesh(new THREE.SphereGeometry(0.38, 16, 12),
      new THREE.MeshBasicMaterial({ color: 0xffa552 }));
    scene.add(periMark); scene.add(apoMark);
  }
  periMark.material.color.set(hp < 0 ? 0xff6b6b : 0x7dffa8);
  periMark.visible = apsisAvailable;
  apoMark.visible = apsisAvailable;
  if (periScene) { periMark.position.copy(periScene); }
  if (apoScene) { apoMark.position.copy(apoScene); }

  /* 赤道网格（半径随远地点自适应；CR3BP 模式覆盖为地月距离尺度） */
  gridR = isCr3bp
    ? Math.min(420, X3.A_EM * KM * 1.25)
    : Math.min(420, Math.max(10, els.a * (1 + els.e) * KM * 1.15));
  /* 赤道面网格用 THREE.PolarGridHelper（GL_LINES，硬件恒定 1px），与 3px 粗的
     轨迹线明确区分（2026-10-16 用户反馈：FatLine 化后网格与轨迹同为 3px 难以分辨） */
  if (grid) { scene.remove(grid); grid.geometry.dispose(); grid.material.dispose(); }
  grid = new THREE.PolarGridHelper(gridR, 24, 6, 96, 0x2a3b52, 0x1b2838);
  grid.material.transparent = true; grid.material.opacity = 0.85;
  scene.add(grid);

  /* 升交点节线（沿 ECI (cosΩ, sinΩ, 0)；CR3BP 隐藏） */
  const nodeDir = new THREE.Vector3(Math.cos(mkElsDeg.Om * C.D2R), 0, -Math.sin(mkElsDeg.Om * C.D2R));
  if (nodeLine) { scene.remove(nodeLine); nodeLine.geometry.dispose(); }
  nodeLine = new THREE.FatLine(
    new THREE.BufferGeometry().setFromPoints(
      [new THREE.Vector3(), nodeDir.clone().multiplyScalar(gridR)]),
    { color: 0xb07cff, dashed: true, dashSize: 0.8, gapSize: 0.6, transparent: true, opacity: 0.8 }
  );
  nodeLine.computeLineDistances();
  nodeLine.visible = !isCr3bp;
  scene.add(nodeLine);

  /* 标签 */
  for (const k of Object.keys(labels)) { scene.remove(labels[k]); delete labels[k]; }
  labels.peri = makeLabel(body === 'moon' ? '近月点' : '近地点', '#7dffa8');
  labels.apo = makeLabel(body === 'moon' ? '远月点' : '远地点', '#ffa552');
  labels.aries = makeLabel('♈ 春分点 (X)', '#66ddaa');
  labels.north = makeLabel('Z 北极', '#7dffa8');
  labels.node = makeLabel('☊ 升交点', '#b07cff');
  labels.sc = makeLabel('航天器', '#4ec9ff');
  /* 地月距离远超月心附近取景范围，故月球选项卡恒有「地球」标签（每帧锚在真实地球球体上沿） */
  if (body === 'moon') labels.earth = makeLabel('地球', '#6fa8ff');
  if (apsisAvailable) {
    placeLabel(labels.peri, periMark.position);
    placeLabel(labels.apo, apoMark.position);
    /* 锚点登记：静态标签由 updateScene 每帧按 zoomK 重放（sc/earth/ν/L点标签本已每帧跟随） */
    labels.peri.userData.anchor = { pos: periMark.position.clone() };
    labels.apo.userData.anchor = { pos: apoMark.position.clone() };
  } else {
    /* 无拱点（周期轨道族）：标签不可见，且不登记 anchor——
       否则 updateScene 的重放循环会把它们留在原点（月心） */
    labels.peri.visible = labels.apo.visible = false;
    delete labels.peri.userData.anchor;
    delete labels.apo.userData.anchor;
  }
  placeLabel(labels.aries, new THREE.Vector3(9.6, 0, 0), 0.4);
  labels.aries.userData.anchor = { pos: new THREE.Vector3(9.6, 0, 0), dy: 0.4 };
  placeLabel(labels.north, new THREE.Vector3(0, 9.6, 0), 0.4);
  labels.north.userData.anchor = { pos: new THREE.Vector3(0, 9.6, 0), dy: 0.4 };
  const nodeEnd = nodeDir.clone().multiplyScalar(gridR);
  placeLabel(labels.node, nodeEnd, 0.4);
  labels.node.userData.anchor = { pos: nodeEnd.clone(), dy: 0.4 };

  /* 角度标注弧：ω（升交点→近地点，紫）、ν（近地点→航天器，黄）
     摄动模式下 ω 弧用瞬时密切根数；ν 弧由 updateScene 每帧覆写；CR3BP 隐藏 */
  const periSt2 = isCr3bp ? null : periSt;
  const hE = periSt2 ? new THREE.Vector3(periSt2.rE[0], periSt2.rE[1], periSt2.rE[2])
    .cross(new THREE.Vector3(periSt2.vE[0], periSt2.vE[1], periSt2.vE[2])).normalize() : null;
  const hScene = hE ? new THREE.Vector3(hE.x, hE.z, -hE.y) : null;
  const periDir = isCr3bp ? null : periMark.position.clone().normalize();

  /* 弧半径随轨道尺度自适应：永远大于中心天体半径、不超过网格的 60% */
  const bodyRUnits = body === 'moon' ? MOON_R : EARTH_R;
  const periUnits = mkElsDeg.a * (1 - mkElsDeg.e) * KM;
  const rW = Math.max(bodyRUnits * 1.6, Math.min(gridR * 0.6, periUnits * 1.6));
  const rV = Math.max(bodyRUnits * 1.15, rW * 0.72);

  if (wArc) { scene.remove(wArc); wArc.geometry.dispose(); wArc.material.dispose(); }
  if (nuArc) { scene.remove(nuArc); nuArc.geometry.dispose(); nuArc.material.dispose(); }

  if (!isCr3bp && hScene) {
    const vW = new THREE.Vector3().crossVectors(hScene, nodeDir);
    const thetaW = angleAbout(nodeDir, periDir, hScene);
    wArc = new THREE.FatLine(
      new THREE.BufferGeometry().setFromPoints(arcPositions(rW, nodeDir, vW, thetaW)),
      { color: 0xb07cff, transparent: true, opacity: 0.95 }
    );

    const craftDir = eciToScene(prop.stateAtTrueAnomaly(nu0).rE).normalize();
    const vV = new THREE.Vector3().crossVectors(hScene, periDir);
    const thetaNu = angleAbout(periDir, craftDir, hScene);
    nuArc = new THREE.FatLine(
      new THREE.BufferGeometry().setFromPoints(arcPositions(rV, periDir, vV, thetaNu)),
      { color: 0xffd166, transparent: true, opacity: 0.95 }
    );
    scene.add(wArc); scene.add(nuArc);
    angleInfo = { h: hScene, periDir, v: vV, rV };

    labels.omega = makeLabel('ω', '#b07cff');
    labels.nu = makeLabel('ν', '#ffd166');
    placeArcLabel(labels.omega, arcMidDir(nodeDir, vW, thetaW), rW);
    labels.omega.userData.arcAnchor = { midDir: arcMidDir(nodeDir, vW, thetaW), r: rW };
    placeArcLabel(labels.nu, arcMidDir(periDir, vV, thetaNu), rV);
  } else {
    wArc = null; nuArc = null;
    angleInfo = null;
  }

  /* 重建后按当前显示开关刷新可见性（否则新弧会无视复选框状态） */
  applyDisplay();

  /* 派生量 + 警告 */
  const T = prop.periodSeconds();
  $('d-T').textContent = fmtT(T);
  if (isCr3bp) {
    const mp = prop.missionParams();
    $('d-hp').textContent = `${fmt(mp.periluneKm)} km`;
    $('d-ha').textContent = `${fmt(mp.apoluneKm)} km`;
  } else {
    $('d-hp').textContent = `${fmt(hp)} km`;
    $('d-ha').textContent = `${fmt(mkElsDeg.a * (1 + mkElsDeg.e) - bodyR)} km`;
  }
  if (body === 'moon') {
    $('warn-atmo').classList.remove('show');
    $('warn-touch').classList.remove('show');
    $('warn-moon-touch').classList.toggle('show', !isCr3bp && hp < 0);
  } else {
    $('warn-moon-touch').classList.remove('show');
    $('warn-atmo').classList.toggle('show', hp >= 0 && hp < 120);
    $('warn-touch').classList.toggle('show', hp < 0);
  }
  $('hint-ecirc').classList.toggle('show', !isCr3bp && mkElsDeg.e < 0.005);

  /* 相位从当前 ν 重新锚定（保持连续性） */
  const tPrev = tSincePeri;
  tSincePeri = isCr3bp ? 0 : prop.timeAtTrueAnomaly(nu0);
  /* UTC 历元：首次对齐到真实当前 UTC；之后随 tSincePeri 重锚定等量平移，
     保证仿真时钟读数跨重建（改根数/切预设/切模型）保持连续不跳变 */
  if (epochUTCms === null) syncEpochToNow();
  else epochUTCms += (tPrev - tSincePeri) * 1000;

  /* 摄动模式专属：拖尾 + 90 d 序列 + 图表/表格 */
  rebuildTrail();
  buildOscSeries();
  const showEvol = !!(prop && prop.isPerturbed && !isCr3bp);
  $('card-evol').style.display = showEvol ? '' : 'none';
  /* 时间进度条：摄动与 CR3BP 共用（都是 90 d 窗口、isPerturbed 契约），
     CR3BP 也可按时间自由拖动 */
  $('ph-prog').style.display = (showEvol || isCr3bp) ? '' : 'none';
  $('ph-hint').style.display = showEvol ? '' : 'none';
  if (showEvol) { drawCharts(); updateOscTable(null); }
  if (showEvol || isCr3bp) updateProgBar();

  /* CR3BP / 月心模式专属卡片 */
  $('card-cr3bp').style.display = isCr3bp ? '' : 'none';
  $('card-mission').style.display = isCr3bp ? '' : 'none';
  /* L 点标记（CR3BP 时显示） */
  lpointGroup.visible = isCr3bp;
  /* 受力分析勾选行仅 CR3BP 模式显示（箭头可见性由 applyDisplay 门控） */
  $('row-forces').style.display = isCr3bp ? '' : 'none';
  /* L2 会合系按钮仅 CR3BP 可用；离开三体模式自动退出 L2 视角。
     L2 视角下隐藏惯性系拖尾（会合系中无意义） */
  $('btn-l2').disabled = !isCr3bp;
  if (!isCr3bp && l2View) setL2View(false, true);
  /* 俯视·地球静止视角仅月球选项卡可用；切回地球自动退出 */
  $('btn-top').disabled = body !== 'moon';
  if (body !== 'moon' && topView) setTopView(false, true);
  if (trailObj) trailObj.visible = !l2View;
}

/* ── 每帧场景更新 ────────────────────────────────────────────────── */
const Z_AXIS = new THREE.Vector3(0, 0, 1);
function updateScene(st) {
  /* 航天器位置/速度一律取实时显示系读数（rotStateToDisplay 已含 ω·t 旋转）。
     L2 视角下会合系环按同一 θ 刚性旋转（见下方 orbitLineSyn），二者恒等，
     故无需再对 st.rot 做 θ=0 冻结特判。 */
  const p = eciToScene(st.rE);
  scGroup.position.copy(p);

  /* ── 天体位置与自转（按中心天体分派） ── */
  if (body === 'earth') {
    /* 月球迷你天体：圆化解析轨道 + 同步自转（潮汐锁定） */
    const rM = window.ElementsPerturb.moonPosEci(tSincePeri);
    moonMesh.position.copy(eciToScene(rM));
    moonMesh.rotation.y = 0;   /* 潮汐锁定：正面朝地球（贴图 UV 对应面向地球侧） */
    earthMini.visible = false;
  } else {
    /* 月球选项卡：月球在原点自转（潮汐锁定 ω_EM），地球迷你球沿圆化轨道实时公转 */
    moonMesh.position.set(0, 0, 0);
    const rE = X3.earthPosMoonTab(tSincePeri);
    const ePos = eciToScene(rE);
    earthMini.position.copy(ePos);
    earthMini.visible = true;
    /* 俯视·地球居中：地心轨迹线 / 月球轨道环的几何原点在「地心」，
       每帧压到实时地球位置，屏幕上即环绕地球的静态路径（相机随 target 锁地球） */
    if (topView) {
      if (geoLine) geoLine.position.copy(ePos);
      if (moonRing) moonRing.position.copy(ePos);
    }
    /* 潮汐锁定：月球自转角 = 地球方位角（绕显示系 y 轴），随公转实时推进 */
    const th = X3.THETA_E0 + X3.OMEGA_EM * tSincePeri;
    moonMesh.rotation.y = -th + Math.PI;
    /* L1/L2 标记随会合系实时旋转 */
    if (lpointGroup.visible) {
      const tL = tSincePeri;
      for (const mk of lpointGroup.children) {
        if (mk.name === 'L1' || mk.name === 'L2') {
          const rL = X3.lpointDisplay(mk.name, tL);
          mk.position.copy(eciToScene(rL));
        } else if (mk.name.endsWith('-label')) {
          const pt = mk.name.split('-')[0];
          const rL = X3.lpointDisplay(pt, tL);
          placeLabel(mk, eciToScene(rL), 0.5);
        }
      }
    }
    /* 「地球」标签直接锚在真实地球球体上沿（球体在 384 400 km ≈ 384.4 场景单位外，
       视角放大到月心附近时在视野外属正常；拉远即见）。dy=0 绕过 placeLabel 的
       labelScale 缩放：球体是固定尺寸，偏移必须固定，否则拉远后标签飘离球体 */
    if (labels.earth) placeLabel(labels.earth,
      ePos.clone().add(new THREE.Vector3(0, EARTH_MINI_R + 1.1, 0)), 0);
  }

  /* L2 会合系视角：会合系静态环（θ=0 锚定）按 FRAME_AXIS 实时刚性旋转，环心恒压在
     实时 L2 上（恒等式 disp(θ) = R(轴, ω·t)·disp(0) 已数值验证，残差 <0.1 km） */
  if (orbitLineSyn) orbitLineSyn.quaternion.setFromAxisAngle(FRAME_AXIS, X3.OMEGA_EM * tSincePeri);

  /* ν 角度弧每帧跟随：近地点方向 → 当前航天器方向（仅覆写顶点，不重建对象） */
  if (angleInfo && nuArc) {
    const { h, periDir, v, rV } = angleInfo;
    const theta = angleAbout(periDir, p.clone().normalize(), h);
    const pa = nuArc.geometry.attributes.position;
    for (let k = 0; k <= ARC_SEG; k++) {
      const t = theta * k / ARC_SEG, ct = Math.cos(t), st2 = Math.sin(t);
      pa.setXYZ(k,
        rV * (periDir.x * ct + v.x * st2),
        rV * (periDir.y * ct + v.y * st2),
        rV * (periDir.z * ct + v.z * st2));
    }
    pa.needsUpdate = true;
    if (labels.nu) placeArcLabel(labels.nu, arcMidDir(periDir, v, theta), rV);
  }

  /* 摄动模式：瞬时密切椭圆随 tNow 岁差进动，每帧覆写 orbitLine 顶点
     （拖尾只显 tNow 前 7 d 尾迹，当前/未来轨道靠这条椭圆补全，避免画面空）
     CR3BP 模式轨道线是完整周期轨迹，无需覆写 */
  if (prop && prop.isPerturbed && !prop.isCr3bp && orbitLine && st.osc) {
    const ol = orbitLine.geometry.attributes.position;
    const mu = body === 'moon' ? X3.MU_MOON : undefined;
    const br = body === 'moon' ? X3.R_MOON : undefined;
    for (let k = 0; k <= 256; k++) {
      const nu = (k / 256) * C.TWO_PI;
      const q = eciToScene(C.classicalToState(st.osc, nu, mu, br).rE);
      ol.setXYZ(k, q.x, q.y, q.z);
    }
    ol.needsUpdate = true;
  }

  const vScene = eciToScene(st.vE).normalize();   /* 实时显示系速度方向（L2 视角同理） */
  const rScene = p.clone().normalize();
  scState = { p, vScene, rScene };   // 供第一视角相机取用

  /* 缩放自适应：航天器 / 近远月点标记 / 箭头随相机视距同比缩放（透视下屏占比恒定，
     标签由 labelScale 内的 zoomK 覆盖；水平盘尺寸由箭头长度 len 派生自动跟随） */
  const zk = zoomK();
  scGroup.scale.setScalar(zk);
  if (periMark) { periMark.scale.setScalar(zk); apoMark.scale.setScalar(zk); }
  /* 静态标签（rebuild 时登记锚点）每帧重放，scale 与偏移同步跟随 zoomK */
  for (const key in labels) {
    const sp = labels[key];
    if (!sp) continue;
    if (sp.userData.anchor) placeLabel(sp, sp.userData.anchor.pos, sp.userData.anchor.dy);
    else if (sp.userData.arcAnchor) placeArcLabel(sp, sp.userData.arcAnchor.midDir, sp.userData.arcAnchor.r);
  }

  /* 航天器 caption 跟随（浮在航天器上方） */
  if (labels.sc) placeLabel(labels.sc, p, 1.35);

  const len = Math.max(0.35, st.v * 0.45) * zk;

  velArrow.position.copy(p);
  velArrow.setDirection(vScene);
  velArrow.setLength(len, Math.max(0.14 * zk, len * 0.16), Math.max(0.08 * zk, len * 0.08));
  rArrow.position.copy(p);
  rArrow.setDirection(rScene);
  rArrow.setLength(len, Math.max(0.14 * zk, len * 0.16), Math.max(0.08 * zk, len * 0.08));

  /* CR3BP 受力分析三箭头（可见性由 applyDisplay 集中管理） */
  if (forceGroup.visible) updateForces(p, st, zk);

  disc.position.copy(p);
  disc.quaternion.setFromUnitVectors(Z_AXIS, rScene);
  const dr = Math.max(0.9 * zk, len * 0.55);
  disc.scale.setScalar(dr);

  const pos = posLine.geometry.attributes.position;
  pos.setXYZ(0, 0, 0, 0);
  pos.setXYZ(1, p.x, p.y, p.z);
  pos.needsUpdate = true;
  posLine.computeLineDistances();
}

/* ── 文本与推导面板 ──────────────────────────────────────────────── */
function mtxTable(m) {
  let h = '<table class="mtx">';
  for (let r = 0; r < 3; r++) {
    h += '<tr>';
    for (let c = 0; c < 3; c++) h += `<td>${m[r * 3 + c].toFixed(5)}</td>`;
    h += '</tr>';
  }
  return h + '</table>';
}

function derivPos(st, el) {
  const rName = body === 'moon' ? '月心距' : '地心距';
  return `
  <div class="step"><b>① 半通径</b>　p = a(1 − e²) = ${fmt(el.a)} × (1 − ${fmt(el.e, 4)}²)
    = <span class="res">${fmt(st.p)} km</span></div>
  <div class="step"><b>② ${rName}</b>　r = p / (1 + e·cos ν) = ${fmt(st.p)} / (1 + ${fmt(el.e, 4)}·cos ${fmtDeg(el.nu)})
    = <span class="res">${fmt(st.r)} km</span></div>
  <div class="step"><b>③ 近焦点系 (PQW)</b>　r<sub>PQW</sub> = (r·cos ν, r·sin ν, 0)
    = <span class="res">${fmtVec(st.rP, 'km')}</span></div>
  <div class="step"><b>④ 旋转矩阵</b>　R = R₃(−Ω)·R₁(−i)·R₃(−ω)（PQW → ECI）${mtxTable(st.matrix)}</div>
  <div class="step"><b>⑤ 惯性系位置</b>　r<sub>ECI</sub> = R·r<sub>PQW</sub>
    = <span class="res">${fmtVec(st.rE, 'km')}</span></div>
  <div class="note">${prop && prop.isPerturbed ? '摄动模式：上式为瞬时密切根数下的二体位置（密切椭圆在 r<sub>ECI</sub> 处与真实摄动轨迹相切）。' : '二体：r<sub>ECI</sub> 即真实位置（Δr ≡ 0）。'}</div>`;
}

function derivAlt(st, el) {
  const isMoon = body === 'moon';
  const RB = isMoon ? X3.R_MOON : C.R_EARTH;      /* 中心天体半径（曾硬编码 R⊕，月球选项卡数值错误） */
  const RName = isMoon ? 'R☾' : 'R⊕';
  const hp = el.a * (1 - el.e) - RB;
  const ha = el.a * (1 + el.e) - RB;
  return `
  <div class="step"><b>① 球体近似高度</b>　h = r − ${RName} = ${fmt(st.r)} − ${fmt(RB, 3)}
    = <span class="res">${fmt(st.alt)} km</span></div>
  <div class="step"><b>② 常数</b>　${RName} = ${fmt(RB, 3)} km（${isMoon ? '月球平均半径' : 'WGS-84 半长轴的球体近似'}；
    未做椭球 / 大地水准高度修正）</div>
  <div class="step"><b>③ 拱点高度</b>　h<sub>p</sub> = a(1 − e) − ${RName} = <span class="res">${fmt(hp)} km</span>，
    h<sub>a</sub> = a(1 + e) − ${RName} = <span class="res">${fmt(ha)} km</span></div>`;
}

function derivVel(st, el) {
  const isMoon = body === 'moon';
  const mu = isMoon ? X3.MU_MOON : C.MU;
  const bodyZh = isMoon ? '月球' : '地球';
  return `
  <div class="step"><b>① 活力公式</b>　v = √( μ(2/r − 1/a) )
    = √( ${mu} × (2/${fmt(st.r)} − 1/${fmt(el.a)}) )
    = <span class="res">${fmt(st.v)} km/s</span></div>
  <div class="step"><b>② 径向 / 横向分解</b>　v<sub>r</sub> = √(μ/p)·e·sin ν
    = <span class="res">${fmt(st.vr)} km/s</span>（沿 r̂，远离${bodyZh}为正）<br>
    v<sub>t</sub> = √(μ/p)·(1 + e·cos ν) = <span class="res">${fmt(st.vt)} km/s</span>（沿飞行方向横向）</div>
  <div class="step"><b>③ 校验</b>　√(v<sub>r</sub>² + v<sub>t</sub>²) = ${fmt(Math.hypot(st.vr, st.vt))} ≡ v ✓</div>
  <div class="step"><b>④ 近焦点系速度</b>　v<sub>PQW</sub> = √(μ/p)·(−sin ν, e + cos ν, 0)
    = ${fmtVec(st.vP, 'km/s')}</div>
  <div class="step"><b>⑤ 惯性系速度</b>　v<sub>ECI</sub> = R·v<sub>PQW</sub>
    = <span class="res">${fmtVec(st.vE, 'km/s')}</span></div>
  <div class="note">${prop && prop.isPerturbed ? '摄动模式：上式为瞬时密切椭圆的二体速度，与真实 v<sub>ECI</sub> 相切。' : '二体：v<sub>ECI</sub> 即真实速度（Δv ≡ 0）。'}</div>`;
}

function derivGam(st, el) {
  const g = st.gamma * C.R2D;
  return `
  <div class="step"><b>① 定义式</b>　γ = atan2( e·sin ν, 1 + e·cos ν )
    = atan2( ${fmt(el.e * Math.sin(el.nu * C.D2R), 4)}, ${fmt(1 + el.e * Math.cos(el.nu * C.D2R), 4)} )
    = <span class="res">${(g >= 0 ? '+' : '') + fmt(g, 3)}°</span></div>
  <div class="step"><b>② 等价关系</b>　tan γ = v<sub>r</sub> / v<sub>t</sub>
    = ${fmt(st.vr)} / ${fmt(st.vt)}（速度矢量与当地水平面 ⊥ r 的夹角）</div>
  <div class="step"><b>③ 符号约定</b>　γ &gt; 0：离开${body === 'moon' ? '月球' : '地球'}（升弧段 ν∈(0°,180°)）；
    γ &lt; 0：返回${body === 'moon' ? '月球' : '地球'}（降弧段）；拱点（ν = 0° / 180°）处 γ = 0。
    3D 视图中黄色速度箭头相对蓝色圆盘的倾角即 γ。</div>`;
}

/* ── 摄动加速度推导（仅摄动模式可见） ──────────────────────────── */
function m2s(v) { return v < 1e-12 ? '<1e-12' : v.toExponential(2); }

function derivZonal(st, el, b) {
  const EP = window.ElementsPerturb;
  const rm = st.r, rE = st.rE;
  const s = rE[2] / rm;                 // sinφ = z/r
  const Re_r = C.R_EARTH / rm;
  const s2 = s * s;
  const P2 = (3 * s2 - 1) / 2;
  const P3 = (5 * s2 * s - 3 * s) / 2;
  const P4 = (35 * s2 * s2 - 30 * s2 + 3) / 8;
  const t2 = Re_r * Re_r, t3 = t2 * Re_r, t4 = t3 * Re_r;
  const sumP = EP.J2 * t2 * P2 + EP.J3 * t3 * P3 + EP.J4 * t4 * P4;
  const a = b.zonal;
  const mag = Math.hypot(a[0], a[1], a[2]) * 1000;
  return `
  <div class="step"><b>① 位函数</b>　U = μ/r·[1 − Σ Jₙ(R⊕/r)ⁿ Pₙ(sinφ)]<br>
    sinφ = z/r = ${fmt(rE[2], 3)}/${fmt(rm, 3)} = ${fmt(s, 5)}（φ 为地心纬度）</div>
  <div class="step"><b>② 带谐系数</b>　J2 = ${EP.J2.toExponential(4)} · J3 = ${EP.J3.toExponential(4)}
    · J4 = ${EP.J4.toExponential(4)}（EGM2008 未归一化）<br>
    R⊕/r = ${fmt(C.R_EARTH, 3)}/${fmt(rm, 3)} = ${fmt(Re_r, 5)}</div>
  <div class="step"><b>③ Legendre 多项式</b>　P₂(s) = (3s²−1)/2 = ${fmt(P2, 5)}<br>
    P₃(s) = (5s³−3s)/2 = ${fmt(P3, 5)} · P₄(s) = (35s⁴−30s²+3)/8 = ${fmt(P4, 5)}</div>
  <div class="step"><b>④ 摄动位</b>　Σ Jₙ(R⊕/r)ⁿ Pₙ = J2·(R⊕/r)²·P2 + J3·(R⊕/r)³·P3 + J4·(R⊕/r)⁴·P4
    = <span class="res">${sumP.toExponential(4)}</span></div>
  <div class="step"><b>⑤ 加速度</b>　a = ∇U − a<sub>中心</sub>（μr̂/r²），ECI 分量
    = <span class="res">${fmtVec([a[0] * 1000, a[1] * 1000, a[2] * 1000], 'm/s²')}</span><br>
    |a<sub>zonal</sub>| = <span class="res">${m2s(mag)}</span></div>
  <div class="note">J2 项主导，J3/J4 量级 1e-6。梯度对 r 与 s 分别求偏导后相加。</div>`;
}

function derivDrag(st, el, b) {
  const EP = window.ElementsPerturb;
  const rm = st.r, rE = st.rE, vE = st.vE;
  const altKm = rm - C.R_EARTH;
  const rho = EP.density(altKm);        // kg/m³
  const OMEGA = EP.OMEGA_E;
  const wxr = [-OMEGA * rE[1], OMEGA * rE[0], 0];   // ω⊕×r，km/s
  const vrel = [vE[0] - wxr[0], vE[1] - wxr[1], vE[2] - wxr[2]]; // km/s
  const vrel_m = vrel.map(x => x * 1000);            // m/s
  const vrelm = Math.hypot(vrel_m[0], vrel_m[1], vrel_m[2]);
  const o = prop.o;
  const Cd = o.Cd, SM = o.SM;
  const k = 0.5 * Cd * SM * rho * vrelm;            // m/s² per (m/s)
  const a = b.drag;
  const mag = Math.hypot(a[0], a[1], a[2]) * 1000;
  return `
  <div class="step"><b>① 大气密度</b>　ρ(h) = ρ₀·exp(−(h−h₀)/H)（分段指数模型）<br>
    h = ${fmt(altKm, 1)} km → ρ = <span class="res">${rho < 1e-15 ? '<1e-15' : rho.toExponential(3)}</span> kg/m³</div>
  <div class="step"><b>② 相对速度</b>　v<sub>rel</sub> = v − ω⊕×r（大气随地球刚性自转）<br>
    ω⊕ = (0, 0, ${OMEGA.toExponential(3)}) rad/s<br>
    v<sub>rel</sub> = <span class="res">${fmtVec(vrel_m, 'm/s')}</span>，|v<sub>rel</sub>| = ${fmt(vrelm, 1)} m/s</div>
  <div class="step"><b>③ 阻力公式</b>　a = −½·C<sub>D</sub>·(S/m)·ρ·|v<sub>rel</sub>|·v<sub>rel</sub><br>
    C<sub>D</sub> = ${Cd} · S/m = ${SM} m²/kg<br>
    ½·C<sub>D</sub>·(S/m)·ρ·|v<sub>rel</sub>| = ${k.toExponential(4)} m/s² per (m/s)</div>
  <div class="step"><b>④ ECI 加速度</b>　a<sub>drag</sub> = <span class="res">${fmtVec([a[0] * 1000, a[1] * 1000, a[2] * 1000], 'm/s²')}</span><br>
    |a<sub>drag</sub>| = <span class="res">${m2s(mag)}</span></div>
  <div class="note">方向与 v<sub>rel</sub> 反向（负号）；结果 ÷1000 化为 km/s² 后叠加进 RK4 右端。</div>`;
}

function derivSrp(st, el, b) {
  const EP = window.ElementsPerturb;
  const o = prop.o;
  const tSec = prop.tNow;
  const sHat = EP.sunDir(tSec, o.lam0);
  const rE = st.rE;
  const nu = EP.shadowFactor(rE, sHat);
  const rdotS = rE[0] * sHat[0] + rE[1] * sHat[1] + rE[2] * sHat[2];
  const r2 = rE[0] * rE[0] + rE[1] * rE[1] + rE[2] * rE[2];
  const d2 = r2 - rdotS * rdotS;
  const a = b.srp;
  const mag = Math.hypot(a[0], a[1], a[2]) * 1000;
  const PSUN = EP.P_SUN;
  return `
  <div class="step"><b>① 太阳方向</b>　λ = λ₀ + 2πt/年（太阳视运动匀速圆轨道，黄赤交角 23.44°）<br>
    ŝ = <span class="res">${fmtVec(sHat, '')}</span>（ECI 单位矢量）</div>
  <div class="step"><b>② 圆柱形地影</b>　r·ŝ = ${fmt(rdotS, 1)} km（&lt; 0 为背阳侧）<br>
    到日地轴垂距² = r² − (r·ŝ)² = ${fmt(d2, 1)} km²，R⊕² = ${fmt(C.R_EARTH * C.R_EARTH, 1)} km²<br>
    受晒系数 ν = <span class="res">${nu}</span>（本影内 0，外部 1，无半影过渡）</div>
  <div class="step"><b>③ 光压常数</b>　p⊙ = ${PSUN.toExponential(2)} N/m²<br>
    C<sub>R</sub> = ${o.CR}（= 1 + η，η = ${(o.CR - 1).toFixed(2)}）· S/m = ${o.SM} m²/kg</div>
  <div class="step"><b>④ 光压加速度</b>　a = −ν·(S/m)·C<sub>R</sub>·p⊙·ŝ<br>
    a<sub>SRP</sub> = <span class="res">${fmtVec([a[0] * 1000, a[1] * 1000, a[2] * 1000], 'm/s²')}</span><br>
    |a<sub>SRP</sub>| = <span class="res">${m2s(mag)}</span></div>
  <div class="note">方向背向太阳（负号）；地影内 ν=0 时光压为零。结果 ÷1000 化为 km/s²。</div>`;
}

/* ── CR3BP 雅可比积分推导 ───────────────────────────────────── */
function derivCj(st, mp) {
  const X = X3;
  const rot = st.rot;
  const r1 = st.r1, r2 = st.r2;
  const v2 = rot.vx * rot.vx + rot.vy * rot.vy + rot.vz * rot.vz;
  const Om = X.omegaPot(rot.x, rot.y, rot.z);
  const cj = 2 * Om - v2;
  return `
  <div class="step"><b>① 无量纲化</b>　特征长度 l = a_EM = ${fmt(X.A_EM, 0)} km · 特征时间 1/ω_EM（ω_EM = √(G(m⊕+m☾)/l³)）<br>
    质量比 μ = m☾/(m⊕+m☾) = ${X.MU_RATIO.toFixed(8)} · P₁(地球) 在 (−μ,0,0)，P₂(月球) 在 (1−μ,0,0)</div>
  <div class="step"><b>② 有效势</b>　Ω = (x²+y²)/2 + (1−μ)/r₁ + μ/r₂<br>
    r₁ = √((x+μ)²+y²+z²) = ${r1.toFixed(6)} · r₂ = √((x−1+μ)²+y²+z²) = ${r2.toFixed(6)}<br>
    Ω = ${Om.toFixed(6)}</div>
  <div class="step"><b>③ 雅可比积分</b>　C = 2Ω − v²（会合系速度）<br>
    v² = ${v2.toFixed(8)} · C = 2×${Om.toFixed(6)} − ${v2.toFixed(6)} = <span class="res">${cj.toFixed(6)}</span></div>
  <div class="step"><b>④ 零速度面对照</b>　平动点阈值：C(L1) = ${mp.cjL1.toFixed(4)} · C(L2) = ${mp.cjL2.toFixed(4)} · C(L3) = ${mp.cjL3.toFixed(4)}<br>
    当前 C ${cj < mp.cjL2 ? '&lt; C(L2)：L2 颈口已开——可经 L2 通往月球背侧外部空间' : '≥ C(L2)：L2 颈口闭合'}<br>
    当前 C ${cj < mp.cjL1 ? '&lt; C(L1)：L1 颈口已开——地月两侧空间连通' : '≥ C(L1)：L1 颈口闭合'}</div>
  <div class="step"><b>⑤ 平动点坐标</b>　共线点解 ∂Ω/∂x = 0（五次方程二分）：<br>
    L1 x = ${mp.xL1.toFixed(6)} l（距月心 ${fmt((mp.xL1 - 1 + X.MU_RATIO) * X.A_EM, 0)} km）<br>
    L2 x = ${mp.xL2.toFixed(6)} l（距月心 ${fmt((mp.xL2 - 1 + X.MU_RATIO) * X.A_EM, 0)} km）</div>
  <div class="note">会合系中 C 沿轨迹守恒（数值验证漂移 &lt; 1e-8）；惯性系能量不守恒（旋转系非惯性）。</div>`;
}

/* ── CR3BP 速度推导（三体无守恒二体能量，活力公式不适用） ─────────── */
function derivVelCr3bp(st) {
  const X = X3;
  const lw = X.A_EM * X.OMEGA_EM;              // 特征速度 l·ω，km/s
  const rot = st.rot;
  const vrot = Math.hypot(rot.vx, rot.vy, rot.vz);
  return `
  <div class="step"><b>① 为什么不是活力公式</b>　三体问题无守恒的二体能量，
    v = √(μ(2/r − 1/a)) 不成立（瞬时密切椭圆的 a 可正可负，仅相切不代表能量）。
    速度由会合系状态经旋转变换直接求得。</div>
  <div class="step"><b>② 会合系读数</b>　ρ = (${rot.x.toFixed(4)}, ${rot.y.toFixed(4)}, ${rot.z.toFixed(4)}) l<br>
    ρ̇ = (${rot.vx.toFixed(4)}, ${rot.vy.toFixed(4)}, ${rot.vz.toFixed(4)}) l/T（|ρ̇| = ${vrot.toFixed(4)}）</div>
  <div class="step"><b>③ 惯性系变换</b>　v = l·ω_EM·[R′(θ)·ρ + R(θ)·ρ̇] − v<sub>moon</sub>(θ)<br>
    特征速度 l·ω_EM = ${fmt(lw, 4)} km/s · θ = ω_EM·t（月球公转相位）<br>
    v<sub>ECI</sub> = <span class="res">${fmtVec(st.vE, 'km/s')}</span></div>
  <div class="step"><b>④ 速度大小</b>　|v| = <span class="res">${fmt(st.v)} km/s</span>（月心惯性系）<br>
    径向 v<sub>r</sub> = ${fmt(st.vr)} km/s · 横向 v<sub>t</sub> = ${fmt(st.vt)} km/s（沿月心 r̂ 分解）</div>
  <div class="note">量级参考：NRHO 近月点 |v| ≈ 0.8 km/s；L2 Halo 月心惯性速度约 0.1–0.5 km/s。
    雅可比常数卡片中的 v_rot 为会合系速度（另一坐标系读数，不可直接比较）。</div>`;
}

function updateTexts(st, nuRad) {
  const isCr3bp = !!(prop && prop.isCr3bp);
  const isMoon = body === 'moon';
  const bodyName = isMoon ? '月' : '地';
  const bodyRName = isMoon ? 'R☾' : 'R⊕';
  $('pos-val').textContent = fmtVec(st.rE, 'km');
  $('pos-sub').textContent = `${isMoon ? '月心距' : '地心距'} |r| = ${fmt(st.r)} km`;
  $('alt-val').textContent = `${fmt(st.alt)} km`;
  $('alt-sub').textContent = st.alt < 0
    ? `⚠ 当前位置在${isMoon ? '月球' : '地球'}内部（r < ${bodyRName}）`
    : `球体近似（r − ${bodyRName}）· ${isMoon ? '月心距' : '地心距'} r = ${fmt(st.r)} km`;

  /* 月面高度卡片：按半长轴 a 推荐传播模型（仅月球选项卡、根数模式；
     CR3BP 预设轨道的 els.a 为占位值无判据意义，隐藏）。
     分界：a < 100+R☾ → 普通二体；100+R☾ < a < 10000+R☾ → 摄动二体；a > 10000+R☾ → CR3BP。
     当前模型与推荐不符时换警示色提示。 */
  const modelHint = $('alt-model-hint');
  /* ELFO 冻结轨道的设计核心即地球三体引力：摄动二体（月球 J2–J4 + 地球三体）
     与 CR3BP 两档模型精度均足够，不再弹「推荐模型」警示（2026-10-15 用户要求）；
     普通二体仍提示升级。仅当仍标称 ELFO 预设时抑制（拖根数后 moonPresetKey 清空）。 */
  if (isMoon && !cr3bpState.active
      && !(moonPresetKey === 'elfo' && ['perturbed-moon', 'cr3bp'].includes($('sel-prop').value))) {
    const RM = X3.R_MOON;
    const [recId, recName] = els.a < 100 + RM ? ['kepler-moon', '普通二体']
      : els.a < 10000 + RM ? ['perturbed-moon', '摄动二体']
      : ['cr3bp', 'CR3BP'];
    const cur = $('sel-prop').value;
    const ok = cur === recId;
    modelHint.style.display = '';
    modelHint.style.color = ok ? '' : 'var(--warn)';
    modelHint.innerHTML = `${ok ? '✓' : '⚠'} a = ${fmt(els.a, 0)} km（距月面 ≈ ${fmt(els.a - RM, 0)} km）→ 推荐：${recName}`
      + (ok ? '' : `（当前：${C.propagators[cur] ? C.propagators[cur].label : cur}）`);
  } else {
    modelHint.style.display = 'none';
  }
  $('vel-val').textContent = `${fmt(st.v)} km/s`;
  $('vel-sub').innerHTML = `v = ${fmtVec(st.vE, 'km/s')} · v<sub>r</sub> = ${fmt(st.vr)} · v<sub>t</sub> = ${fmt(st.vt)} km/s`;
  const g = st.gamma * C.R2D;
  $('gam-val').textContent = `${g >= 0 ? '+' : ''}${fmt(g, 3)}°`;
  $('gam-sub').innerHTML = `${g > 0.05 ? `升弧段：正在离开${isMoon ? '月球' : '地球'}` : g < -0.05 ? `降弧段：正在返回${isMoon ? '月球' : '地球'}` : '拱点附近：γ ≈ 0'} · v<sub>r</sub> = ${fmt(st.vr)} km/s`;

  /* 偏近点角 / 平近点角（随当前 ν 实时求解，归约到 0–360°；CR3BP 无定义隐藏） */
  if (isCr3bp || !isFinite(st.E)) {
    $('d-E').textContent = '—';
    $('d-M').textContent = '—';
  } else {
    $('d-E').textContent = fmtDeg(C.wrap2pi(st.E) * C.R2D);
    $('d-M').textContent = fmtDeg(C.wrap2pi(st.M) * C.R2D);
  }

  /* CR3BP 卡片：雅可比积分 + 会合系读数 + 使命轨道参数 */
  if (isCr3bp && st.rot) {
    const cj = st.cj;
    const mp = prop.missionParams();
    const open1 = cj < mp.cjL1, open2 = cj < mp.cjL2;
    $('cj-val').textContent = `C_J = ${cj.toFixed(6)}`;
    $('cj-sub').innerHTML =
      `会合系坐标 (x,y,z) = (${st.rot.x.toFixed(4)}, ${st.rot.y.toFixed(4)}, ${st.rot.z.toFixed(4)}) l<br>`
      + `r₁ = ${(st.r1).toFixed(4)} l · r₂ = ${(st.r2).toFixed(4)} l · v_rot = ${Math.hypot(st.rot.vx, st.rot.vy, st.rot.vz).toFixed(4)} l/T<br>`
      + `零速度面：L1 阈值 C=${mp.cjL1.toFixed(4)}（${open1 ? '已开' : '闭'}） · L2 阈值 C=${mp.cjL2.toFixed(4)}（${open2 ? '已开' : '闭'}）`;
    $('mission-sub').innerHTML =
      `周期 T = ${(mp.T_s / 86400).toFixed(3)} d · 近月点 ${fmt(mp.periluneKm)} km · 远月点 ${fmt(mp.apoluneKm)} km<br>`
      + `振幅 Az = ${fmt(mp.azKm)} km · Ax = ${fmt(mp.axKm)} km · Ay = ${fmt(mp.ayKm)} km<br>`
      + `平动点 L1 x = ${mp.xL1.toFixed(6)} l · L2 x = ${mp.xL2.toFixed(6)} l（无量纲质心系）`;
    if ($('cj-deriv').classList.contains('show')) {
      $('cj-deriv').innerHTML = derivCj(st, mp);
    }
  }

  /* 相位链（暂停 = 从 ν 反算；播放 = 从 t 正算；CR3BP 显示相位代理） */
  const n = isCr3bp
    ? C.TWO_PI / prop.periodSeconds()          /* 三体：由周期轨道 T 直接给平均角速率（els 为占位根数） */
    : C.meanMotion(els.a, isMoon ? X3.MU_MOON : undefined);
  $('ph-n').textContent = `${fmt(n * C.R2D, 4)} °/s`;
  $('ph-T').textContent = fmtT(prop.periodSeconds());
  const arrow = '<span class="chain-arrow">→</span>';
  if (isCr3bp) {
    const frac = ((tSincePeri / prop.periodSeconds()) % 1 + 1) % 1;
    $('ph-chain').innerHTML =
      `三体模式：会合系相位 = t/T·2π ${arrow} 无开普勒链（非二体）<br>`
      + `相位 = ${fmtDeg(frac * C.TWO_PI * C.R2D)} · T = ${fmtT(prop.periodSeconds())}`;
  } else if (playing) {
    const sol = C.mToE(n * tSincePeri, els.e);
    $('ph-chain').innerHTML =
      `正向链：M = n·Δt ${arrow} E（牛顿迭代 ${sol.iters} 次）${arrow} ν<br>`
      + `M = ${fmtDeg(C.wrap2pi(st.M) * C.R2D)} · E = ${fmtDeg(st.E * C.R2D)} · ν = ${fmtDeg(els.nu)}`;
  } else {
    $('ph-chain').innerHTML =
      `反算链：ν ${arrow} E（半角公式）${arrow} M = E − e·sinE<br>`
      + `M = ${fmtDeg(C.wrap2pi(st.M) * C.R2D)} · E = ${fmtDeg(st.E * C.R2D)} · ν = ${fmtDeg(els.nu)}`;
  }
  const tMod = ((tSincePeri % prop.periodSeconds()) + prop.periodSeconds()) % prop.periodSeconds();
  const tShow = (prop && prop.isPerturbed) ? prop.tNow : tMod;
  /* 时间显示统一用小时（90 d 窗口 = 2160 h） */
  $('ph-t').textContent = (prop && prop.isPerturbed)
    ? `历元起算 t = ${(tShow / 3600).toFixed(1)} h / ${(prop.windowSeconds() / 3600).toFixed(0)} h`
    : `距上次过${isMoon ? '近月点' : '近地点'} Δt = ${fmtT(tShow)}`;
  if (prop && prop.isPerturbed && !isCr3bp) {
    updateOscTable(st);
    /* 图表指针随 tNow 同步（播放/拖动进度条时）；tNow 未变则跳过，省重绘 */
    if (prop.tNow !== lastChartTNow) { drawCharts(); lastChartTNow = prop.tNow; }
  }
  /* 进度条：摄动与 CR3BP 都随 tNow 刷新 */
  if (prop && prop.isPerturbed) updateProgBar();

  /* 摄动加速度卡片（仅地球摄动模式可见）。
     bd = 当前 tNow 处的加速度分解（与 st 同源 stateVecAt(tNow)，数值一致）。
     必须排除月心/CR3BP：Cr3bpPropagator 未设 this.o，本块曾以
     prop.o.lam0 抛 TypeError 使 tick 主循环死亡（3D 画面永久冻结）。 */
  const isP = !!(prop && prop.isPerturbed && body === 'earth');
  $('card-fz').style.display = isP ? '' : 'none';
  $('card-fd').style.display = isP ? '' : 'none';
  $('card-fs').style.display = isP ? '' : 'none';
  let bd = null;
  if (isP) {
    const EP = window.ElementsPerturb;
    bd = EP.accelBreakdown(prop.stateVecAt(prop.tNow), prop.tNow, prop.o);
    const gravOn = $('pp-grav').checked, dragOn = $('pp-drag').checked, srpOn = $('pp-srp').checked;
    const fzMag = gravOn ? Math.hypot(bd.zonal[0], bd.zonal[1], bd.zonal[2]) * 1000 : 0;
    const fdMag = dragOn ? Math.hypot(bd.drag[0], bd.drag[1], bd.drag[2]) * 1000 : 0;
    const fsMag = srpOn ? Math.hypot(bd.srp[0], bd.srp[1], bd.srp[2]) * 1000 : 0;
    const rho = EP.density(st.r - C.R_EARTH);
    const nu = EP.shadowFactor(st.rE, EP.sunDir(prop.tNow, prop.o.lam0));
    $('fz-val').textContent = gravOn ? m2s(fzMag) + ' m/s²' : '（关闭）';
    $('fz-sub').innerHTML = gravOn
      ? `J2+J3+J4 · ECI = ${fmtVec([bd.zonal[0] * 1000, bd.zonal[1] * 1000, bd.zonal[2] * 1000], 'm/s²')}`
      : '地球非球形引力已关闭';
    $('fd-val').textContent = dragOn ? m2s(fdMag) + ' m/s²' : '（关闭）';
    $('fd-sub').innerHTML = dragOn
      ? `½·C_D·(S/m)·ρ·|v_rel| · ρ(${fmt(st.r - C.R_EARTH, 0)} km) = ${rho < 1e-15 ? '<1e-15' : rho.toExponential(2)} kg/m³`
      : '大气阻力已关闭';
    $('fs-val').textContent = srpOn ? m2s(fsMag) + ' m/s²' : '（关闭）';
    $('fs-sub').innerHTML = srpOn
      ? `C_R·(S/m)·p⊙·ν · 受晒系数 ν = ${nu}`
      : '太阳光压已关闭';
  }

  /* 推导面板（仅展开时渲染）。
     摄动模式下 pos/alt/vel/gam 推导需用密切根数重算二体状态（p/rP/vP/matrix），
     因为 PerturbedPropagator.stateAtTime 只返回 rE/vE 等不含 PQW 字段；
     直接传 st 会使 mtxTable(undefined) 抛 TypeError。
     此处对月心 Cowell / CR3BP 同样生效（不能复用上面的 isP，其已限地球摄动）。 */
  let derivSt = st, derivEl = els;
  if (prop && prop.isPerturbed && st.osc) {
    derivSt = C.classicalToState(st.osc, st.osc.nu);
    derivEl = { a: st.osc.a, e: st.osc.e, i: st.osc.i * C.R2D,
                Om: st.osc.Om * C.R2D, w: st.osc.w * C.R2D, nu: st.osc.nu * C.R2D };
  }
  if ($('pos-deriv').classList.contains('show')) $('pos-deriv').innerHTML = derivPos(derivSt, derivEl);
  if ($('alt-deriv').classList.contains('show')) $('alt-deriv').innerHTML = derivAlt(derivSt, derivEl);
  if ($('vel-deriv').classList.contains('show')) {
    /* CR3BP：活力公式不适用（无守恒二体能量），用会合系→惯性系变换推导 */
    $('vel-deriv').innerHTML = isCr3bp ? derivVelCr3bp(st) : derivVel(derivSt, derivEl);
  }
  if ($('gam-deriv').classList.contains('show')) $('gam-deriv').innerHTML = derivGam(derivSt, derivEl);
  if (isP && bd) {
    if ($('fz-deriv').classList.contains('show')) $('fz-deriv').innerHTML = derivZonal(st, derivEl, bd);
    if ($('fd-deriv').classList.contains('show')) $('fd-deriv').innerHTML = derivDrag(st, derivEl, bd);
    if ($('fs-deriv').classList.contains('show')) $('fs-deriv').innerHTML = derivSrp(st, derivEl, bd);
  }

  /* 播放时同步 ν 输入控件 */
  if (playing) syncInputs('nu');
}

/* ── 控件接线 ────────────────────────────────────────────────────── */
const rows = {};
document.querySelectorAll('.el-row').forEach((row) => {
  const key = row.dataset.key;
  const [min, max, step] = RANGES[key];
  const rng = row.querySelector('input[type="range"]');
  const num = row.querySelector('input[type="number"]');
  rng.min = min; rng.max = max; rng.step = step;
  num.min = min; num.max = max;
  rows[key] = { rng, num };

  function apply(val, fromNum) {
    if (key === 'e' && val >= 1) {          // 双曲线/抛物线拒绝
      num.classList.add('invalid');
      return;
    }
    num.classList.remove('invalid');
    const [cMin, cMax] = RANGES[key];   /* clamp 用当前天体范围（地/月不同，随选项卡切换） */
    const v = Math.min(cMax, Math.max(cMin, val));
    if (key === 'nu' && playing) { setPlaying(false); }
    els[key] = v;
    /* 手动改根数形状（ν 仅相位）即脱离标称预设：ELFO 专属警示抑制解除 */
    if (key !== 'nu') moonPresetKey = null;
    if (!fromNum) num.value = fmt(v, 4).replace(/,/g, ''); else rng.value = v;
    if (key === 'nu') { tSincePeri = prop.timeAtTrueAnomaly(v * C.D2R); }
    else rebuild();
  }
  rng.addEventListener('input', () => apply(parseFloat(rng.value), false));
  num.addEventListener('input', () => {
    const s = num.value.trim();
    if (s === '' || s === '-' || s === '.') return;
    const v = parseFloat(s);
    if (!isFinite(v)) return;
    apply(v, true);
  });
  num.addEventListener('blur', () => syncInputs(key));
});
function syncInputs(key) {
  const { rng, num } = rows[key];
  rng.value = els[key];
  num.value = +els[key].toFixed(4);
  num.classList.remove('invalid');
}
function syncAll() { Object.keys(rows).forEach(syncInputs); }

/* 选项卡切换：把当前天体范围整体刷到六根数滑块 / 数值框（min/max/step）。
   此前月球选项卡沿用地球范围（a 最大 20 万 km），相对月球过大——
   根因即 RANGES 固化在 earth 条目、BODY_RANGES.moon 从未生效。 */
function applyBodyRanges() {
  RANGES = BODY_RANGES[body];
  document.querySelectorAll('.el-row').forEach((row) => {
    const [min, max, step] = RANGES[row.dataset.key];
    const rng = row.querySelector('input[type="range"]');
    const num = row.querySelector('input[type="number"]');
    rng.min = min; rng.max = max; rng.step = step;
    num.min = min; num.max = max;
  });
}

/* 播放控制 */
function setPlaying(on) {
  playing = on;
  $('btn-play').textContent = on ? '⏸ 暂停' : '▶ 播放';
  $('btn-play').classList.toggle('active', on);
  /* 迷你条播放键同步（面板隐藏时与主播放键同态） */
  const bm = $('btn-play-mini');
  if (bm) { bm.textContent = on ? '⏸ 暂停' : '▶ 播放'; bm.classList.toggle('active', on); }
}
$('btn-play').addEventListener('click', () => setPlaying(!playing));
$('btn-play-mini').addEventListener('click', () => setPlaying(!playing));

/* ── 展开/隐藏控制面板（全局视角） ───────────────────────────────
   隐藏：body.panel-hidden 隐藏左/右面板、顶部说明、UTC 时钟；
   把 #ph-prog 进度条移入左下角 #mini-ctrl（scrub/--prog 监听随元素走，
   无需改动）；恢复时移回 #phase-box 原位。 */
{
  const prog = $('ph-prog');
  const phaseBox = $('phase-box');
  const miniCtrl = $('mini-ctrl');
  const toggleBtn = $('btn-panel-toggle');
  /* 记录 #ph-prog 在 #phase-box 中的原始后继节点，便于精确还原位置 */
  const progNext = prog ? prog.nextSibling : null;
  function setPanelHidden(hidden) {
    document.body.classList.toggle('panel-hidden', hidden);
    toggleBtn.textContent = hidden ? '⤢ 展开面板' : '⤢ 隐藏面板';
    toggleBtn.classList.toggle('active', hidden);
    if (prog && miniCtrl && phaseBox) {
      if (hidden) {
        miniCtrl.appendChild(prog);          /* 移入迷你条（flex 占满剩余宽度） */
      } else {
        phaseBox.insertBefore(prog, progNext); /* 还原到原位置 */
      }
    }
  }
  toggleBtn.addEventListener('click', () => {
    setPanelHidden(!document.body.classList.contains('panel-hidden'));
  });
}
$('sel-speed').addEventListener('change', (ev) => { speedX = +ev.target.value; });
/* 仿真时钟重新对齐到当前实际 UTC（保持当前轨道相位） */
$('btn-utc-sync').addEventListener('click', () => { syncEpochToNow(); updateUtcClock(); });
$('btn-peri').addEventListener('click', () => {
  els.nu = 0; tSincePeri = 0; syncInputs('nu');
});
$('btn-fit').addEventListener('click', fitView);
function setFpv(on) {
  fpv = on;
  $('btn-fpv').classList.toggle('active', fpv);
  $('btn-ov').classList.toggle('active', !fpv);
  applyDisplay();                 /* 第一视角下隐藏贴身元素（标签/本体/矢量/水平盘/位置线） */
  if (fpv) {
    fpvCam.az = 0; fpvCam.el = 0;
    setPlaying(true);             /* 进入第一视角即开始播放 */
  } else {
    camera.fov = 50; camera.updateProjectionMatrix();   /* 还原概览视场 */
  }
}
$('btn-fpv').addEventListener('click', () => setFpv(true));
$('btn-ov').addEventListener('click', () => setFpv(false));

/* L2 会合系视角开关（仅 CR3BP）：惯性/会合系轨道线互换，相机对准轨迹包围球 */
function setL2View(on, skipFit) {
  l2View = on && !!(prop && prop.isCr3bp);
  if (l2View && topView) setTopView(false, true);   /* 与俯视·地球静止视角互斥（共用同一自转状态源） */
  if (!l2View) l2SpinFrozen = null;   /* 退出 L2 视角时清理拖动冻结态（共转随 l2View 关闭自动失效） */
  $('btn-l2').classList.toggle('active', l2View);
  $('hint-l2').classList.toggle('show', l2View);
  if (orbitLine) orbitLine.visible = !l2View;
  if (orbitLineSyn) orbitLineSyn.visible = l2View;
  if (trailObj) trailObj.visible = !l2View;
  applyDisplay();   /* 信标 / 地球方位标签的可见性随视角与天体变化 */
  if (!skipFit) fitView();
}
$('btn-l2').addEventListener('click', () => setL2View(!l2View));

/* 俯视·地球居中视角（全部月球轨道）：地心惯性俯视——applyCamera 每帧把相机
   target 压到实时地球位置，地球静止于画面中心，月球与航天器环绕地球运动；
   不随会合系共转（不再复用 L2 共转机构），所见即 CR3BP 轨迹环绕地球的运动。
   与 L2 视角互斥；仅相机行为 + 地心轨迹线/月球轨道环显隐，不互换轨迹线。
   el=1.45 = 拖动弹回上限（≈83° 俯角，避免 el=π/2 时 up∥视线导致 lookAt 退化）。 */
function setTopView(on, skipFit) {
  topView = on && body === 'moon';
  if (topView && l2View) setL2View(false, true);
  if (!topView) {
    l2SpinFrozen = null;   /* 退出时清理拖动冻结态 */
    /* 退出俯视：蓝线（月心惯性轨迹）与渐隐拖尾恢复常规可见性——
       俯视下可能被勾选框隐藏，退出后不得永久丢失（同蓝线 2026-10-15 教训） */
    if (orbitLine) orbitLine.visible = !l2View;
    if (trailObj) trailObj.visible = !l2View;
  }
  $('btn-top').classList.toggle('active', topView);
  $('hint-top').classList.toggle('show', topView);
  if (topView) { camCtl.el = 1.45; camCtl.az = 0; }
  applyDisplay();   /* geoLine / moonRing 可见性随视角切换 */
  if (!skipFit) fitView();
}
$('btn-top').addEventListener('click', () => setTopView(!topView));

/* 预设 */
document.querySelectorAll('[data-preset]').forEach((btn) => {
  btn.addEventListener('click', () => {
    els = { ...PRESETS[btn.dataset.preset] };
    cr3bpState.active = false;
    moonPresetKey = null;
    syncAll();
    rebuild();
    fitView();
  });
});

/* ── 中心天体选项卡切换 ───────────────────────────────────── */
function setBody(b) {
  if (b === body) return;
  body = b;
  applyBodyRanges();
  $('tab-earth').classList.toggle('active', b === 'earth');
  $('tab-moon').classList.toggle('active', b === 'moon');
  /* 传播模型下拉随天体切换：地球前两个 / 月球后三个 */
  syncPropOptions();
  /* 面板显隐 */
  $('presets-earth').style.display = b === 'earth' ? '' : 'none';
  $('h-presets').style.display = b === 'earth' ? '' : 'none';
  $('presets-moon').style.display = b === 'moon' ? '' : 'none';
  $('h-moon-presets').style.display = b === 'moon' ? '' : 'none';
  $('warn-atmo').classList.toggle('show', false);
  updateBodyLabels();
  /* 场景对象切换 */
  earth.visible = b === 'earth';
  glow.visible = b === 'earth';
  moonMesh.visible = b === 'moon';
  earthMini.visible = false;   /* 由 updateScene 每帧控制 */
  lpointGroup.visible = false;
  if (b === 'earth') {
    if (!moonOrbitRing) { moonOrbitRing = buildMoonOrbitRing(); scene.add(moonOrbitRing); }
    moonOrbitRing.visible = true;
    moonMesh.visible = true;          /* 地球选项卡：月球迷你天体 */
    /* 恢复地球默认预设 */
    els = { ...PRESETS.molniya };
    cr3bpState.active = false;
    moonPresetKey = null;
    $('sel-prop').value = 'kepler';
    syncPropUI();
    syncAll();
    rebuild();
    fitView();
    setElSlidersEnabled(true);
    updateStatusLine();
  } else {
    if (moonOrbitRing) moonOrbitRing.visible = false;
    /* 月球默认预设：环月圆轨道 CLO（Halo 等三体预设手动点选） */
    applyMoonPreset('clo');
    updateStatusLine();
  }
}
$('tab-earth').addEventListener('click', () => setBody('earth'));
$('tab-moon').addEventListener('click', () => setBody('moon'));

/* 中心天体相关文案：近地点/远地点 ↔ 近月点/远月点（左侧控件 + 按钮） */
function updateBodyLabels() {
  const m = body === 'moon';
  const peri = m ? '近月点' : '近地点', apo = m ? '远月点' : '远地点';
  $('el-name-w').innerHTML = `${peri}幅角 <b>ω</b>`;
  $('d-hp-label').innerHTML = `${peri}高度 h<sub>p</sub>`;
  $('d-ha-label').innerHTML = `${apo}高度 h<sub>a</sub>`;
  $('btn-peri').innerHTML = `⟲ 回${peri}`;
  $('btn-peri').title = `相位回到 ν=0（过${peri}）`;
}

function updateStatusLine() {
  const p = $('status-body');
  if (body === 'earth') {
    p.innerHTML = '拖动左侧根数实时求解位置 / 高度 / 速度 / 飞行路径角 · ECI 惯性系（地球随仿真时间自转）'
      + ' · μ = 398600.4418 km³/s² · R⊕ = 6378.137 km'
      + ' · 摄动模式：J2–J4 + 大气阻力 + 太阳光压 数值积分（90 d 窗口）<span id="status-note"></span>';
  } else {
    p.innerHTML = '月心坐标系（ẑ = 月球北极，SPICE MOON_PA 历元锁定 2026-10-13）'
      + ' · μ☾ = 4902.800 km³/s² · R☾ = 1737.4 km'
      + ' · 月心摄动：J2–J4（GL0660B）+ 地球三体 + 光压 · 三体模式：地月会合系 CR3BP'
      + '<span id="status-note"></span>';
  }
  /* 顶栏错误信息重挂（status-note 重建后保留） */
  if (window.__errs && window.__errs.length) {
    const note = $('status-note');
    if (note) note.textContent = ' · 脚本错误：' + window.__errs.slice(-2).join(' | ');
  }
}

/* ── 月球预设 ──────────────────────────────────────────────── */
function applyMoonPreset(key) {
  moonPresetKey = key;
  document.querySelectorAll('[data-mpreset]').forEach((b) =>
    b.classList.toggle('active', b.dataset.mpreset === key));
  cr3bpState.active = ['halo', 'nrho', 'dro', 'dro31'].includes(key);
  /* DRO 按钮内含共振比下拉（dro/dro31）：切来时保留上次所选共振比，
     从非 DRO 预设切来则回落 2:1 */
  if (key === 'dro') {
    if (!['dro', 'dro31'].includes(cr3bpState.preset)) { cr3bpState.preset = 'dro'; }
  } else {
    cr3bpState.preset = cr3bpState.active ? key : null;
  }
  /* 轨道修正仅鹊桥一号 Halo 提供：切走即关闭并隐藏勾选框 */
  $('row-stationkeep').style.display = key === 'halo' ? '' : 'none';
  if (key !== 'halo') { $('chk-stationkeep').checked = false; cr3bpState.stationKeeping = false; }
  const desc = {
    clo: '环月圆轨道：12 h 周期，倾角可调。演示月球非球形与地球三体摄动下的长期演化（偏心率年增长 ~0.001 量级）。',
    ell: '环月椭圆轨道：近月点 15 km × 远月点 200 km（a ≈ 1845 km，e ≈ 0.05），类 LRO 低轨侦查/着陆转移停泊轨道。低轨受月球非球形引力场（GL0660B）与地球三体显著影响，默认月心数值摄动传播。',
    elfo: '鹊桥二号使命轨道：24 h 周期太阳同步冻结测控回归轨道。冻结条件 e=√(1−5cos²i/3)，ω 冻结在 90°，拱线长期稳定（在轨实际周期 24.84 h、e=0.7824、i=119.2°）。默认地月三体（CR3BP）传播——冻结设计核心即地球三体引力。',
    halo: '鹊桥一号使命轨道：地月 L2 南族晕（Halo）轨道，Az≈13000 km、T≈14 d，为嫦娥四号月背中继运行 5 年以上。三体会合系传播。',
    nrho: 'Gateway 类近直线晕轨道（NRHO）：L2 南族 9:2 共振，T=6.56 d、近月 ~1500 km，兼顾月背覆盖与地球通信。三体会合系传播。',
    dro: '远距离逆行轨道（DRO）：地月会合系下绕月逆行周期轨道族，按共振比 n:1（月球公转 n 圈期间航天器绕月 n 圈）分档，长期近似稳定。无 1:1 成员——1:1 意味着会合系中静止，即平衡点（L4/L5）而非轨道；族延拓在 T→恒星月 处终止（数值打靶验证）。三体会合系传播。',
    dro31: 'DRO 3:1：共振比 3（T=恒星月/3≈9.11 d），月心距 4.9万–5.8万 km，DRO 族延拓打靶解。近月小环本征不稳定，按理想周期参考轨道折叠显示（每圈精确闭合）。三体会合系传播。',
  }[key === 'dro' ? cr3bpState.preset : key] || '—';
  $('mp-desc').textContent = desc;
  buildMoonPresetParams(key);
  if (!cr3bpState.active) {
    /* 根数型预设：CLO 默认二体（看纯开普勒椭圆）；ELL（15×200 km 低轨）默认月心摄动；
       ELFO 默认三体模式 —— S13 冻结设计核心即地球三体（1/2 阶）+月球 J2，CR3BP 直接传播其根数初值；
       perturbed-moon 仍可在下拉手动切换（mp-grav/mp-earth3b/mp-srp 仅对其生效） */
    els = { ...MOON_PRESETS[key] };
    $('sel-prop').value = key === 'elfo' ? 'cr3bp' : key === 'ell' ? 'perturbed-moon' : 'kepler-moon';
    syncPropUI();
    syncAll();
    rebuild();
    fitView();
  } else {
    /* CR3BP 型预设 */
    els = { a: 10000, e: 0.6, i: 90, Om: 0, w: 90, nu: 0 };   /* 占位（根数面板隐藏） */
    $('sel-prop').value = 'cr3bp';
    syncPropUI();
    rebuild();
    fitView();
  }
  setElSlidersEnabled(!cr3bpState.active);
}
document.querySelectorAll('[data-mpreset]').forEach((btn) => {
  btn.addEventListener('click', () => applyMoonPreset(btn.dataset.mpreset));
});

/* 月球预设可调参数（每个轨道独立） */
function buildMoonPresetParams(key) {
  const box = $('mp-params');
  box.innerHTML = '';
  $('mp-frozen').style.display = 'none';
  const mkRow = (label, min, max, step, val, unit, onInput) => {
    const row = document.createElement('div');
    row.className = 'param-row';
    row.innerHTML = `<label style="min-width:88px">${label}</label>`
      + `<input type="range" min="${min}" max="${max}" step="${step}" value="${val}">`
      + `<span class="mp-val">${(+val).toFixed(1)} ${unit}</span>`;
    const rng = row.querySelector('input');
    rng.addEventListener('input', () => {
      row.querySelector('.mp-val').textContent = `${(+rng.value).toFixed(1)} ${unit}`;
      onInput(+rng.value);
    });
    box.appendChild(row);
    return rng;
  };
  if (key === 'clo') {
    /* CLO 不提供额外滑块：周期 T ↔ 半长轴 a、倾角 i 与左侧根数面板完全重复 */
  } else if (key === 'elfo') {
    /* ELFO 同 CLO：周期 T ↔ 半长轴 a、倾角 i 与左侧根数面板完全重复，不提供额外滑块；
       仅显示冻结关系说明（B10 式(3)，由预设根数 i/a 计算） */
    const i = els.i * C.D2R;
    const eFroz = Math.sqrt(Math.max(0, 1 - 5 * Math.cos(i) ** 2 / 3));
    const frozen = $('mp-frozen');
    frozen.style.display = '';
    frozen.innerHTML = `冻结关系（B10 式(3)）：e = √(1 − 5cos²i/3) → i = ${els.i.toFixed(2)}° ⇒ e = ${eFroz.toFixed(4)}`
      + ` · 近月点高度 h<sub>p</sub> = ${(els.a * (1 - eFroz) - X3.R_MOON).toFixed(0)} km（S13 鹊桥二号使命轨道：逆行 i≈119.2°，ω 冻结 90°）`;
  } else if (key === 'halo') {
    mkRow('Az 振幅', 5000, 40000, 500, cr3bpState.azKm, 'km', (v) => {
      cr3bpState.azKm = v; rebuild();
    });
    mkRow('起始相位', 0, 1, 0.01, 0, '', (v) => {
      cr3bpState.phase = v; rebuild();
    });
  } else if (key === 'dro') {
    /* DRO 共振比下拉（2:1 / 3:1）：切换 cr3bpState.preset 并重建。
       1:1 不提供——会合系静止 = 平衡点，非轨道（族延拓 T→2π 处断链验证）；
       4:1 已按需求移除（2026-10-15） */
    const row = document.createElement('div');
    row.className = 'param-row';
    row.innerHTML = `<label style="min-width:88px">共振比</label>`
      + `<select class="ctl" style="flex:1">`
      + `<option value="dro">2:1（T=13.66 d，稳定）</option>`
      + `<option value="dro31">3:1（T=9.11 d）</option>`
      + `</select>`;
    const selDro = row.querySelector('select');
    selDro.value = ['dro', 'dro31'].includes(cr3bpState.preset) ? cr3bpState.preset : 'dro';
    selDro.addEventListener('change', () => {
      cr3bpState.preset = selDro.value;
      /* 共振比切换 → 更新描述（对应 desc 子条目） */
      const d = {
        dro: 'DRO 2:1：共振比 2（T=恒星月/2≈13.66 d），月心距 6.9万–9.3万 km，DRO 族最外稳定成员，长期近似稳定无需轨控。三体会合系传播。',
        dro31: 'DRO 3:1：共振比 3（T=恒星月/3≈9.11 d），月心距 4.9万–5.8万 km，DRO 族延拓打靶解。近月小环本征不稳定，按理想周期参考轨道折叠显示（每圈精确闭合）。三体会合系传播。',
      }[cr3bpState.preset];
      if (d) { $('mp-desc').textContent = d; }
      rebuild();
    });
    box.appendChild(row);
    mkRow('起始相位', 0, 1, 0.01, 0, '', (v) => {
      cr3bpState.phase = v; rebuild();
    });
  } else if (key === 'nrho') {
    mkRow('起始相位', 0, 1, 0.01, 0, '', (v) => {
      cr3bpState.phase = v; rebuild();
    });
  }
}

/* 显示开关 */
function applyDisplay() {
  const lb = $('chk-labels').checked;
  const nodeOn = $('chk-node').checked;
  const fpvHide = fpv;   /* 第一视角：隐藏贴身元素，只留地球/轨道/弧线等远景信息 */
  velArrow.visible = rArrow.visible = $('chk-vel').checked && !fpvHide;
  /* CR3BP 受力分析：仅三体模式 + 勾选 + 非第一视角；其标签随总标签开关 */
  forceGroup.visible = $('chk-forces').checked && !!(prop && prop.isCr3bp) && !fpvHide;
  for (const k in fLabels) fLabels[k].visible = lb;
  disc.visible = $('chk-disc').checked && !fpvHide;
  scGroup.visible = posLine.visible = !fpvHide;
  if (nodeLine) nodeLine.visible = nodeOn;
  if (wArc) wArc.visible = nodeOn;     /* ω 弧随「升交点节线」开关 */
  /* ν 弧始终显示（标注当前轨道相位），仅其文字受标签开关控制 */
  const _gridOn = $('chk-grid').checked;
  if (grid) grid.visible = _gridOn;
  axesGroup.visible = _gridOn;
  /* 地球方位信标仅在月球选项卡显示（地球选项卡下地球本体就在视野中心） */
  /* 2026-10-14：地球信标点已移除，标签直接锚在真实地球球体上（updateScene 每帧） */
  for (const [k, sp] of Object.entries(labels)) {
    let on = lb && !fpvHide && ((k !== 'node' && k !== 'omega') || nodeOn);
    if (k === 'earth') on = on && body === 'moon';
    /* 拱点标签：当前模式无拱点标记时（周期轨道族 halo/nrho/dro*）恒隐藏。
       updateCr3bpReveal/rebuild 只在有拱点时登记 anchor，若无此门控，
       标签开关会把它们点亮在未定位的默认原点（月心）——即用户所见 bug */
    if ((k === 'peri' || k === 'apo') && !apsisAvailable) on = false;
    sp.visible = on;
  }
  /* L1/L2 标记的文字标签不在 labels 注册表（lpointGroup 子对象），随同一开关隐藏 */
  for (const mk of lpointGroup.children) {
    if (mk.name.endsWith('-label')) mk.visible = lb && !fpvHide;
  }
  /* 俯视·地球居中：月球地心轨道环（浅蓝）恒随该视角显示。
     CR3BP 下蓝线（相对月球的月心惯性轨迹）与绿线（相对地球的地心惯性轨迹）
     各有独立勾选——只切 visible 不触发重建；勾选行仅在此视角 + CR3BP 出现。 */
  const topCr3 = topView && !!(prop && prop.isCr3bp);
  $('row-top-trails').style.display = topCr3 ? '' : 'none';
  if (moonRing) moonRing.visible = topView;
  if (geoLine) geoLine.visible = topCr3 && $('chk-trail-earth').checked;
  if (orbitLine && topCr3) orbitLine.visible = $('chk-trail-moon').checked;
  /* 渐隐拖尾与蓝线同一勾选管辖（暂停时帧闸门不推进也要即时生效） */
  if (trailObj && topCr3) trailObj.visible = $('chk-trail-moon').checked;
}
['chk-vel', 'chk-disc', 'chk-node', 'chk-grid', 'chk-labels', 'chk-forces',
 'chk-trail-moon', 'chk-trail-earth'].forEach((id) =>
  $(id).addEventListener('change', applyDisplay));
/* 一键开关所有文字标签：与「文字标签」复选框同一状态源，按钮高亮双向同步 */
function syncLabelsBtn() { $('btn-labels').classList.toggle('active', $('chk-labels').checked); }
$('btn-labels').addEventListener('click', () => {
  $('chk-labels').checked = !$('chk-labels').checked;
  syncLabelsBtn();
  applyDisplay();
});
$('chk-labels').addEventListener('change', syncLabelsBtn);

/* 推导展开开关 */
document.querySelectorAll('[data-deriv]').forEach((btn) => {
  btn.addEventListener('click', () => {
    const panel = $(btn.dataset.deriv + '-deriv');
    const on = !panel.classList.contains('show');
    panel.classList.toggle('show', on);
    btn.classList.toggle('on', on);
    btn.textContent = on ? '推导 ∧' : '推导 ∨';
  });
});

/* 符号速查折叠开关 */
$('sym-toggle').addEventListener('click', () => {
  const on = $('sym-legend').classList.toggle('collapsed');
  $('sym-toggle').textContent = on ? '展开 ∨' : '收起 ∧';
});

/* 传播模型下拉（注册表驱动，按中心天体过滤：地球前两个 / 月球后三个，不交叠显示） */
const PROP_BY_BODY = { earth: ['kepler', 'perturbed'], moon: ['kepler-moon', 'perturbed-moon', 'cr3bp'] };
function syncPropOptions() {
  const sel = $('sel-prop');
  const ids = PROP_BY_BODY[body];
  sel.innerHTML = '';
  for (const id of ids) {
    const o = document.createElement('option');
    o.value = id; o.textContent = C.propagators[id].label;
    sel.appendChild(o);
  }
  if (!ids.includes(sel.value)) sel.value = ids[0];
}
function syncPropUI() {
  const sel = $('sel-prop');
  const cls = C.propagators[sel.value];
  const isP = cls && cls.hasPerturbation;
  const isEarthPert = isP && sel.value === 'perturbed';
  const isMoonPert = isP && sel.value === 'perturbed-moon';
  const isCr3bp = sel.value === 'cr3bp';
  $('pert-ctl').style.display = isEarthPert ? '' : 'none';
  /* mp-grav/mp-earth3b/mp-srp 仅 perturbed-moon 生效：CR3BP 隐含地球三体，面板对其隐藏以免误导 */
  $('pert-ctl-moon').style.display = isMoonPert ? '' : 'none';
  $('prop-note-kepler').style.display = isP ? 'none' : '';
  if (cls) {
    $('prop-note').textContent =
      isCr3bp ? (cr3bpState.active
        ? '当前：地月会合系三体传播（90 d 窗口，周期轨道已预解）。'
        : '当前：地月会合系三体传播（90 d 窗口，根数初值直接积分）。')
      : isMoonPert ? '当前：月心数值摄动（90 d 窗口）。'
      : isP ? '当前：数值摄动（90 d 窗口）。'
      : '当前：开普勒二体（Δr ≡ 0，Δv ≡ 0）。';
  }
}
{
  const sel = $('sel-prop');
  syncPropOptions();
  sel.value = 'kepler';
  $('prop-note').textContent = '当前：开普勒二体（Δr ≡ 0，Δv ≡ 0）。';
  syncPropUI();
  sel.addEventListener('change', () => {
    const cls = C.propagators[sel.value];
    if (!cls) { sel.value = 'kepler'; return; }
    /* CR3BP 只在月球选项卡可用 */
    if (sel.value === 'cr3bp' && body !== 'moon') {
      $('status-note').textContent = '（三体会合系仅月球选项卡可用，已切回）';
      sel.value = 'kepler';
      syncPropUI();
      return;
    }
    syncPropUI();
    rebuild();
  });
  ['pp-grav', 'pp-drag', 'pp-srp', 'pp-moon3b'].forEach((id) =>
    $(id).addEventListener('change', () => { if (C.propagators[sel.value] && C.propagators[sel.value].hasPerturbation) rebuild(); }));
  ['mp-grav', 'mp-earth3b', 'mp-srp'].forEach((id) =>
    $(id).addEventListener('change', () => { if (sel.value === 'perturbed-moon') rebuild(); }));
  ['pp-cd', 'pp-cr', 'pp-sm', 'pp-lam'].forEach((id) =>
    $(id).addEventListener('change', () => { if (C.propagators[sel.value] && C.propagators[sel.value].hasPerturbation) rebuild(); }));

  /* Halo 轨道修正勾选：切换后重建传播器，toast 提示 3 秒 */
  $('chk-stationkeep').addEventListener('change', () => {
    cr3bpState.stationKeeping = $('chk-stationkeep').checked;
    const toast = $('sk-toast');
    const enabled = cr3bpState.stationKeeping;
    toast.classList.toggle('sk-off', !enabled);
    toast.querySelector('.sk-tx').textContent = enabled ? '轨道修正已启用' : '轨道修正已关闭';
    toast.querySelector('.sk-sub').textContent = enabled
      ? `每 T/8 ≈44 h 单脉冲 · 90 d 累计 Δv ~1 mm/s 量级`
      : '无修正时 Halo 约 3–4 圈（≈1200 h）后漂离（L2 本征指数不稳定）';
    toast.hidden = false;
    clearTimeout(window.__skToastT);
    window.__skToastT = setTimeout(() => { toast.hidden = true; }, 3000);
    rebuild();
  });
}

/* 进度条拖动 scrub：拖到任意位置 → tNow = frac×window，
   拖尾 / 进度条 / 右侧图表指针 / 卫星位置全部同步。拖动期间暂停播放。
   摄动与 CR3BP 共用（同一 90 d 窗口）；CR3BP 主时钟 = tSincePeri，暂停后
   每帧仅按它锚定 tNow，不会再被 ν 折叠回退。 */
{
  const bar = $('ph-prog');
  let dragging = false;
  const fracAt = (ev) => {
    const r = bar.getBoundingClientRect();
    return r.width ? Math.min(1, Math.max(0, (ev.clientX - r.left) / r.width)) : 0;
  };
  const scrub = (ev) => {
    if (!(prop && prop.isPerturbed)) return;
    setPlaying(false);
    tSincePeri = fracAt(ev) * prop.windowSeconds();
    const nuRad = prop.trueAnomalyAtTime(tSincePeri);   /* 锚定 prop.tNow / _osc */
    /* CR3BP 相位代理未折叠（可超 2π）；ν 显示折叠到 [0,2π)，主时钟不动 */
    els.nu = (prop.isCr3bp ? C.wrap2pi(nuRad) : nuRad) * C.R2D;
    syncInputs('nu');
    lastChartTNow = -1;                                   /* 强制下一帧重绘图表指针 */
    updateProgBar();                                      /* 立即刷新，跟手 */
    updateTrailReveal();
    if (prop.isCr3bp) updateCr3bpReveal();                /* 三体尾迹同步收缩/延伸 */
  };
  bar.addEventListener('pointerdown', (ev) => {
    if (!(prop && prop.isPerturbed)) return;
    dragging = true;
    try { bar.setPointerCapture(ev.pointerId); } catch (_) {}
    scrub(ev);
    ev.preventDefault();
  });
  bar.addEventListener('pointermove', (ev) => { if (dragging) scrub(ev); });
  const release = (ev) => {
    if (!dragging) return;
    dragging = false;
    try { bar.releasePointerCapture(ev.pointerId); } catch (_) {}
  };
  bar.addEventListener('pointerup', release);
  bar.addEventListener('pointercancel', release);
}

/* ── 页内自检（与 node_check_elements.cjs 同判据） ───────────────── */
function runSelfTest() {
  const badge = $('selftest');
  const detail = $('selftest-detail');
  const D = window.ELEMENTS_SELFTEST;
  if (!D) {
    badge.textContent = '✗ 基准数据缺失（static/elements_selftest.js 未加载）';
    badge.className = 'bad';
    return;
  }
  let pass = 0, fail = 0;
  const lines = [];
  const vecErr = (v, ref) => Math.max(
    Math.abs(v[0] - ref[0]), Math.abs(v[1] - ref[1]), Math.abs(v[2] - ref[2])
  ) / Math.hypot(ref[0], ref[1], ref[2]);

  for (const cs of D.cases) {
    const mu = cs.mu || C.MU;               // 月球用例带自身 μ☾（DE440）
    const bodyR = cs.body_r || C.R_EARTH;   // 月球用例带 R☾
    const el = { a: cs.a, e: cs.e, i: cs.i * C.D2R, Om: cs.Om * C.D2R, w: cs.w * C.D2R };
    const nu = cs.nu * C.D2R;
    const st = C.classicalToState(el, nu, mu, bodyR);
    const checks = [
      ['r_ECI', vecErr(st.rE, cs.r_eci_km) < 1e-9],
      ['v_ECI', vecErr(st.vE, cs.v_eci_kms) < 1e-9],
      ['r', Math.abs(st.r - cs.r_km) / cs.r_km < 1e-9],
      ['v', Math.abs(st.v - cs.v_kms) / cs.v_kms < 1e-9],
      ['alt', Math.abs(st.alt - cs.alt_km) < 1e-6],
      ['γ', C.angDist(st.gamma, cs.gamma_deg * C.D2R) < 1e-8],
      ['E/M', C.angDist(st.E, cs.E_deg * C.D2R) < 1e-8 && C.angDist(st.M, cs.M_deg * C.D2R) < 1e-8],
      ['T', Math.abs(C.period(cs.a, mu) - cs.period_s) / cs.period_s < 1e-9],
      ['链路闭环', C.angDist(C.eToNu(C.mToE(C.nuToM(nu, cs.e).M, cs.e).E, cs.e), nu) < 1e-10],
    ];
    const bad = checks.filter((c) => !c[1]).map((c) => c[0]);
    if (bad.length === 0) { pass += checks.length; lines.push(`✓ ${cs.name}`); }
    else { fail += bad.length; pass += checks.length - bad.length; lines.push(`✗ ${cs.name}：${bad.join(', ')} 不符`); }
  }
  const total = pass + fail;
  if (fail === 0) {
    badge.textContent = `✓ 自检 ${pass}/${total} · 基准 poliastro ${D.meta.poliastro_version}（点击查看明细）`;
    badge.className = 'ok';
  } else {
    badge.textContent = `✗ 自检 ${pass}/${total} 有失败项（点击查看明细）`;
    badge.className = 'bad';
  }
  lines.push('', `基准: ${D.meta.benchmark}`,
    `poliastro ${D.meta.poliastro_version} · astropy ${D.meta.astropy_version}`,
    `生成时双路最大相对误差: r=${D.meta.max_rel_err_r} · v=${D.meta.max_rel_err_v} · T=${D.meta.max_rel_err_T}`);
  detail.textContent = lines.join('\n');
  badge.addEventListener('click', () => detail.classList.toggle('show'));
  detail.addEventListener('click', () => detail.classList.remove('show'));
}

/* ── 主循环 ──────────────────────────────────────────────────────── */
let lastFrame = performance.now();
let lastText = 0;
let lastActivity = performance.now();
let _frameIdx = 0;
const RENDER_EVERY = SOFT_GL ? 4 : 1;   /* 软件渲染：~15fps 上限，避免吃满 CPU */
/* 任意用户交互刷新活跃时间戳（用于空闲停帧） */
['pointerdown', 'pointermove', 'pointerup', 'wheel', 'keydown', 'input', 'change'].forEach((t) =>
  window.addEventListener(t, () => { lastActivity = performance.now(); }, { capture: true, passive: true }));
window.addEventListener('visibilitychange', () => { if (!document.hidden) { lastFrame = performance.now(); lastActivity = lastFrame; } });

function tick(now) {
  if (document.hidden) { requestAnimationFrame(tick); return; }   /* 标签隐藏：不渲染 */
  const dt = Math.min(0.1, (now - lastFrame) / 1000);
  lastFrame = now;

  let nuRad;
  if (playing) {
    tSincePeri += dt * speedX;
    if (prop && prop.isPerturbed) {
      // 摄动窗口播放：钳制到 [0, window]，播完自动停
      if (tSincePeri >= prop.windowSeconds()) {
        tSincePeri = prop.windowSeconds();
        setPlaying(false);
      }
      nuRad = prop.trueAnomalyAtTime(tSincePeri);
      /* CR3BP 的 nuRad 是未折叠相位代理（t/T·2π，随圈数递增）；
         els.nu 只作相位显示，折叠到 [0,2π)，主时钟仍是 tSincePeri */
      els.nu = (prop.isCr3bp ? C.wrap2pi(nuRad) : nuRad) * C.R2D;
    } else {
      nuRad = prop.trueAnomalyAtTime(tSincePeri);
      els.nu = nuRad * C.R2D;
    }
  } else {
    nuRad = els.nu * C.D2R;
    if (prop && prop.isPerturbed) {
      /* Cowell/CR3BP 的 stateAtTrueAnomaly 忽略 ν、只读 prop.tNow；暂停时
         必须把 tSincePeri（播放 / 进度条拖动 / ν 滑条都写它，是唯一主时钟）
         回灌为 prop.tNow，否则拖 ν 滑条航天器冻结。绝不能由 ν 反推
         tSincePeri：ν 折叠到 [0,2π)，反推会丢掉整周期数——播放越过一圈
         后暂停，进度回退到起点。 */
      prop.trueAnomalyAtTime(tSincePeri);
    }
  }
  const st = prop.stateAtTrueAnomaly(nuRad);
  /* 地球随仿真时间绕极轴（场景 Y，即 ECI z）做恒星自转，方向与公转/真近点角
     推进一致：地球在 ECI 中绕 Y 轴东向自转，angle = ω⊕·t。仅在播放或拖动 nu
     改变 tSincePeri 时随时间推进，暂停静止时停转。月球选项卡下地球迷你球
     不自转（视觉简化），月球本体自转在 updateScene 中处理（潮汐锁定）。 */
  if (body === 'earth') earth.rotation.y = EARTH_OMEGA * tSincePeri;
  updateScene(st);
  /* 摄动拖尾：每帧按 tNow 滑动 7 d 窗口（drawRange[start,count]），
     不画整条 90 d 轨迹，避免 LEO 千余圈糊成密带看不清 */
  if (prop && prop.isPerturbed) updateTrailReveal();
  /* CR3BP：轨道线由航天器随 tSincePeri 逐段画出（不预渲染整周期） */
  if (prop && prop.isCr3bp) updateCr3bpReveal();
  /* 每帧把相机状态（拖拽旋转 / 滚轮与捏合缩放 / 适应视角按钮）写入相机，
     否则 pointermove 只改了 camCtl 数值，画面永远不会动 */
  applyCamera();
  if (now - lastText > 80) { updateTexts(st, nuRad); updateUtcClock(); lastText = now; }

  /* 暂停且无交互超 5s → 跳过渲染（场景静止，省电；软件渲染尤其受益） */
  const idle = !playing && (now - lastActivity > 5000);
  if (!idle && (_frameIdx++ % RENDER_EVERY) === 0) {
    try {
      renderer.render(scene, camera);
    } catch (err) {
      /* 抛错几乎都是贴图被 WebGL 拒绝（典型：file:// 直开，本地图 tainted）。
         地球与月球都要回退，否则一个抛错每帧复发、另一个永远黑球。 */
      if (!texFailed) useProceduralEarth('file:// 跨源拦截');
      if (!moonTexFailed) useProceduralMoon('file:// 跨源拦截（浏览器禁止本地贴图进 WebGL，请通过 http 访问）');
      renderer.render(scene, camera);
    }
  }
  requestAnimationFrame(tick);
}

window.addEventListener('resize', () => {
  camera.aspect = window.innerWidth / window.innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(window.innerWidth, window.innerHeight);
  THREE.setFatLinePixelRatio(renderer.getPixelRatio());
  THREE.setFatLineResolution(renderer.domElement.width, renderer.domElement.height);
  drawCharts();   // 画布分辨率随窗口重绘
});

/* ── 启动 ────────────────────────────────────────────────────────── */
syncAll();
rebuild();
fitView();
applyDisplay();
applyCamera();
runSelfTest();
requestAnimationFrame(tick);

/* 自诊断：?selftest=1 时自动派发 wheel + 模拟 pointer 拖拽，把相机变化写进 DOM，
 * 便于用 headless --dump-dom 一次性读取（判断事件是否到达 canvas、相机是否更新）。 */
if (new URLSearchParams(location.search).has('selftest')) {
  setTimeout(() => {
    const note = $('status-note');
    const canvas = document.querySelector('#stage canvas');
    if (!canvas) { note.textContent = ' [selftest] 无 canvas'; return; }
    /* 暴露相机对象供外部（诊断/测试）读取 */
    window.__camCtl = camCtl;
    const r = canvas.getBoundingClientRect();
    const cx = r.left + r.width / 2, cy = r.top + r.height / 2;
    const before = { az: +camCtl.az.toFixed(4), el: +camCtl.el.toFixed(4), dist: +camCtl.dist.toFixed(3) };
    canvas.dispatchEvent(new WheelEvent('wheel', { deltaY: 240, bubbles: true, cancelable: true, clientX: cx, clientY: cy }));
    const distAfterWheel = +camCtl.dist.toFixed(3);
    const mk = (type, x, y) => new PointerEvent(type, {
      pointerId: 7, bubbles: true, cancelable: true, isPrimary: true,
      clientX: x, clientY: y, button: 0, buttons: 1,
    });
    canvas.dispatchEvent(mk('pointerdown', cx, cy));
    canvas.dispatchEvent(mk('pointermove', cx + 60, cy + 20));
    canvas.dispatchEvent(mk('pointerup', cx + 60, cy + 20));
    const after = { az: +camCtl.az.toFixed(4), el: +camCtl.el.toFixed(4), dist: +camCtl.dist.toFixed(3) };
    note.textContent = ` [selftest] before=${JSON.stringify(before)} ` +
      `distAfterWheel=${distAfterWheel} afterDrag=${JSON.stringify(after)}` +
      ` topAtCenter=${document.elementFromPoint ? (document.elementFromPoint(cx, cy) || {}).tagName : 'n/a'}`;
  }, 700);
}
