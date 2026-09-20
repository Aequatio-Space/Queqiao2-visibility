/* orbits.js - Yang 论文 24 条地月转移轨道 XYZ 三维可视化
 *
 * 移植自 AgentMoons tools/interactive_viewer/orbits_3d.html 的内联脚本
 * （逻辑未改动，仅：外置以满足 CSP script-src 'self'；meta 状态栏改用
 * textContent 避免 innerHTML 注入；追加访问统计上报）。
 *
 * 坐标系：MEP (Moon-Earth Perifocal) 冻结于到达时刻 t_arrive
 *   X = 月地连线，由月->地；Y / Z 右手系完备（Z = 轨道面法向）
 *   原点 = 地心；单位约定：1 three.js unit = 1000 km
 * 数据：GET /api/yang_orbits_stream (SSE, meta/orbit/done)
 */

/* eslint-disable */
/* =================================================================
 * Yang 24-orbits XYZ viewer
 * -----------------------------------------------------------------
 * 坐标系：MEP (Moon-Earth Perifocal) 冻结于到达时刻 t_arrive
 *   X = 月地连线，由月→地（到达时月球在 +X 反方向即 -X）
 *   Y / Z = 右手系完备（Z = 轨道面法向）
 *   原点 = 地心
 * 单位约定：1 three.js unit = 1000 km（Mm）
 * =================================================================*/

(() => {
  'use strict';

  // ---------------------- constants --------------------------------
  const KM_PER_UNIT = 1000.0;          // 1 unit = 1000 km
  const R_EARTH_UNITS = 6378.137 / KM_PER_UNIT;
  const R_MOON_UNITS  = 1737.4    / KM_PER_UNIT;
  const DEFAULT_MOON_DIST_UNITS = 384400 / KM_PER_UNIT;  // 仅用于备份参考

  // ---------------------- scene bootstrap --------------------------
  const stage = document.getElementById('stage');
  const scene = new THREE.Scene();
  scene.background = new THREE.Color(0x0a0e1a);

  const camera = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 1, 2000);
  camera.position.set(-600, 0, 0);   // default from-moon POV: X→月，+X指向地球，-X指向月球

  const renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.setPixelRatio(window.devicePixelRatio);
  renderer.setSize(window.innerWidth, window.innerHeight);
  stage.appendChild(renderer.domElement);

  // simple in-house orbit controls (the vendored GLTFLoader doesn't ship one)
  const controls = makeOrbitControls(camera, renderer.domElement);

  // ---------------------- lights -----------------------------------
  scene.add(new THREE.AmbientLight(0xffffff, 0.55));
  const sun = new THREE.DirectionalLight(0xffffff, 0.85);
  sun.position.set(200, 200, 200);
  scene.add(sun);
  // Fill light from the opposite side so the Moon's night side isn't pitch black.
  const fill = new THREE.DirectionalLight(0x8899bb, 0.30);
  fill.position.set(-200, -150, -200);
  scene.add(fill);
  // Soft sky/ground bounce for the bodies.
  scene.add(new THREE.HemisphereLight(0x8899bb, 0x223344, 0.35));

  // ---------------------- reference bodies & frame -----------------
  const layers = {};
  layers.earth         = makeEarth();
  layers.moon          = makeMoon();
  layers.axis          = makeAxes();
  layers.moonline      = makeEarthMoonLine();
  layers.moontrack     = makeEmptyLine();
  layers.grid          = makeGrid();
  layers.earthgraticule= makeEarthGraticule();  // R-1V-e follow-up: 经纬度网格
  layers.moongraticule = makeMoonGraticule();   // R-1V-e follow-up: 经纬度网格
  layers.labels        = new THREE.Group();   // orbit # labels
  layers.endpoints     = new THREE.Group();   // start/end markers

  // Graticules must follow their parent body so they stay glued to the Moon
  // (which is pinned at the frozen arrival position in tick()).  The Earth is
  // at the origin and does not move, so the Earth graticule stays in the
  // scene root.
  layers.moon.add(layers.moongraticule);

  for (const k of Object.keys(layers)) {
    layers[k].name = 'layer_' + k;
    scene.add(layers[k]);
  }

  // 月球自转轴 + 南北极标注：挂在月球组下（随月球移动），不加入 scene 根。
  layers.moonaxis = makeMoonAxis();
  layers.moon.add(layers.moonaxis);

  // ---------------------- state ------------------------------------
  /** @type {Array<OrbitInfo>} */
  let ORBITS = [];
  let PAYLOAD = null;        // raw /api/yang_orbits response
  let SELECTED = null;        // family + orbit index
  const layerVisible = {
    earth: true, moon: true, axis: true, moonline: true,
    moontrack: true, labels: true,
    earthgraticule: true, moongraticule: true, moonaxis: true,
    grid: false, endpoints: false,
  };

  // ---------------------- load data --------------------------------
  // 进度条 DOM 引用（EventSource handler 与 finishProgress 共用）
  const progressEl    = document.getElementById('progress');
  const progressFill  = document.getElementById('progress-fill');
  const progressCur   = document.getElementById('progress-current');
  const progressTot   = document.getElementById('progress-total');
  const progressLabel = progressEl.querySelector('.label > span:first-child');

  /** 立刻把进度条切到"加载数据"状态（fetch 阶段就显示，避免空白页） */
  function showProgressStage(labelText, current, total) {
    progressLabel.textContent = labelText;
    progressCur.textContent   = String(current);
    progressTot.textContent   = String(total);
    progressFill.style.width  = total > 0 ? (current / total) * 100 + '%' : '8%';
    progressEl.classList.remove('hidden');
  }
  showProgressStage('加载轨道数据 /api/yang_orbits_stream …', 0, 0);

  // ----- Server-Sent Events loader ----------------------------------------
  // The endpoint emits one ``orbit`` event per trajectory so we can build
  // three.js geometry and advance the progress bar as each record lands,
  // instead of blocking on a single 430 KB JSON blob.
  const metaExpected = {
    t_arrive_utc: '', h_p_km: 200.0, H_m_km: 300.0,
    n_samples_per_orbit: 121, moon_model: 'circular', n_orbits_expected: 24,
  };
  PAYLOAD = { n_orbits: 0, orbits: [], ...metaExpected };
  const streamStartTs = performance.now();
  const MIN_VISIBLE_MS = 350;
  let received = 0;
  let streamClosed = false;

  function finishProgress(ok) {
    if (streamClosed) return;
    streamClosed = true;
    progressFill.style.width = '100%';
    if (!ok) {
      progressLabel.textContent = '加载失败';
      progressFill.style.background = 'var(--error)';
    }
    const wait = Math.max(0, MIN_VISIBLE_MS - (performance.now() - streamStartTs));
    setTimeout(() => progressEl.classList.add('hidden'), wait);
  }

  function appendOneOrbit(rec) {
    ORBITS.push(rec);
    if (rec && !rec.error) {
      try {
        const info = buildOneOrbit(rec);
        orbitObjects.push(info);
        addControlRow(info);
      } catch (e) {
        console.error('[orbits_3d] 构建轨道失败:', rec, e);
        return false;
      }
    }
    received += 1;
    progressCur.textContent = String(received);
    // Use n_orbits_expected (set by the meta event) for the fill ratio;
    // PAYLOAD.n_orbits is set later by the done event, so it would
    // stay 0 throughout the load and leave the bar at its initial 8%.
    const total = PAYLOAD.n_orbits_expected || 1;
    progressFill.style.width = (received / total) * 100 + '%';
    return true;
  }

  const es = new EventSource('/api/yang_orbits_stream');

  es.addEventListener('meta', ev => {
    let meta;
    try { meta = JSON.parse(ev.data); } catch (e) { meta = {}; }
    PAYLOAD = { ...PAYLOAD, ...meta };
    document.getElementById('status-msg').textContent =
      `${meta.n_orbits_expected ?? 24} 条轨道 · 每条 ${meta.n_samples_per_orbit} 个采样点 · ` +
      `t_arrive = ${meta.t_arrive_utc} · H_m = ${meta.H_m_km} km · 月模型 ${meta.moon_model}`;
    showProgressStage('正在构建轨道…', 0, meta.n_orbits_expected ?? 24);
  });

  es.addEventListener('orbit', ev => {
    let rec;
    try { rec = JSON.parse(ev.data); } catch (e) {
      console.error('[orbits_3d] 无法解析 orbit 事件', ev.data, e);
      return;
    }
    appendOneOrbit(rec);
  });

  es.addEventListener('done', ev => {
    es.close();
    PAYLOAD.n_orbits = received;
    finishProgress(true);
  });

  es.addEventListener('error', () => {
    // EventSource auto-retries by default; only flag if the connection has
    // actually dropped and no orbit has been received yet.
    if (received === 0 && !streamClosed) {
      const el = document.getElementById('status-msg');
      el.parentElement.parentElement.classList.add('error');
      el.textContent = '加载失败：EventSource 连接错误';
      es.close();
      finishProgress(false);
    }
  });

  // ===================================================================
  // builders
  // ===================================================================

  function makeEarth() {
    const g = new THREE.Group();
    const sphere = new THREE.Mesh(
      new THREE.SphereGeometry(R_EARTH_UNITS, 48, 32),
      new THREE.MeshPhongMaterial({ color: 0x2a78d2, emissive: 0x0a1830 })
    );
    g.add(sphere);
    return g;
  }

  function makeMoon() {
    const g = new THREE.Group();
    const sphere = new THREE.Mesh(
      new THREE.SphereGeometry(R_MOON_UNITS, 48, 32),
      new THREE.MeshPhongMaterial({
        map: makeMoonTexture(),
        color: 0xffffff,
        emissive: 0x101010,
        specular: 0x222222,
        shininess: 4,
      })
    );
    g.add(sphere);
    return g;
  }

  // Procedural equirectangular Moon texture (1024×512): grey highlands,
  // darker maria patches and scattered craters.  No external asset needed.
  function makeMoonTexture() {
    const W = 1024, H = 512;
    const canvas = document.createElement('canvas');
    canvas.width = W; canvas.height = H;
    const ctx = canvas.getContext('2d');

    // ---- base: mottled grey highlands ----
    const base = ctx.createLinearGradient(0, 0, 0, H);
    base.addColorStop(0, '#b8b4ae');
    base.addColorStop(0.5, '#c9c5bf');
    base.addColorStop(1, '#b0aca6');
    ctx.fillStyle = base;
    ctx.fillRect(0, 0, W, H);

    // fine grain noise
    const grain = ctx.createImageData(W, H);
    for (let i = 0; i < grain.data.length; i += 4) {
      const v = 128 + (Math.random() - 0.5) * 24;
      grain.data[i] = v; grain.data[i + 1] = v; grain.data[i + 2] = v;
      grain.data[i + 3] = 255;
    }
    ctx.putImageData(grain, 0, 0);

    // ---- maria: darker blotches (equirectangular → ellipse distortion) ----
    const maria = [
      { x: 0.30, y: 0.42, rx: 0.10, ry: 0.07, a: 0.35 },  // Mare Imbrium
      { x: 0.38, y: 0.55, rx: 0.08, ry: 0.06, a: 0.30 },  // Mare Serenitatis
      { x: 0.24, y: 0.60, rx: 0.09, ry: 0.05, a: 0.28 },  // Mare Tranquillitatis
      { x: 0.18, y: 0.48, rx: 0.07, ry: 0.05, a: 0.25 },  // Oceanus Procellarum
      { x: 0.48, y: 0.52, rx: 0.06, ry: 0.05, a: 0.22 },  // Mare Fecunditatis
      { x: 0.55, y: 0.40, rx: 0.05, ry: 0.04, a: 0.20 },  // Mare Crisium
      { x: 0.70, y: 0.45, rx: 0.06, ry: 0.05, a: 0.18 },  // Mare Australe
    ];
    for (const m of maria) {
      const g2 = ctx.createRadialGradient(
        m.x * W, m.y * H, 0,
        m.x * W, m.y * H, Math.max(m.rx * W, m.ry * H));
      g2.addColorStop(0, `rgba(70,66,62,${m.a})`);
      g2.addColorStop(1, 'rgba(70,66,62,0)');
      ctx.save();
      ctx.translate(m.x * W, m.y * H);
      ctx.scale(1, m.ry / m.rx);
      ctx.translate(-m.x * W, -m.y * H);
      ctx.fillStyle = g2;
      ctx.beginPath();
      ctx.arc(m.x * W, m.y * H, m.rx * W, 0, Math.PI * 2);
      ctx.fill();
      ctx.restore();
    }

    // ---- craters: darker floor + lighter rim ----
    for (let i = 0; i < 90; i++) {
      const cx = Math.random() * W;
      const cy = Math.random() * H;
      const r = 2 + Math.random() * 9;
      const shade = 0.5 + Math.random() * 0.5;
      ctx.beginPath();
      ctx.arc(cx, cy, r, 0, Math.PI * 2);
      ctx.fillStyle = `rgba(${90 * shade | 0},${88 * shade | 0},${84 * shade | 0},0.55)`;
      ctx.fill();
      ctx.beginPath();
      ctx.arc(cx, cy, r * 0.75, 0, Math.PI * 2);
      ctx.fillStyle = `rgba(${150 * shade | 0},${148 * shade | 0},${144 * shade | 0},0.35)`;
      ctx.fill();
    }

    const tex = new THREE.CanvasTexture(canvas);
    tex.wrapS = THREE.RepeatWrapping;
    tex.wrapT = THREE.ClampToEdgeWrapping;
    tex.anisotropy = 4;
    tex.needsUpdate = true;
    return tex;
  }

  function makeAxes() {
    const g = new THREE.Group();
    const len = 20;     // 20 units = 20,000 km for axis hint
    g.add(new THREE.ArrowHelper(new THREE.Vector3(1,0,0), new THREE.Vector3(), len, 0xff5577, 3, 2));
    g.add(new THREE.ArrowHelper(new THREE.Vector3(0,1,0), new THREE.Vector3(), len, 0x55ff77, 3, 2));
    g.add(new THREE.ArrowHelper(new THREE.Vector3(0,0,1), new THREE.Vector3(), len, 0x7799ff, 3, 2));
    // ---- axis labels (X→月 / Y→东 / Z→北) ----
    // Paper convention (Yang 1991 §3):
    //   X axis = 月地连线 (Moon→Earth), so +X points to Earth and -X points to Moon.
    //   Y axis = east in the equatorial plane (right-handed with X+Z north).
    //   Z axis = orbital-plane normal = north.
    // Place each sprite a bit beyond the arrow tip.
    const labelOffset = 1.15;
    const colorX = new THREE.Color(0xff5577);
    const colorY = new THREE.Color(0x55ff77);
    const colorZ = new THREE.Color(0x7799ff);
    const sx = makeLabelSprite('→月', colorX);
    sx.position.set(-len * labelOffset, 0, 0);
    sx.scale.set(8, 4, 1);
    g.add(sx);
    const sy = makeLabelSprite('→东', colorY);
    sy.position.set(0,  len * labelOffset, 0);
    sy.scale.set(8, 4, 1);
    g.add(sy);
    const sz = makeLabelSprite('→北', colorZ);
    sz.position.set(0, 0,  len * labelOffset);
    sz.scale.set(8, 4, 1);
    g.add(sz);
    return g;
  }

  function makeEarthMoonLine() {
    // Earth-Moon line: a thin dashed cylinder along -X of length ~DEFAULT_MOON_DIST_UNITS.
    const g = new THREE.Group();
    const dist = DEFAULT_MOON_DIST_UNITS;
    const mat = new THREE.LineDashedMaterial({
      color: 0x4ec9ff, dashSize: 4, gapSize: 4, opacity: 0.35, transparent: true,
    });
    const pts = [new THREE.Vector3(0, 0, 0), new THREE.Vector3(-dist, 0, 0)];
    const geo = new THREE.BufferGeometry().setFromPoints(pts);
    const line = new THREE.Line(geo, mat);
    line.computeLineDistances();
    g.add(line);
    g.name = 'earth_moon_line';
    g.userData.dist = dist;   // store for moon placement
    return g;
  }

  function makeEmptyLine() {
    // The Moon position evolves along the sidereal-mean-motion circle.  We
    // build a fixed placeholder; the per-frame update() places the Moon
    // mesh on it.  Here we just create a thin circle outline for context.
    const g = new THREE.Group();
    const segments = 96;
    const radius = DEFAULT_MOON_DIST_UNITS;
    const pts = [];
    for (let i = 0; i <= segments; i++) {
      const a = (i / segments) * Math.PI * 2;
      pts.push(new THREE.Vector3(Math.cos(a) * radius, 0, Math.sin(a) * radius));
    }
    const geo = new THREE.BufferGeometry().setFromPoints(pts);
    const mat = new THREE.LineDashedMaterial({
      color: 0x4ec9ff, dashSize: 2, gapSize: 4, opacity: 0.18, transparent: true,
    });
    const line = new THREE.Line(geo, mat);
    line.computeLineDistances();
    g.add(line);
    return g;
  }

  function makeGrid() {
    const grid = new THREE.GridHelper(800, 16, 0x335577, 0x223344);
    grid.position.y = 0;          // XZ plane grid by default (Y up)
    grid.rotation.x = Math.PI / 2; // put into XY plane (Z up)
    grid.material.opacity = 0.4;
    grid.material.transparent = true;
    return grid;
  }

  // ---------------------- body graticules (经纬度网格) ----------------
  // R-1V-e follow-up: 给 Earth / Moon 加经纬度网格，让"东西南北"在 3D
  // 视角里一目了然（之前用户反映"东西方向变成了南北"——核心原因是
  // 缺少方向参考系）。经纬度约定遵循论文 §3：
  //   λ (东向经度) 从 +X 轴在 XY-plane 内度量度（绕 +Z 轴）
  //   φ (北纬) 从赤道面（XY）量起向北（+Z 方向）
  // 单位球面 → MEP XYZ 的转换（Yang 1991 §3 / PRD §S-0）:
  //   x = R cos(φ) cos(λ),  y = R cos(φ) sin(λ),  z = R sin(φ)
  // 参数 lat_deg ∈ [-90, 90] (北正), lon_deg ∈ [-180, 180] (东正).
  function _latLonToXyz(lat_deg, lon_deg, R) {
    const lat = lat_deg * Math.PI / 180.0;
    const lon = lon_deg * Math.PI / 180.0;
    return new THREE.Vector3(
      R * Math.cos(lat) * Math.cos(lon),
      R * Math.cos(lat) * Math.sin(lon),
      R * Math.sin(lat),
    );
  }

  function _makeGraticule(R, opts = {}) {
    // R-1V-e: 通用经纬度网格生成器.
    // 返回一个 THREE.Group，包含：
    //   1) 赤道 (φ=0)
    //   2) 0° 子午线 (λ=0)
    //   3) ±90° / 180° / 270° 子午线
    //   4) ±30° / ±60° 纬圈
    // 颜色默认: opts.equatorHex / meridian0Hex / meridianHex / latitudeHex.
    // opts.equatorHex   default '#4ec9ff' (亮蓝)
    // opts.meridian0Hex default '#4ec9ff'
    // opts.meridianHex  default '#2a78d2' (中蓝)
    // opts.latitudeHex  default '#335577' (浅蓝灰)
    const o = {
      equatorHex:   opts.equatorHex   ?? '#4ec9ff',
      meridian0Hex: opts.meridian0Hex ?? '#4ec9ff',
      meridianHex:  opts.meridianHex  ?? '#2a78d2',
      latitudeHex:  opts.latitudeHex  ?? '#335577',
      nSeg:         opts.nSeg         ?? 96,
    };
    const g = new THREE.Group();

    function pushLine(points, hex, opacity) {
      const geo = new THREE.BufferGeometry().setFromPoints(points);
      const mat = new THREE.LineBasicMaterial({
        color: new THREE.Color(hex), transparent: true, opacity,
      });
      g.add(new THREE.Line(geo, mat));
    }

    // ---- 1. 赤道 (φ=0) -------------------------------------------------
    const equator = [];
    for (let i = 0; i <= o.nSeg; i++) {
      const lon = (i / o.nSeg) * 2 * Math.PI;
      equator.push(_latLonToXyz(0, lon * 180 / Math.PI, R));
    }
    pushLine(equator, o.equatorHex, 0.85);

    // ---- 2. 0° 子午线（从北极到南极，亮色）----------------------------
    const meridian0 = [];
    for (let i = 0; i <= o.nSeg; i++) {
      const lat = -90 + (i / o.nSeg) * 180;  // -90° → +90°
      meridian0.push(_latLonToXyz(lat, 0, R));
    }
    pushLine(meridian0, o.meridian0Hex, 0.85);

    // ---- 3. ±90° / 180° / 270° 子午线（中等亮度）-----------------------
    for (const lonDeg of [90, 180, 270]) {
      const m = [];
      for (let i = 0; i <= o.nSeg; i++) {
        const lat = -90 + (i / o.nSeg) * 180;
        m.push(_latLonToXyz(lat, lonDeg, R));
      }
      pushLine(m, o.meridianHex, 0.55);
    }

    // ---- 4. ±30° / ±60° 纬圈（浅色）-----------------------------------
    for (const latDeg of [30, 60, -30, -60]) {
      const lat = [];
      for (let i = 0; i <= o.nSeg; i++) {
        const lon = (i / o.nSeg) * 2 * Math.PI;
        lat.push(_latLonToXyz(latDeg, lon * 180 / Math.PI, R));
      }
      pushLine(lat, o.latitudeHex, 0.45);
    }

    return g;
  }

  function makeEarthGraticule() {
    // R-1V-e: 经纬度网格 (赤道 + 0/90/180/270° 子午线 + ±30°/±60° 纬圈).
    // convention per Yang 1991 §3 / PRD §S-0.
    // 半径 +0.5% 让网格浮在球体表面之上避免 z-fighting.
    return _makeGraticule(R_EARTH_UNITS * 1.005, {
      nSeg: 96,
    });
  }

  function makeMoonGraticule() {
    // R-1V-e follow-up: 月球经纬度网格.
    // 月球半径较小, nSeg 64 已满足视觉精度.
    // 注意: 这只是辅助网格, 论文 (Yang 1991 §5.1) 假设月球赤道面与月球轨道面重合,
    // 所以月面经度原点定义为面向地球的子子点; 这里用同一种 mapping.
    return _makeGraticule(R_MOON_UNITS * 1.005, {
      nSeg: 64,
    });
  }

  // 月球自转轴 + 南北极标注。MEP 系中 Z 轴 = 轨道面法向 = 北（论文 §3 约定，
  // 月球赤道面与轨道面重合），所以自转轴沿 Z 方向，北极在 +Z、南极在 -Z。
  function makeMoonAxis() {
    const g = new THREE.Group();
    const len = R_MOON_UNITS * 1.6;   // 伸出球面之外
    const axisPts = [new THREE.Vector3(0, 0, -len), new THREE.Vector3(0, 0, len)];
    const axisGeo = new THREE.BufferGeometry().setFromPoints(axisPts);
    const axisMat = new THREE.LineBasicMaterial({ color: 0xffd166, transparent: true, opacity: 0.9 });
    g.add(new THREE.Line(axisGeo, axisMat));

    // 北极（+Z）：黄色圆点 + “北” 标签
    const north = new THREE.Mesh(
      new THREE.SphereGeometry(0.25, 12, 12),
      new THREE.MeshBasicMaterial({ color: 0xffd166 })
    );
    north.position.set(0, 0, len);
    g.add(north);
    const nLabel = makeLabelSprite('北', new THREE.Color(0xffd166), 5);
    nLabel.position.set(0, 0, len + 1.4);
    g.add(nLabel);

    // 南极（-Z）：蓝色圆点 + “南” 标签
    const south = new THREE.Mesh(
      new THREE.SphereGeometry(0.25, 12, 12),
      new THREE.MeshBasicMaterial({ color: 0x66ccff })
    );
    south.position.set(0, 0, -len);
    g.add(south);
    const sLabel = makeLabelSprite('南', new THREE.Color(0x66ccff), 5);
    sLabel.position.set(0, 0, -len - 1.4);
    g.add(sLabel);

    return g;
  }

  // ---------------------- orbit objects ----------------------------
  /**
   * @typedef {Object} OrbitInfo
   * @property {string} family
   * @property {number} orbit
   * @property {THREE.Line} line
   * @property {THREE.Points} dots
   * @property {THREE.Sprite} label
   * @property {THREE.Group} endpoints
   * @property {string} colorHex
   * @property {Object} meta
   */
  /** @type {OrbitInfo[]} */
  const orbitObjects = [];

  // Family → list-element id (left-panel section).
  const FAMILY_LIST_ID = {
    coplanar_table1:       'list-coplanar',
    non_coplanar_table2:   'list-noncoplanar',
    same_direction_table3: 'list-samedir',
    polar_table4:          'list-polar',
  };
  // Family → short letter used in orbit labels (#C1, #N1, #S1, #P1).
  const FAMILY_SHORT = {
    coplanar_table1:       'C',
    non_coplanar_table2:   'N',
    same_direction_table3: 'S',
    polar_table4:          'P',
  };
  // Family → human-readable prefix shown before the orbit number.
  const FAMILY_LABEL_PREFIX = {
    coplanar_table1:       (m) => `T${m.T_h}h · `,
    non_coplanar_table2:   (m) => `i₀=${m.i0_deg}° · `,
    same_direction_table3: (m) => `i₀=${m.i0_deg}° · `,
    polar_table4:          (m) => `i₀=${m.i0_deg}° · `,
  };

  function buildOneOrbit(rec) {
    const colorHex = familyColor(rec.family, rec.orbit);
    const color = new THREE.Color(colorHex);

    // trajectory line
    const rMep = rec.r_mep.map(([x, y, z]) => new THREE.Vector3(x / KM_PER_UNIT, y / KM_PER_UNIT, z / KM_PER_UNIT));
    const lineGeo = new THREE.BufferGeometry().setFromPoints(rMep);
    const lineMat = new THREE.LineBasicMaterial({ color, linewidth: 1, opacity: 0.85, transparent: true });
    const line = new THREE.Line(lineGeo, lineMat);

    // sample dots (small)
    const dotGeo = new THREE.BufferGeometry().setFromPoints(rMep);
    const dotMat = new THREE.PointsMaterial({ color, size: 1.6, sizeAttenuation: false });
    const dots = new THREE.Points(dotGeo, dotMat);

    // ---- surface markers (replaces the old end-marker yellow ball) ----
    // perigee on Earth surface (geocentric lon/lat from PRD §S-0).
    // Table 1 coplanar 12 entries have φ_p_ref = 0 (the paper omitted the
    // launch-site latitude) — fall back to Xichang φ = 28.1° N.
    // Table 4 polar entries already have φ_p_ref ∈ ±1.3°.
    const familyShort = FAMILY_SHORT[rec.family] || '?';
    const orbitLabel  = `${familyShort}${rec.orbit}`;
    const perigeeGroup = new THREE.Group();
    const phi_p_ref = (rec.perigee_lat_ref_deg != null) ? rec.perigee_lat_ref_deg : 0.0;
    const phi_p_use = (Math.abs(phi_p_ref) < 0.5) ? 28.1 : phi_p_ref;     // Xichang fallback
    const lambda_p_use = (rec.perigee_lon_ref_deg != null) ? rec.perigee_lon_ref_deg : rec.perigee_lon_sol_deg;
    const perigeePos = _latLonToXyz(phi_p_use, lambda_p_use, R_EARTH_UNITS * 1.005);
    perigeeGroup.add(new THREE.Mesh(
      new THREE.SphereGeometry(0.8, 16, 12),
      new THREE.MeshBasicMaterial({ color })
    ).clone().translateX(perigeePos.x).translateY(perigeePos.y).translateZ(perigeePos.z));
    perigeeGroup.children[perigeeGroup.children.length - 1].position.copy(perigeePos);
    // orbit # label at perigee (small sprite so it doesn't cover the Earth)
    perigeeGroup.add(makeLabelSprite(`#${orbitLabel}`, color, 3).clone()
      .translateX(perigeePos.x).translateY(perigeePos.y + 3).translateZ(perigeePos.z));
    perigeeGroup.children[perigeeGroup.children.length - 1].position.copy(
      new THREE.Vector3(perigeePos.x, perigeePos.y + 3, perigeePos.z));
    // lon/lat sprite: show Xichang label when fallback was applied
    const isXichang = Math.abs(phi_p_ref) < 0.5;
    const phiLabel = isXichang ? `${phi_p_use.toFixed(1)}°N 西昌` : `${phi_p_use.toFixed(2)}°`;
    const periSprite = makeLabelSprite(
      `λ=${lambda_p_use.toFixed(2)}°\nφ=${phiLabel}`,
      new THREE.Color(isXichang ? 0xffd166 : 0xeeeeee), 3);
    periSprite.position.copy(new THREE.Vector3(perigeePos.x, perigeePos.y - 4, perigeePos.z));
    perigeeGroup.add(periSprite);
    layers.earth.add(perigeeGroup);

    // perilune on Moon surface (selenocentric lon/lat).
    // Use the paper/test-case target (perilune_lon_ref_deg / perilune_lat_ref_deg)
    // so the endpoint marker lands exactly on the specified近月点经纬度.
    // For the four i_m=90° polar orbits the solved values collapse onto one
    // point (see mode_b_validation.json), so the ref target is authoritative.
    // Fall back to the solved lon/lat only when a ref value is missing
    // (e.g. coplanar Table-1 rows, which publish no lunar lon/lat).
    const periluneGroup = new THREE.Group();
    const lambda_m_use = (rec.perilune_lon_ref_deg != null)
      ? rec.perilune_lon_ref_deg
      : ((rec.perilune_lon_sol_deg != null) ? rec.perilune_lon_sol_deg : 0.0);
    const phi_m_use = (rec.perilune_lat_ref_deg != null)
      ? rec.perilune_lat_ref_deg
      : ((rec.perilune_lat_sol_deg != null) ? rec.perilune_lat_sol_deg : 0.0);
    const perilunePos = _latLonToXyz(phi_m_use, lambda_m_use, R_MOON_UNITS * 1.005);
    periluneGroup.add(new THREE.Mesh(
      new THREE.SphereGeometry(0.45, 12, 12),
      new THREE.MeshBasicMaterial({ color })
    ).clone().translateX(perilunePos.x).translateY(perilunePos.y).translateZ(perilunePos.z));
    periluneGroup.children[periluneGroup.children.length - 1].position.copy(perilunePos);
    periluneGroup.add(makeLabelSprite(`#${orbitLabel}`, color, 3).clone()
      .translateX(perilunePos.x).translateY(perilunePos.y + 1.6).translateZ(perilunePos.z));
    periluneGroup.children[periluneGroup.children.length - 1].position.copy(
      new THREE.Vector3(perilunePos.x, perilunePos.y + 1.6, perilunePos.z));
    layers.moon.add(periluneGroup);

    // ---- floating sprite label at the trajectory midpoint ----
    // Smaller than before (8 units ≈ 8,000 km) so it no longer covers the
    // Moon/Earth disks when the camera is far away.
    const label = makeLabelSprite(`#${orbitLabel}`, color, 8);

    const group = new THREE.Group();
    group.add(line);
    group.add(dots);
    group.add(label);
    scene.add(group);

    return {
      family: rec.family, orbit: rec.orbit,
      line, dots, label, group,
      surfaceMarkers: { earth: perigeeGroup, moon: periluneGroup },
      colorHex, meta: rec,
    };
  }

  function familyColor(family, orbit) {
    // Coplanar: warm hues (orange-amber), 4 hues × 3 times
    // Non-coplanar: green-mint hues
    // Same-direction: yellow-gold hues
    // Polar: cool hues (magenta-cyan), 4 distinct
    const coplanarPal = ['#ffa552', '#ff7b3a', '#ffc266', '#ff9c40',
                         '#ffb060', '#ff8a4a', '#ffd07a', '#ff9f55',
                         '#ffba70', '#ff955c', '#ffdc8e', '#ffaa66'];
    if (family === 'coplanar_table1') {
      return coplanarPal[orbit - 1] || '#ffaa66';
    }
    if (family === 'non_coplanar_table2') {
      const pal = ['#66ddaa', '#4ec99c', '#7be8b8', '#3dbf90'];
      return pal[orbit - 1] || '#66ddaa';
    }
    if (family === 'same_direction_table3') {
      const pal = ['#ffe066', '#ffd54a', '#ffe98f', '#ffcc33'];
      return pal[orbit - 1] || '#ffe066';
    }
    // polar
    const polarPal = ['#b07cff', '#8a6cff', '#6c8aff', '#7cb6ff'];
    return polarPal[orbit - 1] || '#9090ff';
  }

  function makeLabelSprite(text, color, size = 20) {
    // Multi-line aware: split on '\n' and size the canvas to fit all lines.
    // depthTest:true so labels are occluded by the Earth/Moon spheres instead
    // of floating through them; depthWrite:false so sprites never occlude
    // other geometry.
    const lines = String(text).split('\n');
    const lineH = 34;                 // px per line
    const padX = 14, padY = 10;
    const canvas = document.createElement('canvas');
    canvas.width = 128;
    canvas.height = padY * 2 + lineH * lines.length;
    const ctx = canvas.getContext('2d');
    ctx.fillStyle = 'rgba(8,12,22,0.85)';
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.strokeStyle = '#' + color.getHexString();
    ctx.lineWidth = 2;
    ctx.strokeRect(1, 1, canvas.width - 2, canvas.height - 2);
    ctx.fillStyle = '#' + color.getHexString();
    ctx.font = 'bold 30px -apple-system, sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    lines.forEach((ln, i) => {
      ctx.fillText(ln, canvas.width / 2, padY + lineH * (i + 0.5));
    });
    const tex = new THREE.CanvasTexture(canvas);
    tex.needsUpdate = true;
    const mat = new THREE.SpriteMaterial({ map: tex, transparent: true, depthTest: true, depthWrite: false });
    const sprite = new THREE.Sprite(mat);
    // Keep the world-space size proportional to the canvas aspect ratio so
    // multi-line labels don't get squashed.
    const baseW = size;
    sprite.scale.set(baseW, baseW * (canvas.height / canvas.width), 1);
    sprite.position.set(0, 0, 0);   // positioned per-frame in updateLabels()
    sprite.userData.baseText = text;
    return sprite;
  }

  // ---------------------- control panel -----------------------------
  /** @type {HTMLInputElement[]} */
  const orbitCheckboxes = [];

  function addControlRow(info) {
    const target = document.getElementById(FAMILY_LIST_ID[info.family] || 'list-coplanar');
    const row = document.createElement('label');
    row.className = 'row';
    const labelPrefix = (FAMILY_LABEL_PREFIX[info.family] || (() => ''))(info.meta);
    const labelShort = FAMILY_SHORT[info.family] + info.orbit;
    row.innerHTML = `
      <input type="checkbox" checked data-key="${info.family}::${info.orbit}">
      <span class="swatch" style="background:${info.colorHex}"></span>
      <span class="name">${labelPrefix}
        #${labelShort}
        <span class="meta">v_p=${info.meta.vp_ref_km_s.toFixed(3)} km/s</span>
      </span>`;
    const checkbox = row.querySelector('input');
    orbitCheckboxes.push(checkbox);
    checkbox.addEventListener('change', e => {
      updateOrbitVisibility(info, e.target.checked);
    });
    row.addEventListener('click', e => {
      if (e.target.tagName === 'INPUT') return;
      SELECTED = { family: info.family, orbit: info.orbit };
      showMeta(info);
    });
    target.appendChild(row);
  }

  // 全选 / 全不选：同步所有轨道 checkbox 与 3D 可见性。
  function setAllOrbitsVisible(visible) {
    orbitCheckboxes.forEach(cb => { cb.checked = visible; });
    orbitObjects.forEach(o => updateOrbitVisibility(o, visible));
  }
  document.getElementById('select-all').addEventListener('click', () => setAllOrbitsVisible(true));
  document.getElementById('select-none').addEventListener('click', () => setAllOrbitsVisible(false));

  // 左侧面板收起/展开：只隐藏轨道列表，保留标题栏与「展开」按钮。
  const controlsPanel = document.getElementById('controls');
  const collapseBtn = document.getElementById('collapse-controls');
  collapseBtn.addEventListener('click', () => {
    const collapsed = controlsPanel.classList.toggle('collapsed');
    collapseBtn.textContent = collapsed ? '»' : '«';
    collapseBtn.title = collapsed ? '展开面板' : '收起面板';
  });

  // 单一入口：根据「轨道 checkbox 状态 + 图层开关」计算一条轨道的完整可见性。
  // 轨道组（轨迹线/采样点/中点标签）由 checkbox 控制；地表标记（近地点/近月点
  // 圆点及其标签）挂在 layers.earth / layers.moon 下，因此还要叠加对应天体图层
  // 与 labels 开关。
  function updateOrbitVisibility(info, orbitOn) {
    info.group.visible = orbitOn;
    info.label.visible = layerVisible.labels && orbitOn;
    if (info.surfaceMarkers) {
      info.surfaceMarkers.earth.visible = layerVisible.earth && orbitOn;
      info.surfaceMarkers.moon.visible = layerVisible.moon && orbitOn;
      info.surfaceMarkers.earth.traverse(obj => {
        if (obj.isSprite) obj.visible = layerVisible.labels && orbitOn;
      });
      info.surfaceMarkers.moon.traverse(obj => {
        if (obj.isSprite) obj.visible = layerVisible.labels && orbitOn;
      });
    }
  }

  function showMeta(info) {
    const m = info.meta;
    const t = (k) => Array.isArray(m[k]) ? `[${m[k].length} pts]` : String(m[k]);
    const fmt = (v, d=2) => (typeof v === 'number') ? v.toFixed(d) : String(v);
    const text =
      `family    : ${m.family}\n` +
      `orbit     : #${m.orbit}\n` +
      `T         : ${m.T_h} h\n` +
      `i₀        : ${fmt(m.i0_deg, 2)}°\n` +
      `v_p(ref)  : ${fmt(m.vp_ref_km_s, 6)} km/s\n` +
      `i_m(ref)  : ${fmt(m.im_ref_deg, 1)}°\n` +
      `λ_p(ref)  : ${fmt(m.perigee_lon_ref_deg, 3)}°\n` +
      `λ_m(ref)  : ${m.perilune_lon_ref_deg !== null ? fmt(m.perilune_lon_ref_deg, 3) + '°' : 'n/a'}\n` +
      `φ_m(ref)  : ${m.perilune_lat_ref_deg !== null ? fmt(m.perilune_lat_ref_deg, 3) + '°' : 'n/a'}\n` +
      `p_sol     : (${m.p_sol.map(x => fmt(x,4)).join(', ')})\n` +
      `q         : (${m.q.map(x => fmt(x,6)).join(', ')})\n` +
      `i_m(sol)  : ${fmt(m.i_m_sol_deg, 4)}°\n` +
      `v_perilune: ${fmt(m.v_perilune_km_s, 6)} km/s\n` +
      `C3        : ${fmt(m.C3, 4)} km²/s²\n` +
      `samples   : ${t('t_s')} per orbit`;
    document.getElementById('meta-body').textContent = text;
  }

  // ---------------------- layer toggles ----------------------------
  document.querySelectorAll('#meta button[data-view]').forEach(btn => {
    btn.addEventListener('click', () => setView(btn.dataset.view));
  });
  document.querySelectorAll('#meta button[data-layer]').forEach(btn => {
    btn.addEventListener('click', () => {
      const k = btn.dataset.layer;
      layerVisible[k] = !layerVisible[k];
      btn.classList.toggle('on', layerVisible[k]);
      applyLayerVisibility();
    });
  });

  function applyLayerVisibility() {
    layers.earth.visible = layerVisible.earth;
    layers.moon.visible = layerVisible.moon;
    layers.axis.visible = layerVisible.axis;
    layers.moonline.visible = layerVisible.moonline;
    layers.moontrack.visible = layerVisible.moontrack;
    layers.grid.visible = layerVisible.grid;
    // Graticules: earth graticule is a scene child; moon graticule is a child of
    // layers.moon (so moon visibility acts as a master switch for it).
    layers.earthgraticule.visible = layerVisible.earthgraticule;
    layers.moongraticule.visible = layerVisible.moongraticule && layerVisible.moon;
    // Moon spin axis + N/S labels: child of layers.moon, so it also needs the
    // moon master switch.
    layers.moonaxis.visible = layerVisible.moonaxis && layerVisible.moon;
    // Recompute each orbit's full visibility from its checkbox + layer toggles.
    orbitObjects.forEach((o, i) => {
      const orbitOn = orbitCheckboxes[i] ? orbitCheckboxes[i].checked : o.group.visible;
      updateOrbitVisibility(o, orbitOn);
    });
  }

  function setView(kind) {
    const dist = 600;
    if (kind === 'from_moon_pov') { camera.position.set(-dist, 0, 0); controls.target.set(0,0,0); }
    else if (kind === 'from_earth_pov') { camera.position.set( dist, 0, 0); controls.target.set(0,0,0); }
    else if (kind === 'xy')        { camera.position.set(0, 0, dist); controls.target.set(0,0,0); }
    else if (kind === 'xz')   { camera.position.set(0, dist, 0); controls.target.set(0,0,0); }
    else if (kind === 'yz')   { camera.position.set(dist, 0, 0); controls.target.set(0,0,0); }
    else if (kind === 'iso')  { camera.position.set(450, 350, 450); controls.target.set(0,0,0); }
    else                       { /* free: keep current */ }
    camera.lookAt(controls.target);
    // Re-sync the orbit-controls internal spherical state so the next
    // wheel/drag continues from the preset view instead of snapping back
    // to the pre-preset zoom/angle.
    controls.syncFromCamera();
    controls.update();
  }

  // ---------------------- per-frame --------------------------------
  const moonDist = DEFAULT_MOON_DIST_UNITS;

  function tick() {
    requestAnimationFrame(tick);
    // The MEP frame is frozen at t_arrive, and every sampled trajectory ends
    // at the Moon's arrival position.  Pin the Moon at the frozen anti-X
    // arrival endpoint instead of animating it away, so the transfer-orbit
    // endpoints always land on (or just above) the Moon's surface.
    layers.moon.position.set(-moonDist, 0, 0);
    controls.update();
    renderer.render(scene, camera);
    updateLabels();
  }

  function updateLabels() {
    orbitObjects.forEach(o => {
      if (!o.group.visible) return;
      // Position the label at the midpoint of the trajectory, projected to
      // a constant screen-relative offset so it doesn't overlap the line.
      const m = o.meta;
      const mid = m.r_mep[Math.floor(m.r_mep.length / 2)];
      o.label.position.set(mid[0] / KM_PER_UNIT, mid[1] / KM_PER_UNIT, mid[2] / KM_PER_UNIT);
    });
  }

  tick();

  // ---------------------- resize -----------------------------------
  window.addEventListener('resize', () => {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth, window.innerHeight);
  });

  // ===================================================================
  // Minimal OrbitControls (avoids needing an extra vendored file).
  // Implements: left-drag rotate, right-drag pan, wheel zoom.
  // ===================================================================
  function makeOrbitControls(cam, dom) {
    const target = new THREE.Vector3(0, 0, 0);
    const spherical = new THREE.Spherical();
    spherical.setFromVector3(cam.position.clone().sub(target));

    let dragging = false, panning = false, lastX = 0, lastY = 0;

    const onDown = (e) => {
      dragging = true;
      panning = (e.button === 2 || e.shiftKey);
      lastX = e.clientX; lastY = e.clientY;
    };
    const onMove = (e) => {
      if (!dragging) return;
      const dx = e.clientX - lastX;
      const dy = e.clientY - lastY;
      lastX = e.clientX; lastY = e.clientY;
      if (panning) {
        const panSpeed = spherical.radius * 0.0015;
        const right = new THREE.Vector3().crossVectors(cam.position.clone().sub(target), new THREE.Vector3(0,0,1)).normalize();
        const up = new THREE.Vector3(0, 0, 1);
        target.addScaledVector(right, -dx * panSpeed);
        target.addScaledVector(up,    dy * panSpeed);
      } else {
        spherical.theta -= dx * 0.005;
        spherical.phi   -= dy * 0.005;
        spherical.phi = Math.max(0.05, Math.min(Math.PI - 0.05, spherical.phi));
      }
      updateCam();
    };
    const onUp = () => { dragging = false; panning = false; };
    const onWheel = (e) => {
      e.preventDefault();
      const factor = Math.exp(e.deltaY * 0.001);
      spherical.radius = Math.max(20, Math.min(1500, spherical.radius * factor));
      updateCam();
    };

    function updateCam() {
      const v = new THREE.Vector3().setFromSpherical(spherical);
      cam.position.copy(target.clone().add(v));
      cam.lookAt(target);
    }

    // Re-derive the internal spherical state from the current camera position.
    // Needed after a view preset moves the camera directly (setView), so the
    // next wheel/drag event continues from the preset instead of snapping
    // back to the old zoom/angle.
    function syncFromCamera() {
      spherical.setFromVector3(cam.position.clone().sub(target));
    }

    dom.addEventListener('mousedown', onDown);
    window.addEventListener('mousemove', onMove);
    window.addEventListener('mouseup', onUp);
    dom.addEventListener('wheel', onWheel, { passive: false });
    dom.addEventListener('contextmenu', e => e.preventDefault());

    return {
      get target() { return target; },
      set target(v) { target.copy(v); updateCam(); },
      update() { /* no-op; updated on every event */ },
      syncFromCamera,
    };
  }
})();

/* ---- 访问统计上报（与 viewer/static/app.js 同口径） ---- */
(function () {
  var url = '/api/visit?path=' + encodeURIComponent(location.pathname);
  if (navigator.sendBeacon) navigator.sendBeacon(url);
  else fetch(url, { method: 'POST', keepalive: true }).catch(function () {});
})();
