/* orbits.js — 杨维廉论文 24 条地月转移轨道 XYZ 三维可视化
 *
 * 数据：GET /api/yang_orbits_stream (SSE: meta → 24 × orbit → done)
 * 坐标系：MEP (Moon-Earth Perifocal) 冻结于到达时刻 t_arrive
 *   X = 月地连线，由月→地；Y / Z 右手系完备（Z = 轨道面法向）
 *   原点 = 地心；单位约定：1 three.js unit = 1000 km
 *
 * 本文件在原 AgentMoons orbits_3d.html 移植版基础上重写界面：
 *   三栏布局（轨道列表 / 3D 视口 / 详情与图表）、搜索过滤、族谱分组、
 *   任务时间轴动画（航天器沿采样轨迹推进、月球沿真实采样轨迹移动）、
 *   C3 柱状图（可点击选择）、聚焦模式、相机聚焦、键盘快捷键。
 */
(() => {
  'use strict';

  // =====================================================================
  // 常量与元数据
  // =====================================================================
  const KM_PER_UNIT = 1000.0;
  const R_EARTH_UNITS = 6378.137 / KM_PER_UNIT;
  const R_MOON_UNITS = 1737.4 / KM_PER_UNIT;
  const DEFAULT_MOON_DIST_UNITS = 384400 / KM_PER_UNIT;
  const ANIM_DURATION_S = 12.0;   // 1× 速度下跑完全程的秒数

  const FAMILIES = {
    coplanar_table1: {
      short: 'C', label: '表 1 · 共面转移', cls: 'coplanar', color: '#ffa552',
      palette: ['#ffa552', '#ff7b3a', '#ffc266', '#ff9c40', '#ffb060', '#ff8a4a',
                '#ffd07a', '#ff9f55', '#ffba70', '#ff955c', '#ffdc8e', '#ffaa66'],
    },
    non_coplanar_table2: {
      short: 'N', label: '表 2 · 非共面转移', cls: 'noncoplanar', color: '#66ddaa',
      palette: ['#66ddaa', '#4ec99c', '#7be8b8', '#3dbf90'],
    },
    same_direction_table3: {
      short: 'S', label: '表 3 · 同向发射', cls: 'samedir', color: '#ffe066',
      palette: ['#ffe066', '#ffd54a', '#ffe98f', '#ffcc33'],
    },
    polar_table4: {
      short: 'P', label: '表 4 · 极月卫星', cls: 'polar', color: '#b07cff',
      palette: ['#b07cff', '#8a6cff', '#6c8aff', '#7cb6ff'],
    },
  };
  const FAMILY_ORDER = ['coplanar_table1', 'non_coplanar_table2',
                        'same_direction_table3', 'polar_table4'];

  const LAYERS = [
    ['earth', '地球'], ['moon', '月球'], ['axis', 'XYZ 轴'],
    ['moonline', '地月连线'], ['moontrack', '月球轨迹'], ['labels', '轨道编号'],
    ['earthgraticule', '地球经纬'], ['moongraticule', '月球经纬'],
    ['moonaxis', '月球自转轴'], ['grid', 'XYZ 网格'], ['endpoints', '起末点'],
  ];

  const layerVisible = {
    earth: true, moon: true, axis: true, moonline: true, moontrack: true,
    labels: true, earthgraticule: true, moongraticule: true, moonaxis: true,
    grid: false, endpoints: true,
  };

  // =====================================================================
  // 场景
  // =====================================================================
  const stage = document.getElementById('stage');
  const scene = new THREE.Scene();
  scene.background = new THREE.Color(0x05070c);

  const camera = new THREE.PerspectiveCamera(45, 1, 0.05, 4000);
  camera.position.set(-600, 0, 0);

  const renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
  renderer.setSize(stage.clientWidth || 1, stage.clientHeight || 1);
  stage.appendChild(renderer.domElement);

  const controls = makeOrbitControls(camera, renderer.domElement);

  scene.add(new THREE.AmbientLight(0xffffff, 0.55));
  const sun = new THREE.DirectionalLight(0xffffff, 0.85);
  sun.position.set(200, 200, 200);
  scene.add(sun);
  const fill = new THREE.DirectionalLight(0x8899bb, 0.30);
  fill.position.set(-200, -150, -200);
  scene.add(fill);
  scene.add(new THREE.HemisphereLight(0x8899bb, 0x223344, 0.35));

  // ---- 参考天体与图层 ----
  const layers = {};
  layers.earth = makeEarth();
  layers.moon = makeMoon();
  layers.axis = makeAxes();
  layers.moonline = makeEarthMoonLine();
  layers.moontrack = makeEmptyLine();
  layers.grid = makeGrid();
  layers.earthgraticule = makeEarthGraticule();
  layers.moongraticule = makeMoonGraticule();
  layers.labels = new THREE.Group();
  layers.endpoints = new THREE.Group();

  layers.moon.add(layers.moongraticule);
  layers.moonaxis = makeMoonAxis();
  layers.moon.add(layers.moonaxis);

  for (const k of Object.keys(layers)) {
    layers[k].name = 'layer_' + k;
    scene.add(layers[k]);
  }

  // ---- 动画对象：航天器标记 / 已飞轨迹 / 选中轨道的月球轨迹 ----
  const scMarker = new THREE.Mesh(
    new THREE.SphereGeometry(1.1, 16, 12),
    new THREE.MeshBasicMaterial({ color: 0xffffff })
  );
  scMarker.visible = false;
  scene.add(scMarker);

  const trailLine = new THREE.Line(
    new THREE.BufferGeometry(),
    new THREE.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.95 })
  );
  trailLine.visible = false;
  scene.add(trailLine);

  const moonTrackSel = new THREE.Line(
    new THREE.BufferGeometry(),
    new THREE.LineBasicMaterial({ color: 0x58a6ff, transparent: true, opacity: 0.5 })
  );
  moonTrackSel.visible = false;
  scene.add(moonTrackSel);

  // =====================================================================
  // 状态
  // =====================================================================
  let META = { t_arrive_utc: '', n_samples_per_orbit: 121, n_orbits_expected: 24 };
  let ORBITS = [];
  const orbitObjects = [];
  let selected = null;
  let focusMode = true;   // 默认聚焦选中轨道，其余变暗
  const anim = { playing: false, frac: 1, lastTs: 0 };

  // =====================================================================
  // 数据加载（SSE）
  // =====================================================================
  const progressEl = document.getElementById('progress');
  const progressFill = document.getElementById('progress-fill');
  const progressCur = document.getElementById('progress-current');
  const progressTot = document.getElementById('progress-total');
  const progressLabel = document.getElementById('progress-label');
  const statusEl = document.getElementById('status-msg');

  function showProgress(label, cur, total) {
    progressLabel.textContent = label;
    progressCur.textContent = String(cur);
    progressTot.textContent = String(total);
    progressFill.style.width = total > 0 ? (cur / total) * 100 + '%' : '6%';
    progressEl.classList.remove('hidden');
  }
  showProgress('加载轨道数据…', 0, 0);

  let received = 0;
  let streamClosed = false;

  function finishProgress(ok, msg) {
    if (streamClosed) return;
    streamClosed = true;
    progressFill.style.width = '100%';
    if (!ok) {
      progressLabel.textContent = msg || '加载失败';
      progressFill.style.background = 'var(--error)';
      statusEl.textContent = msg || '加载失败';
    }
    setTimeout(() => progressEl.classList.add('hidden'), 400);
  }

  const es = new EventSource('/api/yang_orbits_stream');

  es.addEventListener('meta', (ev) => {
    let meta;
    try { meta = JSON.parse(ev.data); } catch (e) { meta = {}; }
    META = Object.assign({}, META, meta);
    statusEl.innerHTML = '';
    const strong = document.createElement('strong');
    strong.textContent = String(META.n_orbits_expected || 24);
    statusEl.appendChild(strong);
    statusEl.appendChild(document.createTextNode(
      ' 条轨道 · 每条 ' + META.n_samples_per_orbit + ' 采样 · t_arrive = '
      + META.t_arrive_utc + ' · H_m = ' + META.H_m_km + ' km'));
    showProgress('正在构建轨道…', 0, META.n_orbits_expected || 24);
  });

  es.addEventListener('orbit', (ev) => {
    let rec;
    try { rec = JSON.parse(ev.data); } catch (e) {
      console.error('[orbits] 无法解析 orbit 事件', e);
      return;
    }
    ORBITS.push(rec);
    if (!rec.error) {
      try {
        const info = buildOneOrbit(rec);
        orbitObjects.push(info);
      } catch (e) {
        console.error('[orbits] 构建轨道失败:', rec, e);
      }
    }
    received += 1;
    progressCur.textContent = String(received);
    const total = META.n_orbits_expected || 1;
    progressFill.style.width = (received / total) * 100 + '%';
  });

  es.addEventListener('done', () => {
    es.close();
    finishProgress(true);
    renderOrbitLists();
    renderLegend();
    drawC3Chart();
    applyLayerVisibility();
    updateSelCount();
    const first = orbitObjects.find((o) => o.rec.family === 'coplanar_table1' && o.rec.orbit === 1);
    if (first) selectOrbit(first);
  });

  es.addEventListener('error', () => {
    if (received === 0 && !streamClosed) {
      es.close();
      finishProgress(false, '加载失败：EventSource 连接错误');
    }
  });

  // =====================================================================
  // 场景构建函数
  // =====================================================================
  function makeEarth() {
    const g = new THREE.Group();
    const geo = new THREE.SphereGeometry(R_EARTH_UNITS, 96, 64);
    const fallbackMat = new THREE.MeshPhongMaterial({ color: 0x2a78d2, emissive: 0x0a1830 });
    const mesh = new THREE.Mesh(geo, fallbackMat);
    g.add(mesh);

    // 地球昼夜贴图：直接引用 coast_tool 的 NASA Blue Marble 素材
    // （/coast_tool/assets/earth_day.jpg、earth_night.jpg，同源静态文件）
    const loader = new THREE.TextureLoader();
    let loaded = 0;
    const dayTex = loader.load('/coast_tool/assets/earth_day.jpg', onTex, undefined, onErr);
    const nightTex = loader.load('/coast_tool/assets/earth_night.jpg', onTex, undefined, onErr);
    function onTex() {
      if (++loaded < 2) return;
      [dayTex, nightTex].forEach((t) => {
        t.anisotropy = renderer.capabilities.getMaxAnisotropy();
      });
      // 昼夜项按太阳方向混合（与 coast_tool 同一方案）；本页渲染器为线性
      // 输出，故不做额外 gamma 转换，直接显示贴图原始 sRGB 颜色。
      mesh.material = new THREE.ShaderMaterial({
        uniforms: {
          dayMap: { value: dayTex },
          nightMap: { value: nightTex },
          sunDir: { value: new THREE.Vector3(1, 1, 1).normalize() },
        },
        vertexShader: [
          'varying vec2 vUv;',
          'varying vec3 vNormalL;',
          'void main() {',
          '  vUv = uv;',
          '  vNormalL = normal;',
          '  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);',
          '}',
        ].join('\n'),
        fragmentShader: [
          'uniform sampler2D dayMap, nightMap;',
          'uniform vec3 sunDir;',
          'varying vec2 vUv;',
          'varying vec3 vNormalL;',
          'void main() {',
          '  float d = dot(normalize(vNormalL), normalize(sunDir));',
          '  float k = smoothstep(-0.05, 0.05, d);',
          '  vec3 day = texture2D(dayMap, vUv).rgb * (0.22 + 0.95 * clamp(d, 0.0, 1.0));',
          '  vec3 night = texture2D(nightMap, vUv).rgb * vec3(1.05, 0.95, 0.78);',
          '  gl_FragColor = vec4(mix(night, day, k), 1.0);',
          '}',
        ].join('\n'),
      });
    }
    function onErr() { /* 贴图加载失败：保持纯色地球 */ }

    // 大气辉光（coast_tool 同款：背面加性混合）
    g.add(new THREE.Mesh(
      new THREE.SphereGeometry(R_EARTH_UNITS * 1.045, 64, 48),
      new THREE.ShaderMaterial({
        vertexShader: [
          'varying vec3 vN;',
          'void main() {',
          '  vN = normalize(normalMatrix * normal);',
          '  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);',
          '}',
        ].join('\n'),
        fragmentShader: [
          'varying vec3 vN;',
          'void main() {',
          '  float i = pow(clamp(0.66 - dot(vN, vec3(0.0, 0.0, 1.0)), 0.0, 1.0), 5.0);',
          '  gl_FragColor = vec4(0.30, 0.55, 1.0, 1.0) * i * 0.65;',
          '}',
        ].join('\n'),
        side: THREE.BackSide,
        blending: THREE.AdditiveBlending,
        transparent: true,
        depthWrite: false,
      })
    ));
    return g;
  }

  function makeMoon() {
    const g = new THREE.Group();
    g.add(new THREE.Mesh(
      new THREE.SphereGeometry(R_MOON_UNITS, 48, 32),
      new THREE.MeshPhongMaterial({
        map: makeMoonTexture(), color: 0xffffff, emissive: 0x101010,
        specular: 0x222222, shininess: 4,
      })
    ));
    return g;
  }

  function makeMoonTexture() {
    const W = 1024, H = 512;
    const canvas = document.createElement('canvas');
    canvas.width = W; canvas.height = H;
    const ctx = canvas.getContext('2d');
    const base = ctx.createLinearGradient(0, 0, 0, H);
    base.addColorStop(0, '#b8b4ae');
    base.addColorStop(0.5, '#c9c5bf');
    base.addColorStop(1, '#b0aca6');
    ctx.fillStyle = base;
    ctx.fillRect(0, 0, W, H);
    const grain = ctx.createImageData(W, H);
    for (let i = 0; i < grain.data.length; i += 4) {
      const v = 128 + (Math.random() - 0.5) * 24;
      grain.data[i] = v; grain.data[i + 1] = v; grain.data[i + 2] = v;
      grain.data[i + 3] = 255;
    }
    ctx.putImageData(grain, 0, 0);
    const maria = [
      { x: 0.30, y: 0.42, rx: 0.10, ry: 0.07, a: 0.35 },
      { x: 0.38, y: 0.55, rx: 0.08, ry: 0.06, a: 0.30 },
      { x: 0.24, y: 0.60, rx: 0.09, ry: 0.05, a: 0.28 },
      { x: 0.18, y: 0.48, rx: 0.07, ry: 0.05, a: 0.25 },
      { x: 0.48, y: 0.52, rx: 0.06, ry: 0.05, a: 0.22 },
      { x: 0.55, y: 0.40, rx: 0.05, ry: 0.04, a: 0.20 },
      { x: 0.70, y: 0.45, rx: 0.06, ry: 0.05, a: 0.18 },
    ];
    for (const m of maria) {
      const g2 = ctx.createRadialGradient(
        m.x * W, m.y * H, 0, m.x * W, m.y * H, Math.max(m.rx * W, m.ry * H));
      g2.addColorStop(0, 'rgba(70,66,62,' + m.a + ')');
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
    for (let i = 0; i < 90; i++) {
      const cx = Math.random() * W, cy = Math.random() * H;
      const r = 2 + Math.random() * 9;
      const shade = 0.5 + Math.random() * 0.5;
      ctx.beginPath();
      ctx.arc(cx, cy, r, 0, Math.PI * 2);
      ctx.fillStyle = 'rgba(' + (90 * shade | 0) + ',' + (88 * shade | 0)
        + ',' + (84 * shade | 0) + ',0.55)';
      ctx.fill();
      ctx.beginPath();
      ctx.arc(cx, cy, r * 0.75, 0, Math.PI * 2);
      ctx.fillStyle = 'rgba(' + (150 * shade | 0) + ',' + (148 * shade | 0)
        + ',' + (144 * shade | 0) + ',0.35)';
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
    const len = 20;
    g.add(new THREE.ArrowHelper(new THREE.Vector3(1, 0, 0), new THREE.Vector3(), len, 0xff5577, 3, 2));
    g.add(new THREE.ArrowHelper(new THREE.Vector3(0, 1, 0), new THREE.Vector3(), len, 0x55ff77, 3, 2));
    g.add(new THREE.ArrowHelper(new THREE.Vector3(0, 0, 1), new THREE.Vector3(), len, 0x7799ff, 3, 2));
    const sx = makeLabelSprite('→月', new THREE.Color(0xff5577));
    sx.position.set(-len * 1.15, 0, 0); sx.scale.set(8, 4, 1); g.add(sx);
    const sy = makeLabelSprite('→东', new THREE.Color(0x55ff77));
    sy.position.set(0, len * 1.15, 0); sy.scale.set(8, 4, 1); g.add(sy);
    const sz = makeLabelSprite('→北', new THREE.Color(0x7799ff));
    sz.position.set(0, 0, len * 1.15); sz.scale.set(8, 4, 1); g.add(sz);
    return g;
  }

  function makeEarthMoonLine() {
    const g = new THREE.Group();
    const dist = DEFAULT_MOON_DIST_UNITS;
    const mat = new THREE.LineDashedMaterial({
      color: 0x4ec9ff, dashSize: 4, gapSize: 4, opacity: 0.35, transparent: true,
    });
    const pts = [new THREE.Vector3(0, 0, 0), new THREE.Vector3(-dist, 0, 0)];
    const line = new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts), mat);
    line.computeLineDistances();
    g.add(line);
    return g;
  }

  function makeEmptyLine() {
    const g = new THREE.Group();
    const segments = 96;
    const radius = DEFAULT_MOON_DIST_UNITS;
    const pts = [];
    for (let i = 0; i <= segments; i++) {
      const a = (i / segments) * Math.PI * 2;
      pts.push(new THREE.Vector3(Math.cos(a) * radius, 0, Math.sin(a) * radius));
    }
    const mat = new THREE.LineDashedMaterial({
      color: 0x4ec9ff, dashSize: 2, gapSize: 4, opacity: 0.18, transparent: true,
    });
    const line = new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts), mat);
    line.computeLineDistances();
    g.add(line);
    return g;
  }

  function makeGrid() {
    const grid = new THREE.GridHelper(800, 16, 0x335577, 0x223344);
    grid.rotation.x = Math.PI / 2;
    grid.material.opacity = 0.4;
    grid.material.transparent = true;
    return grid;
  }

  function _latLonToXyz(latDeg, lonDeg, R) {
    const lat = latDeg * Math.PI / 180.0;
    const lon = lonDeg * Math.PI / 180.0;
    return new THREE.Vector3(
      R * Math.cos(lat) * Math.cos(lon),
      R * Math.cos(lat) * Math.sin(lon),
      R * Math.sin(lat),
    );
  }

  function _makeGraticule(R, opts) {
    const o = {
      equatorHex: (opts && opts.equatorHex) || '#4ec9ff',
      meridian0Hex: (opts && opts.meridian0Hex) || '#4ec9ff',
      meridianHex: (opts && opts.meridianHex) || '#2a78d2',
      latitudeHex: (opts && opts.latitudeHex) || '#335577',
      nSeg: (opts && opts.nSeg) || 96,
    };
    const g = new THREE.Group();
    function pushLine(points, hex, opacity) {
      const geo = new THREE.BufferGeometry().setFromPoints(points);
      const mat = new THREE.LineBasicMaterial({
        color: new THREE.Color(hex), transparent: true, opacity,
      });
      g.add(new THREE.Line(geo, mat));
    }
    const equator = [];
    for (let i = 0; i <= o.nSeg; i++) {
      equator.push(_latLonToXyz(0, (i / o.nSeg) * 360, R));
    }
    pushLine(equator, o.equatorHex, 0.85);
    const meridian0 = [];
    for (let i = 0; i <= o.nSeg; i++) {
      meridian0.push(_latLonToXyz(-90 + (i / o.nSeg) * 180, 0, R));
    }
    pushLine(meridian0, o.meridian0Hex, 0.85);
    for (const lonDeg of [90, 180, 270]) {
      const m = [];
      for (let i = 0; i <= o.nSeg; i++) {
        m.push(_latLonToXyz(-90 + (i / o.nSeg) * 180, lonDeg, R));
      }
      pushLine(m, o.meridianHex, 0.55);
    }
    for (const latDeg of [30, 60, -30, -60]) {
      const lat = [];
      for (let i = 0; i <= o.nSeg; i++) {
        lat.push(_latLonToXyz(latDeg, (i / o.nSeg) * 360, R));
      }
      pushLine(lat, o.latitudeHex, 0.45);
    }
    return g;
  }

  function makeEarthGraticule() {
    return _makeGraticule(R_EARTH_UNITS * 1.005, { nSeg: 96 });
  }
  function makeMoonGraticule() {
    return _makeGraticule(R_MOON_UNITS * 1.005, { nSeg: 64 });
  }

  function makeMoonAxis() {
    const g = new THREE.Group();
    const len = R_MOON_UNITS * 1.6;
    const axisPts = [new THREE.Vector3(0, 0, -len), new THREE.Vector3(0, 0, len)];
    g.add(new THREE.Line(
      new THREE.BufferGeometry().setFromPoints(axisPts),
      new THREE.LineBasicMaterial({ color: 0xffd166, transparent: true, opacity: 0.9 })
    ));
    const north = new THREE.Mesh(
      new THREE.SphereGeometry(0.25, 12, 12),
      new THREE.MeshBasicMaterial({ color: 0xffd166 })
    );
    north.position.set(0, 0, len);
    g.add(north);
    const nLabel = makeLabelSprite('北', new THREE.Color(0xffd166), 5);
    nLabel.position.set(0, 0, len + 1.4);
    g.add(nLabel);
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

  function makeLabelSprite(text, color, size) {
    const lines = String(text).split('\n');
    const lineH = 34, padX = 14, padY = 10;
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
    lines.forEach((ln, i) => ctx.fillText(ln, canvas.width / 2, padY + lineH * (i + 0.5)));
    const tex = new THREE.CanvasTexture(canvas);
    tex.needsUpdate = true;
    const mat = new THREE.SpriteMaterial({
      map: tex, transparent: true, depthTest: true, depthWrite: false,
    });
    const sprite = new THREE.Sprite(mat);
    const baseW = size || 20;
    sprite.scale.set(baseW, baseW * (canvas.height / canvas.width), 1);
    sprite.userData.baseText = text;
    return sprite;
  }

  // =====================================================================
  // 单条轨道构建
  // =====================================================================
  function familyColor(family, orbit) {
    const fam = FAMILIES[family];
    if (!fam) return '#9090ff';
    return fam.palette[orbit - 1] || fam.color;
  }

  function buildOneOrbit(rec) {
    const colorHex = familyColor(rec.family, rec.orbit);
    const color = new THREE.Color(colorHex);
    const rMep = rec.r_mep.map(([x, y, z]) =>
      new THREE.Vector3(x / KM_PER_UNIT, y / KM_PER_UNIT, z / KM_PER_UNIT));

    const line = new THREE.Line(
      new THREE.BufferGeometry().setFromPoints(rMep),
      new THREE.LineBasicMaterial({ color, linewidth: 1, opacity: 0.85, transparent: true })
    );
    const dots = new THREE.Points(
      new THREE.BufferGeometry().setFromPoints(rMep),
      new THREE.PointsMaterial({ color, size: 1.6, sizeAttenuation: false })
    );

    const familyShort = FAMILIES[rec.family] ? FAMILIES[rec.family].short : '?';
    const orbitLabel = familyShort + rec.orbit;

    // ---- 近地点地表标记（挂在地球图层下） ----
    const perigeeGroup = new THREE.Group();
    const phiRef = (rec.perigee_lat_ref_deg != null) ? rec.perigee_lat_ref_deg : 0.0;
    const phiUse = (Math.abs(phiRef) < 0.5) ? 28.1 : phiRef;   // 表 1 未给纬度 → 西昌
    const lamUse = (rec.perigee_lon_ref_deg != null)
      ? rec.perigee_lon_ref_deg : rec.perigee_lon_sol_deg;
    const perigeePos = _latLonToXyz(phiUse, lamUse, R_EARTH_UNITS * 1.005);
    const perigeeDot = new THREE.Mesh(
      new THREE.SphereGeometry(0.8, 16, 12),
      new THREE.MeshBasicMaterial({ color })
    );
    perigeeDot.position.copy(perigeePos);
    perigeeGroup.add(perigeeDot);
    const perigeeLabel = makeLabelSprite('#' + orbitLabel, color, 3);
    perigeeLabel.position.copy(perigeePos).add(new THREE.Vector3(0, 3, 0));
    perigeeGroup.add(perigeeLabel);
    const isXichang = Math.abs(phiRef) < 0.5;
    const periSprite = makeLabelSprite(
      'λ=' + lamUse.toFixed(2) + '°\nφ=' + phiUse.toFixed(2) + '°'
        + (isXichang ? ' (西昌假设)' : ''),
      new THREE.Color(isXichang ? 0xffd166 : 0xeeeeee), 3);
    periSprite.position.copy(perigeePos).add(new THREE.Vector3(0, -4, 0));
    perigeeGroup.add(periSprite);
    layers.earth.add(perigeeGroup);

    // ---- 近月点地表标记（挂在月球图层下） ----
    const periluneGroup = new THREE.Group();
    const lamM = (rec.perilune_lon_ref_deg != null)
      ? rec.perilune_lon_ref_deg
      : ((rec.perilune_lon_sol_deg != null) ? rec.perilune_lon_sol_deg : 0.0);
    const phiM = (rec.perilune_lat_ref_deg != null)
      ? rec.perilune_lat_ref_deg
      : ((rec.perilune_lat_sol_deg != null) ? rec.perilune_lat_sol_deg : 0.0);
    const perilunePos = _latLonToXyz(phiM, lamM, R_MOON_UNITS * 1.005);
    const periluneDot = new THREE.Mesh(
      new THREE.SphereGeometry(0.45, 12, 12),
      new THREE.MeshBasicMaterial({ color })
    );
    periluneDot.position.copy(perilunePos);
    periluneGroup.add(periluneDot);
    const periluneLabel = makeLabelSprite('#' + orbitLabel, color, 3);
    periluneLabel.position.copy(perilunePos).add(new THREE.Vector3(0, 1.6, 0));
    periluneGroup.add(periluneLabel);
    layers.moon.add(periluneGroup);

    const label = makeLabelSprite('#' + orbitLabel, color, 8);

    const group = new THREE.Group();
    group.add(line);
    group.add(dots);
    group.add(label);
    scene.add(group);

    return {
      key: rec.family + '::' + rec.orbit,
      family: rec.family, orbit: rec.orbit, rec,
      line, dots, label, group,
      markers: { earth: perigeeGroup, moon: periluneGroup },
      colorHex,
      checked: true,
      samples: null,
    };
  }

  // =====================================================================
  // 左栏：轨道列表 / 过滤 / 选择
  // =====================================================================
  const orbitListsEl = document.getElementById('orbit-lists');
  const filterEl = document.getElementById('filter');
  const selCountEl = document.getElementById('sel-count');

  function renderOrbitLists() {
    orbitListsEl.textContent = '';
    for (const fam of FAMILY_ORDER) {
      const meta = FAMILIES[fam];
      const group = document.createElement('div');
      group.className = 'family-group ' + meta.cls;
      group.dataset.family = fam;

      const h3 = document.createElement('h3');
      const title = document.createElement('span');
      title.textContent = meta.label;
      h3.appendChild(title);
      const count = document.createElement('span');
      count.className = 'count';
      count.textContent = '0';
      h3.appendChild(count);
      const famToggle = document.createElement('button');
      famToggle.type = 'button';
      famToggle.className = 'fam-toggle';
      famToggle.textContent = '全选/反选';
      famToggle.addEventListener('click', () => {
        const infos = orbitObjects.filter((o) => o.family === fam);
        const target = !infos.every((o) => o.checked);
        infos.forEach((o) => setOrbitChecked(o, target));
        applyLayerVisibility();
        updateSelCount();
      });
      h3.appendChild(famToggle);
      group.appendChild(h3);

      const list = document.createElement('div');
      list.className = 'family-list';
      for (const info of orbitObjects.filter((o) => o.family === fam)) {
        list.appendChild(buildRow(info));
      }
      group.appendChild(list);
      orbitListsEl.appendChild(group);
    }
    applyFilter();
  }

  function buildRow(info) {
    const rec = info.rec;
    const meta = FAMILIES[info.family];
    const row = document.createElement('div');
    row.className = 'orbit-row';
    row.dataset.key = info.key;
    row.dataset.search = [
      meta.short + info.orbit, info.family, meta.label,
      'T' + rec.T_h, rec.i0_deg, rec.im_ref_deg,
      rec.vp_ref_km_s.toFixed(3),
    ].join(' ').toLowerCase();

    const cb = document.createElement('input');
    cb.type = 'checkbox';
    cb.checked = true;
    cb.addEventListener('change', (e) => {
      info.checked = e.target.checked;
      applyLayerVisibility();
      updateSelCount();
    });
    row.appendChild(cb);

    const sw = document.createElement('span');
    sw.className = 'swatch';
    sw.style.background = info.colorHex;
    row.appendChild(sw);

    const name = document.createElement('span');
    name.className = 'name';
    const sub = document.createElement('span');
    sub.className = 'sub';
    sub.textContent = ' ' + meta.short + info.orbit + ' · ' + rec.T_h + 'h'
      + ' · i₀ ' + rec.i0_deg + '°';
    name.textContent = '#' + meta.short + info.orbit;
    name.appendChild(sub);
    row.appendChild(name);

    const m = document.createElement('span');
    m.className = 'meta';
    m.textContent = rec.vp_ref_km_s.toFixed(3);
    row.appendChild(m);

    row.addEventListener('click', (e) => {
      if (e.target.tagName === 'INPUT') return;
      selectOrbit(info);
    });
    return row;
  }

  function applyFilter() {
    const q = (filterEl.value || '').trim().toLowerCase();
    document.querySelectorAll('.family-group').forEach((g) => {
      let visible = 0;
      g.querySelectorAll('.orbit-row').forEach((r) => {
        const hit = !q || r.dataset.search.indexOf(q) >= 0;
        r.classList.toggle('hidden-by-filter', !hit);
        if (hit) visible += 1;
      });
      g.style.display = visible ? '' : 'none';
      const c = g.querySelector('.count');
      if (c) c.textContent = visible + '';
    });
  }
  filterEl.addEventListener('input', applyFilter);

  function setOrbitChecked(info, on) {
    info.checked = on;
    const row = document.querySelector('.orbit-row[data-key="' + info.key + '"]');
    if (row) {
      const cb = row.querySelector('input');
      if (cb) cb.checked = on;
    }
  }

  function updateSelCount() {
    const n = orbitObjects.filter((o) => o.checked).length;
    selCountEl.textContent = n + ' / ' + orbitObjects.length + ' 显示';
  }

  document.getElementById('select-all').addEventListener('click', () => {
    orbitObjects.forEach((o) => setOrbitChecked(o, true));
    applyLayerVisibility();
    updateSelCount();
  });
  document.getElementById('select-none').addEventListener('click', () => {
    orbitObjects.forEach((o) => setOrbitChecked(o, false));
    applyLayerVisibility();
    updateSelCount();
  });

  const focusBtn = document.getElementById('focus-toggle');
  focusBtn.classList.toggle('on', focusMode);
  focusBtn.addEventListener('click', () => {
    focusMode = !focusMode;
    focusBtn.classList.toggle('on', focusMode);
    applyLayerVisibility();
  });

  // =====================================================================
  // 选择与详情
  // =====================================================================
  const metaBody = document.getElementById('meta-body');

  function fmt(v, d) {
    return (typeof v === 'number' && isFinite(v)) ? v.toFixed(d == null ? 2 : d) : '—';
  }

  function fmtUtcDate(ms) {
    const d = new Date(ms);
    const p = (n) => String(n).padStart(2, '0');
    return d.getUTCFullYear() + '-' + p(d.getUTCMonth() + 1) + '-' + p(d.getUTCDate())
      + ' ' + p(d.getUTCHours()) + ':' + p(d.getUTCMinutes()) + ' UTC';
  }

  function showMeta(info) {
    const m = info.rec;
    const meta = FAMILIES[info.family];
    metaBody.textContent = '';
    function add(label, value, cls, title) {
      const div = document.createElement('div');
      div.className = 'info-row';
      if (title) div.title = title;
      const l = document.createElement('span');
      l.className = 'label'; l.textContent = label;
      const v = document.createElement('span');
      v.className = 'value' + (cls ? ' ' + cls : '');
      v.textContent = value;
      div.appendChild(l); div.appendChild(v);
      metaBody.appendChild(div);
    }
    add('族 / 编号', meta.label + ' / #' + meta.short + m.orbit, 'highlight');
    add('转移时间 T', fmt(m.T_h, 0) + ' h');
    add('初始倾角 i₀', fmt(m.i0_deg, 2) + '°');
    add('参考 vp', fmt(m.vp_ref_km_s, 6) + ' km/s', null, '论文参考近地点速度');
    add('目标 i_m(ref)', fmt(m.im_ref_deg, 1) + '°');
    add('求解 p', m.p_sol.map((x) => fmt(x, 4)).join(', '), null, 'p = (v_p, Ω, ω)：近地点速度 / 升交点经度 / 近地点幅角');
    add('求解 q', m.q.map((x) => fmt(x, 6)).join(', '), null, 'q = (H_m, cos f, cos θ)：近月点高度 / 位置速度夹角余弦 / 轨道面法向');
    add('近月倾角 i_m(sol)', fmt(m.i_m_sol_deg, 4) + '°');
    add('近月点速度', fmt(m.v_perilune_km_s, 6) + ' km/s');
    add('C3', fmt(m.C3, 4) + ' km²/s²', 'highlight');
    add('近地点 λ/φ (ref)',
      (m.perigee_lon_ref_deg != null ? fmt(m.perigee_lon_ref_deg, 3) + '°' : '—')
      + ' / ' + (m.perigee_lat_ref_deg != null ? fmt(m.perigee_lat_ref_deg, 3) + '°' : '—'));
    add('近地点 λ/φ (sol)',
      fmt(m.perigee_lon_sol_deg, 3) + '° / ' + fmt(m.perigee_lat_sol_deg, 3) + '°');
    add('近月点 λ/φ (ref)',
      (m.perilune_lon_ref_deg != null ? fmt(m.perilune_lon_ref_deg, 3) + '°' : '—')
      + ' / ' + (m.perilune_lat_ref_deg != null ? fmt(m.perilune_lat_ref_deg, 3) + '°' : '—'));
    add('近月点 λ/φ (sol)',
      fmt(m.perilune_lon_sol_deg, 3) + '° / ' + fmt(m.perilune_lat_sol_deg, 3) + '°');
    const arriveMs = Date.parse(META.t_arrive_utc);
    if (isFinite(arriveMs)) {
      add('发射时刻 (UTC)', fmtUtcDate(arriveMs - m.T_h * 3600e3));
      add('到达时刻 (UTC)', fmtUtcDate(arriveMs));
    }
    add('采样', m.t_s.length + ' 点 / ' + fmt(m.T_h, 0) + ' h');

    const det = document.createElement('details');
    det.className = 'raw';
    const sum = document.createElement('summary');
    sum.textContent = '原始数据（JSON 摘要）';
    det.appendChild(sum);
    const pre = document.createElement('div');
    pre.textContent = JSON.stringify({
      family: m.family, orbit: m.orbit, T_h: m.T_h, i0_deg: m.i0_deg,
      p_sol: m.p_sol, q: m.q, C3: m.C3, v_perilune_km_s: m.v_perilune_km_s,
      i_m_sol_deg: m.i_m_sol_deg,
      r_mep_first_km: m.r_mep[0].map((x) => Math.round(x)),
      r_mep_last_km: m.r_mep[m.r_mep.length - 1].map((x) => Math.round(x)),
      n_samples: m.t_s.length,
    }, null, 1);
    det.appendChild(pre);
    metaBody.appendChild(det);
  }

  function selectOrbit(info) {
    selected = info;
    document.querySelectorAll('.orbit-row').forEach((r) =>
      r.classList.toggle('active', r.dataset.key === info.key));
    showMeta(info);
    buildAnimation(info);
    anim.frac = 1;
    setPlaying(false);
    updateAnimation();
    updateChartHighlight();
    applyLayerVisibility();
    updateTimelineEnabled();
  }

  // =====================================================================
  // 动画与时间轴
  // =====================================================================
  const playBtn = document.getElementById('tl-play');
  const prevBtn = document.getElementById('tl-prev');
  const nextBtn = document.getElementById('tl-next');
  const speedSel = document.getElementById('tl-speed');
  const tlSlider = document.getElementById('tl-slider');
  const tlTime = document.getElementById('tl-time');
  const tlPhase = document.getElementById('tl-phase');

  function updateTimelineEnabled() {
    const on = !!selected;
    playBtn.disabled = !on;
    prevBtn.disabled = !on;
    nextBtn.disabled = !on;
    tlSlider.disabled = !on;
    if (!on) {
      tlTime.textContent = '--';
      tlPhase.textContent = '选择轨道';
      scMarker.visible = false;
      trailLine.visible = false;
      moonTrackSel.visible = false;
    }
  }

  function buildAnimation(info) {
    const rec = info.rec;
    const pts = rec.r_mep.map((v) => new THREE.Vector3(v[0] / KM_PER_UNIT, v[1] / KM_PER_UNIT, v[2] / KM_PER_UNIT));
    trailLine.geometry.dispose();
    trailLine.geometry = new THREE.BufferGeometry().setFromPoints(pts);
    trailLine.material.color.set(info.colorHex);
    trailLine.visible = true;

    const mpts = rec.r_moon_mep.map((v) => new THREE.Vector3(v[0] / KM_PER_UNIT, v[1] / KM_PER_UNIT, v[2] / KM_PER_UNIT));
    moonTrackSel.geometry.dispose();
    moonTrackSel.geometry = new THREE.BufferGeometry().setFromPoints(mpts);
    moonTrackSel.visible = layerVisible.moontrack;

    info.samples = { t_s: rec.t_s, r_mep: rec.r_mep, r_moon_mep: rec.r_moon_mep };
    scMarker.material.color.set(info.colorHex);
    scMarker.visible = true;
  }

  function lerp3(a, b, w) {
    return [
      a[0] + (b[0] - a[0]) * w,
      a[1] + (b[1] - a[1]) * w,
      a[2] + (b[2] - a[2]) * w,
    ];
  }

  function fmtCountdown(sec) {
    const h = sec / 3600;
    if (h >= 48) return (h / 24).toFixed(1) + ' d';
    return h.toFixed(1) + ' h';
  }

  function updateAnimation() {
    if (!selected || !selected.samples) return;
    const S = selected.samples;
    const n = S.t_s.length;
    const f = anim.frac * (n - 1);
    const i0 = Math.max(0, Math.min(n - 2, Math.floor(f)));
    const w = Math.max(0, Math.min(1, f - i0));
    const p = lerp3(S.r_mep[i0], S.r_mep[i0 + 1], w);
    const m = lerp3(S.r_moon_mep[i0], S.r_moon_mep[i0 + 1], w);
    scMarker.position.set(p[0] / KM_PER_UNIT, p[1] / KM_PER_UNIT, p[2] / KM_PER_UNIT);
    layers.moon.position.set(m[0] / KM_PER_UNIT, m[1] / KM_PER_UNIT, m[2] / KM_PER_UNIT);

    const k = Math.max(1, Math.min(n, Math.round(f) + 1));
    trailLine.geometry.setDrawRange(0, k);

    const t = S.t_s[i0] + w * (S.t_s[i0 + 1] - S.t_s[i0]);
    const dtArr = t - S.t_s[n - 1];
    const arriveMs = Date.parse(META.t_arrive_utc);
    if (isFinite(arriveMs)) {
      tlTime.textContent = fmtUtcDate(arriveMs + dtArr * 1000);
    } else {
      tlTime.textContent = 'T' + (dtArr >= 0 ? '+' : '−') + fmtCountdown(Math.abs(dtArr));
    }
    tlPhase.textContent = dtArr >= -1 ? '到达' : 'T−' + fmtCountdown(-dtArr);
    tlSlider.value = String(Math.round(anim.frac * 1000));
  }

  function setPlaying(on) {
    anim.playing = !!on;
    playBtn.textContent = anim.playing ? '⏸' : '▶';
  }

  playBtn.addEventListener('click', () => {
    if (!selected) return;
    if (anim.frac >= 1) anim.frac = 0;
    setPlaying(!anim.playing);
  });
  tlSlider.addEventListener('input', (e) => {
    if (!selected) return;
    setPlaying(false);
    anim.frac = parseInt(e.target.value, 10) / 1000;
    updateAnimation();
  });
  prevBtn.addEventListener('click', () => stepSample(-1));
  nextBtn.addEventListener('click', () => stepSample(1));

  function stepSample(dir) {
    if (!selected || !selected.samples) return;
    const n = selected.samples.t_s.length;
    setPlaying(false);
    anim.frac = Math.max(0, Math.min(1, anim.frac + dir / (n - 1)));
    updateAnimation();
  }

  document.addEventListener('keydown', (e) => {
    const tag = (e.target && e.target.tagName) || '';
    if (tag === 'INPUT' || tag === 'SELECT' || tag === 'TEXTAREA') return;
    if (e.code === 'Space') {
      e.preventDefault();
      if (selected) playBtn.click();
    } else if (e.code === 'ArrowLeft') {
      e.preventDefault();
      stepSample(-1);
    } else if (e.code === 'ArrowRight') {
      e.preventDefault();
      stepSample(1);
    }
  });

  // =====================================================================
  // 图层可见性
  // =====================================================================
  function refreshOrbitVisual(o) {
    const focused = !focusMode || !selected || o === selected;
    o.group.visible = o.checked;
    o.line.material.opacity = focused ? (o === selected ? 1.0 : 0.85) : 0.12;
    o.dots.material.opacity = focused ? 1 : 0.15;
    if (o.label) o.label.visible = o.checked && layerVisible.labels && focused;
    if (o.markers) {
      o.markers.earth.visible = o.checked && layerVisible.earth && layerVisible.endpoints;
      o.markers.moon.visible = o.checked && layerVisible.moon && layerVisible.endpoints;
      [o.markers.earth, o.markers.moon].forEach((grp) => {
        grp.traverse((obj) => {
          if (obj.isSprite) {
            obj.visible = layerVisible.labels && focused && obj.parent.visible;
          }
        });
      });
    }
  }

  function applyLayerVisibility() {
    layers.earth.visible = layerVisible.earth;
    layers.moon.visible = layerVisible.moon;
    layers.axis.visible = layerVisible.axis;
    layers.moonline.visible = layerVisible.moonline;
    layers.grid.visible = layerVisible.grid;
    layers.earthgraticule.visible = layerVisible.earthgraticule;
    layers.moongraticule.visible = layerVisible.moongraticule && layerVisible.moon;
    layers.moonaxis.visible = layerVisible.moonaxis && layerVisible.moon;
    layers.moontrack.visible = layerVisible.moontrack && !selected;
    moonTrackSel.visible = layerVisible.moontrack && !!selected;
    orbitObjects.forEach(refreshOrbitVisual);
  }

  function renderLayerToggles() {
    const box = document.getElementById('layer-toggles');
    box.textContent = '';
    for (const [key, label] of LAYERS) {
      const l = document.createElement('label');
      const cb = document.createElement('input');
      cb.type = 'checkbox';
      cb.checked = !!layerVisible[key];
      cb.addEventListener('change', () => {
        layerVisible[key] = cb.checked;
        applyLayerVisibility();
      });
      l.appendChild(cb);
      l.appendChild(document.createTextNode(label));
      box.appendChild(l);
    }
  }
  renderLayerToggles();

  // =====================================================================
  // 视角
  // =====================================================================
  function setView(kind) {
    const dist = 600;
    if (kind === 'fit') { fitSelected(); return; }
    if (kind === 'from_moon_pov') { camera.position.set(-dist, 0, 0); controls.target.set(0, 0, 0); }
    else if (kind === 'from_earth_pov') { camera.position.set(dist, 0, 0); controls.target.set(0, 0, 0); }
    else if (kind === 'xy') { camera.position.set(0, 0, dist); controls.target.set(0, 0, 0); }
    else if (kind === 'xz') { camera.position.set(0, dist, 0); controls.target.set(0, 0, 0); }
    else if (kind === 'yz') { camera.position.set(dist, 0, 0); controls.target.set(0, 0, 0); }
    else if (kind === 'iso') { camera.position.set(450, 350, 450); controls.target.set(0, 0, 0); }
    else { return; }
    camera.lookAt(controls.target);
    controls.syncFromCamera();
    controls.update();
  }

  function fitSelected() {
    if (!selected) return;
    const pts = selected.rec.r_mep.map((v) =>
      new THREE.Vector3(v[0] / KM_PER_UNIT, v[1] / KM_PER_UNIT, v[2] / KM_PER_UNIT));
    const box = new THREE.Box3().setFromPoints(pts);
    const center = box.getCenter(new THREE.Vector3());
    const radius = Math.max(10, box.getSize(new THREE.Vector3()).length() / 2);
    // 固定等距视角，避免极轨在旧视角下被"边缘朝前"看成一条直线
    const dir = new THREE.Vector3(1, 0.55, 1).normalize();
    controls.target.copy(center);
    camera.position.copy(center).addScaledVector(dir, radius * 2.4);
    camera.lookAt(center);
    controls.syncFromCamera();
    controls.update();
  }

  document.querySelectorAll('#view-row button[data-view]').forEach((btn) => {
    btn.addEventListener('click', () => setView(btn.dataset.view));
  });

  // =====================================================================
  // C3 柱状图
  // =====================================================================
  const chart = document.getElementById('chart-c3');
  let chartBars = [];

  function drawC3Chart() {
    chart.textContent = '';
    chartBars = [];
    if (!ORBITS.length) return;
    const W = 340, H = 132, padL = 30, padR = 8, padT = 10, padB = 24;
    const maxC3 = Math.max.apply(null, ORBITS.map((o) => o.C3)) * 1.08;
    const groupGap = 8;
    const groupW = (W - padL - padR - groupGap * (FAMILY_ORDER.length - 1))
      / FAMILY_ORDER.length;
    const y0 = H - padB;
    const plotH = H - padT - padB;
    const ns = 'http://www.w3.org/2000/svg';

    function el(tag, attrs) {
      const e = document.createElementNS(ns, tag);
      for (const k in attrs) e.setAttribute(k, attrs[k]);
      return e;
    }

    // 基线 + 刻度
    chart.appendChild(el('line', { x1: padL, y1: y0, x2: W - padR, y2: y0, stroke: '#30363d', 'stroke-width': 1 }));
    chart.appendChild(el('text', { x: 2, y: y0 + 3, fill: '#8b949e', 'font-size': 8 })).textContent = '0';
    const maxLabel = el('text', { x: 2, y: padT + 6, fill: '#8b949e', 'font-size': 8 });
    maxLabel.textContent = maxC3.toFixed(1);
    chart.appendChild(maxLabel);
    const axisTitle = el('text', { x: 2, y: 8, fill: '#556677', 'font-size': 7 });
    axisTitle.textContent = 'C3';
    chart.appendChild(axisTitle);

    FAMILY_ORDER.forEach((fam, fi) => {
      const meta = FAMILIES[fam];
      const recs = ORBITS.filter((o) => o.family === fam);
      const gx = padL + fi * (groupW + groupGap);
      const barGap = 2;
      const barW = Math.max(2, (groupW - barGap * (recs.length - 1)) / recs.length);
      recs.forEach((rec, ri) => {
        const h = (rec.C3 / maxC3) * plotH;
        const x = gx + ri * (barW + barGap);
        const bar = el('rect', {
          x: x.toFixed(2), y: (y0 - h).toFixed(2),
          width: barW.toFixed(2), height: h.toFixed(2),
          fill: familyColor(fam, rec.orbit), opacity: 0.85,
          rx: 1, cursor: 'pointer',
        });
        const title = el('title', {});
        title.textContent = '#' + meta.short + rec.orbit + ' · T=' + rec.T_h
          + ' h · C3=' + rec.C3.toFixed(3) + ' km²/s²';
        bar.appendChild(title);
        bar.addEventListener('click', () => {
          const info = orbitObjects.find((o) => o.rec === rec);
          if (info) selectOrbit(info);
        });
        chart.appendChild(bar);
        chartBars.push({ rec, bar });
      });
      const famLabel = el('text', {
        x: (gx + groupW / 2).toFixed(1), y: H - 8,
        fill: meta.color, 'font-size': 9, 'text-anchor': 'middle',
      });
      famLabel.textContent = meta.short + ' (' + recs.length + ')';
      chart.appendChild(famLabel);
    });
    updateChartHighlight();
  }

  function updateChartHighlight() {
    chartBars.forEach(({ rec, bar }) => {
      const on = selected && selected.rec === rec;
      bar.setAttribute('stroke', on ? '#ffffff' : 'none');
      bar.setAttribute('stroke-width', on ? '1.5' : '0');
      bar.setAttribute('opacity', on ? '1' : '0.8');
    });
  }

  function renderLegend() {
    const box = document.getElementById('legend');
    box.textContent = '';
    for (const fam of FAMILY_ORDER) {
      const meta = FAMILIES[fam];
      const recs = ORBITS.filter((o) => o.family === fam);
      const item = document.createElement('div');
      item.className = 'item';
      const sw = document.createElement('span');
      sw.className = 'swatch';
      sw.style.background = meta.color;
      item.appendChild(sw);
      const label = document.createElement('b');
      label.textContent = meta.label;
      item.appendChild(label);
      const n = document.createElement('span');
      n.className = 'n';
      const Ts = Array.from(new Set(recs.map((r) => r.T_h))).sort((a, b) => a - b);
      n.textContent = recs.length + ' 条 · T=' + Ts.join('/') + 'h';
      item.appendChild(n);
      box.appendChild(item);
    }
  }

  // =====================================================================
  // 每帧循环
  // =====================================================================
  function tick(ts) {
    requestAnimationFrame(tick);
    const dt = anim.lastTs ? (ts - anim.lastTs) / 1000 : 0;
    anim.lastTs = ts;
    if (anim.playing && selected) {
      const speed = parseFloat(speedSel.value) || 1;
      anim.frac = Math.min(1, anim.frac + dt * speed / ANIM_DURATION_S);
      updateAnimation();
      if (anim.frac >= 1) setPlaying(false);
    }
    if (!selected) {
      layers.moon.position.set(-DEFAULT_MOON_DIST_UNITS, 0, 0);
    }
    controls.update();
    renderer.render(scene, camera);
    updateLabels();
  }

  function updateLabels() {
    orbitObjects.forEach((o) => {
      if (!o.group.visible || !o.label.visible) return;
      const mid = o.rec.r_mep[Math.floor(o.rec.r_mep.length / 2)];
      o.label.position.set(mid[0] / KM_PER_UNIT, mid[1] / KM_PER_UNIT, mid[2] / KM_PER_UNIT);
    });
  }

  function onResize() {
    const w = stage.clientWidth || 1;
    const h = stage.clientHeight || 1;
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
    renderer.setSize(w, h);
  }
  window.addEventListener('resize', onResize);
  onResize();
  // 初始：等距视角框住整个地月系统（地心—月心连线的中点为目标）
  controls.target.set(-DEFAULT_MOON_DIST_UNITS / 2, 0, 0);
  camera.position.copy(controls.target).add(new THREE.Vector3(470, 360, 470));
  camera.lookAt(controls.target);
  controls.syncFromCamera();
  tick(0);

  // =====================================================================
  // 内置极简 OrbitControls（左键旋转 / 右键平移 / 滚轮缩放）
  // =====================================================================
  function makeOrbitControls(cam, dom) {
    const target = new THREE.Vector3(0, 0, 0);
    const spherical = new THREE.Spherical();
    spherical.setFromVector3(cam.position.clone().sub(target));
    let dragging = false, panning = false, lastX = 0, lastY = 0;

    function updateCam() {
      const v = new THREE.Vector3().setFromSpherical(spherical);
      cam.position.copy(target.clone().add(v));
      cam.lookAt(target);
    }
    function onDown(e) {
      dragging = true;
      panning = (e.button === 2 || e.shiftKey);
      lastX = e.clientX; lastY = e.clientY;
    }
    function onMove(e) {
      if (!dragging) return;
      const dx = e.clientX - lastX, dy = e.clientY - lastY;
      lastX = e.clientX; lastY = e.clientY;
      if (panning) {
        const panSpeed = spherical.radius * 0.0015;
        const right = new THREE.Vector3()
          .crossVectors(cam.position.clone().sub(target), new THREE.Vector3(0, 0, 1)).normalize();
        target.addScaledVector(right, -dx * panSpeed);
        target.addScaledVector(new THREE.Vector3(0, 0, 1), dy * panSpeed);
      } else {
        spherical.theta -= dx * 0.005;
        spherical.phi -= dy * 0.005;
        spherical.phi = Math.max(0.05, Math.min(Math.PI - 0.05, spherical.phi));
      }
      updateCam();
    }
    function onUp() { dragging = false; panning = false; }
    function onWheel(e) {
      e.preventDefault();
      const factor = Math.exp(e.deltaY * 0.001);
      spherical.radius = Math.max(1.5, Math.min(3000, spherical.radius * factor));
      updateCam();
    }
    function syncFromCamera() {
      spherical.setFromVector3(cam.position.clone().sub(target));
    }

    dom.addEventListener('mousedown', onDown);
    window.addEventListener('mousemove', onMove);
    window.addEventListener('mouseup', onUp);
    dom.addEventListener('wheel', onWheel, { passive: false });
    dom.addEventListener('contextmenu', (e) => e.preventDefault());

    return {
      get target() { return target; },
      set target(v) { target.copy(v); updateCam(); },
      update() { /* 事件驱动 */ },
      syncFromCamera,
    };
  }

  // =====================================================================
  // 访问统计上报（与 viewer/static/app.js 同口径）
  // =====================================================================
  (function () {
    const url = '/api/visit?path=' + encodeURIComponent(location.pathname);
    if (navigator.sendBeacon) navigator.sendBeacon(url);
    else fetch(url, { method: 'POST', keepalive: true }).catch(function () {});
  })();
})();
