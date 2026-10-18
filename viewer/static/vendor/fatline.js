/**
 * FatLine — 屏幕空间等宽粗线（three r147，UMD 全局/CommonJS 双用）
 *
 * 背景：WebGL 规范中 gl.lineWidth 在 macOS / Windows(ANGLE) 等绝大多数平台
 * 只支持 1px，THREE.LineBasicMaterial / LineDashedMaterial 的 linewidth 会被
 * 静默忽略——想要可见的 2px 线条，只能把每段线展开成 2 个三角形的四边形。
 *
 * 实现：每个线段 6 个角点（start± / end± 两三角形），在顶点着色器里把两端点
 * 变换到裁剪空间、取屏幕方向的法线，按 uniform uWidthPx（物理像素）做
 * NDC 偏移；宽随窗口/像素比统一更新，世界缩放下屏宽恒定（与原 GL_LINE 语义
 * 一致，只是更粗）。
 *
 * 兼容契约（供 elements_demo 平滑替换 THREE.Line）：
 *  - geometry.attributes.position 仍是「源点列」（n 个点），外部 setXYZ /
 *    getX / count 语义与 Line 完全一致；每帧改顶点后置 needsUpdate，
 *    onBeforeRender 检测 version 自动重建派生角点（无需显式调用）。
 *  - geometry.setDrawRange(start, vertexCount) 保留 LINE_STRIP 顶点区间语义；
 *    真实 GL drawRange 恒为全量（窗口外线段在着色器中塌缩为零面积三角形），
 *    逻辑区间另存 geometry.userData.logicalRange。
 *  - computeLineDistances()：虚线线长属性重建（LineDashedMaterial 同名方法）。
 *  - material.color 是真 THREE.Color（getHex/set 可用）；Mesh 具备
 *    position/quaternion/visible/scale/scene.remove/dispose 等全部能力。
 *  - 含 logarithmicDepthBuffer 的 logdepth 着色器块（页面开启了对数深度，
 *    不加会与地球/月球写的对数深度错误比较）。
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) {
    const THREE = require('./three.min.js');
    factory(THREE);
    module.exports = THREE;
  } else {
    factory(root.THREE);
  }
})(typeof self !== 'undefined' ? self : this, function (THREE) {
  'use strict';

  /* 所有 FatLine 材质共享的全局宽度/分辨率状态 */
  const FAT_MATERIALS = [];
  let FAT_CSS_PX = 3;                       /* 线宽（CSS 像素），物理像素 = ×pixelRatio */
  let FAT_PIXEL_RATIO = 1;
  /* 当前绘图缓冲分辨率（物理像素）。页面在创建任何 FatLine 之前就调用
     setFatLineResolution，材质必须从该全局值取初值——否则后创建的材质
     uResolution 恒为 (1,1)，宽度偏移 off ∝ 1/uResolution 被算成 ±uWidthPx
     个 NDC（数倍屏宽），线段四边形整体冲出屏幕 */
  let FAT_RES_W = 1, FAT_RES_H = 1;

  function applyWidth(mat) {
    mat.uniforms.uWidthPx.value = mat.userData.fatCss * FAT_PIXEL_RATIO;
  }

  /** 全局：设置绘图缓冲分辨率（物理像素），通常在 renderer.setSize 后调用 */
  THREE.setFatLineResolution = function (w, h) {
    FAT_RES_W = w; FAT_RES_H = h;
    for (const m of FAT_MATERIALS) m.uniforms.uResolution.value.set(w, h);
  };
  /** 全局：设置 CSS→物理像素比（renderer.getPixelRatio()），重设窗口时调用 */
  THREE.setFatLinePixelRatio = function (pr) {
    FAT_PIXEL_RATIO = pr || 1;
    for (const m of FAT_MATERIALS) applyWidth(m);
  };
  /** 全局：调整线宽（CSS 像素），默认 2 */
  THREE.setFatLineCssWidth = function (cssPx) {
    FAT_CSS_PX = cssPx;
    for (const m of FAT_MATERIALS) { m.userData.fatCss = cssPx; applyWidth(m); }
  };

  const VERT = [
    'uniform vec2 uResolution;',
    'uniform float uWidthPx;',
    'uniform float uSegStart;',
    'uniform float uSegEnd;',
    'attribute vec3 aStart;',
    'attribute vec3 aEnd;',
    'attribute float aSide;',     /* -1 / +1：线段法线两侧 */
    'attribute float aCorner;',   /* 0=start 端，1=end 端 */
    'attribute float aSeg;',      /* 线段序号（揭示窗口判定用） */
    'attribute vec3 aCol0;',
    'attribute vec3 aCol1;',
    'attribute float aDist0;',
    'attribute float aDist1;',
    'varying vec3 vCol;',
    'varying float vDist;',
    '#include <common>',
    '#include <logdepthbuf_pars_vertex>',
    'void main() {',
    '  vec4 c0 = projectionMatrix * modelViewMatrix * vec4(aStart, 1.0);',
    '  vec4 c1 = projectionMatrix * modelViewMatrix * vec4(aEnd, 1.0);',
    '  vec2 px0 = c0.xy / c0.w * (uResolution * 0.5);',
    '  vec2 px1 = c1.xy / c1.w * (uResolution * 0.5);',
    '  vec2 dd = px1 - px0;',
    '  float l2 = dot(dd, dd);',
    '  vec2 dn = l2 > 1e-12 ? dd * inversesqrt(l2) : vec2(1.0, 0.0);',
    '  vec2 nrm = vec2(-dn.y, dn.x);',
    '  float corner = clamp(aCorner, 0.0, 1.0);',
    '  vec4 base = mix(c0, c1, corner);',
    '  vec2 off = nrm * (aSide * uWidthPx * 0.5) / (uResolution * 0.5);',
    '  gl_Position = vec4(base.xy + off * base.w, base.z, base.w);',
    '  /* 揭示窗口外：四角全部塌缩回端点 → 两个三角形面积为 0，不产生片元 */',
    '  if (aSeg < uSegStart || aSeg >= uSegEnd) gl_Position = base;',
    '  vCol = mix(aCol0, aCol1, corner);',
    '  vDist = mix(aDist0, aDist1, corner);',
    '  #include <logdepthbuf_vertex>',
    '}',
  ].join('\n');

  const FRAG = [
    'uniform vec3 uColor;',
    'uniform float uOpacity;',
    'uniform float uUseVColor;',
    'uniform float uDashed;',
    'uniform float uDashSize;',
    'uniform float uGapSize;',
    'uniform float uDashScale;',
    'varying vec3 vCol;',
    'varying float vDist;',
    '#include <common>',
    /* 片元必须用 pars_fragment：logdepthbuf_fragment 引用的 uniform logDepthBufFC
       只在 pars_fragment 声明（pars_vertex 在 USE_LOGDEPTHBUF_EXT 路径下只声明
       varying）。页面 logarithmicDepthBuffer:true + WebGL2 → 写错成 pars_vertex
       时片元编译失败，所有 FatLine 整体隐身（2026-10-16 轨迹消失根因） */
    '#include <logdepthbuf_pars_fragment>',
    'void main() {',
    /* 虚线语义对齐 LineDashedMaterial：mod(lineDist * dashScale, size+gap) <= size */
    '  if (uDashed > 0.5) {',
    '    float ph = mod(vDist * uDashScale, uDashSize + uGapSize);',
    '    if (ph > uDashSize) discard;',
    '  }',
    '  vec3 rgb = uUseVColor > 0.5 ? uColor * vCol : uColor;',
    '  gl_FragColor = vec4(rgb, uOpacity);',
    '  #include <logdepthbuf_fragment>',
    '}',
  ].join('\n');

  function makeMaterial(opts) {
    const mat = new THREE.ShaderMaterial({
      uniforms: {
        uResolution: { value: new THREE.Vector2(FAT_RES_W, FAT_RES_H) },
        uWidthPx: { value: FAT_CSS_PX * FAT_PIXEL_RATIO },
        uColor: { value: new THREE.Color(opts.color != null ? opts.color : 0xffffff) },
        uOpacity: { value: opts.opacity != null ? opts.opacity : 1 },
        uUseVColor: { value: opts.vertexColors ? 1 : 0 },
        uDashed: { value: opts.dashed ? 1 : 0 },
        uDashSize: { value: opts.dashSize != null ? opts.dashSize : 1 },
        uGapSize: { value: opts.gapSize != null ? opts.gapSize : 1 },
        uDashScale: { value: opts.dashScale != null ? opts.dashScale : 1 },
        uSegStart: { value: 0 },
        uSegEnd: { value: 1e9 },
      },
      vertexShader: VERT,
      fragmentShader: FRAG,
      side: THREE.DoubleSide,
      transparent: true,
      depthWrite: opts.depthWrite !== false,
    });
    /* material.color 直连 uniform（getHex 等按 THREE.Color 工作） */
    mat.color = mat.uniforms.uColor.value;
    mat.userData.fatCss = opts.widthCss != null ? opts.widthCss : FAT_CSS_PX;
    applyWidth(mat);
    FAT_MATERIALS.push(mat);
    return mat;
  }

  class FatLine extends THREE.Mesh {
    /**
     * @param source 可选，下列之一：
     *   THREE.Vector3[]（折线点列）/ Float32Array(n*3) /
     *   { position: BufferAttribute|Float32Array, color?: BufferAttribute|Float32Array } /
     *   THREE.BufferGeometry（取其 position/color；topology 同其他源默认 'strip'，
     *   需要线段对拓扑时显式传 opts.topology='segments'）
     * @param opts { color, opacity, vertexColors, dashed, dashSize, gapSize,
     *   dashScale, topology: 'strip'|'segments', widthCss }
     */
    constructor(source, opts) {
      opts = opts || {};
      const geometry = new THREE.BufferGeometry();
      const material = makeMaterial(opts);
      super(geometry, material);
      this.type = 'FatLine';
      this.isFatLine = true;
      this._topology = opts.topology || 'strip';
      this._posVersion = -1;
      this._colVersion = -1;

      /* 逻辑揭示区间（顶点语义，同 LINE_STRIP 的 setDrawRange）。
         真实 geometry.drawRange 必须保持全量：渲染的是派生角点缓冲，
         窗口裁剪在着色器按 aSeg 完成。 */
      geometry.userData.logicalRange = { start: 0, count: Infinity };
      const self = this;
      geometry.setDrawRange = function (start, count) {
        geometry.userData.logicalRange = { start, count };
        if (self._topology === 'segments') {
          material.uniforms.uSegStart.value = Math.floor(start / 2);
          material.uniforms.uSegEnd.value = Math.floor(start / 2) + Math.max(0, Math.floor(count / 2));
        } else {
          const nSeg = count > 0 ? count - 1 : 0;   /* c 个点连 c-1 段 */
          material.uniforms.uSegStart.value = start;
          material.uniforms.uSegEnd.value = start + nSeg;
        }
        geometry.drawRange.start = 0;
        geometry.drawRange.count = Infinity;
      };

      this.onBeforeRender = function () { self.syncDerived(); };

      if (source) this.setSource(source);
    }

    get pointCount() {
      const p = this.geometry.attributes.position;
      return p ? p.count : 0;
    }

    /** 等价语义的源点列属性（setXYZ/getX/count 均为源点） */
    get srcPosition() { return this.geometry.attributes.position; }

    setSource(source, topology) {
      if (topology) this._topology = topology;
      const g = this.geometry;
      let posAttr = null;
      let colAttr = null;

      if (Array.isArray(source)) {
        const arr = new Float32Array(source.length * 3);
        for (let k = 0; k < source.length; k++) {
          arr[k * 3] = source[k].x; arr[k * 3 + 1] = source[k].y; arr[k * 3 + 2] = source[k].z;
        }
        posAttr = new THREE.BufferAttribute(arr, 3);
      } else if (source instanceof Float32Array) {
        posAttr = new THREE.BufferAttribute(source, 3);
      } else if (source && source.isBufferAttribute === true) {
        posAttr = source;
      } else if (source instanceof THREE.BufferGeometry) {
        posAttr = source.attributes.position;
        colAttr = source.attributes.color || null;
        /* 不覆盖拓扑：页面用 setFromPoints 构造的线全部继承 THREE.Line 的
           LINE_STRIP 语义（含 setDrawRange 顶点区间揭示），强制 'segments'
           会丢一半连接段并把揭示窗口折半（2026-10-16 轨迹断裂根因之一）。
           segments 拓扑由调用方显式 opts.topology='segments' 指定（如粗网格） */
      } else if (source && (source.position || source.positions)) {
        const p = source.position || source.positions;
        const c = source.color || source.colors;
        posAttr = p && p.isBufferAttribute === true ? p : new THREE.BufferAttribute(p, 3);
        if (c) colAttr = c.isBufferAttribute === true ? c : new THREE.BufferAttribute(c, 3);
      } else {
        throw new Error('FatLine.setSource: unsupported source');
      }

      g.setAttribute('position', posAttr);
      if (colAttr) g.setAttribute('color', colAttr); else g.deleteAttribute('color');
      this._posVersion = -1;
      this._colVersion = -1;
      this.syncDerived(true);
      /* 初始全显（页面随后会按需 setDrawRange 隐藏/揭示） */
      g.setDrawRange(0, this._topology === 'segments' ? posAttr.count : posAttr.count);
      return this;
    }

    /** 从 Vector3 点列重建（便捷方法，语义同 Line 的 setFromPoints） */
    setPoints(points, topology) {
      this.setSource(points, topology || 'strip');
      return this;
    }

    /** LineDashedMaterial 同名方法：源点变更后重建线长（FatLine 会自动重建） */
    computeLineDistances() {
      this.syncDerived(true);
      return this;
    }

    /** 检测源点/顶点色 version，必要时重建派生角点缓冲 */
    syncDerived(force) {
      const g = this.geometry;
      const pos = g.attributes.position;
      if (!pos) return;
      const col = g.attributes.color;
      const posChanged = force === true || pos.version !== this._posVersion;
      const colChanged = col && (force === true || col.version !== this._colVersion);
      if (!posChanged && !colChanged) return;

      const topology = this._topology;
      const n = pos.count;
      const nSeg = topology === 'segments' ? n >> 1 : Math.max(0, n - 1);
      const nv = nSeg * 6;

      const pa = pos.array;
      const ca = col ? col.array : null;

      const ensure = (name, itemSize) => {
        let attr = g.getAttribute(name);
        if (!attr || attr.count !== nv) {
          attr = new THREE.BufferAttribute(new Float32Array(nv * itemSize), itemSize);
          attr.setUsage(THREE.DynamicDrawUsage);
          g.setAttribute(name, attr);
        }
        return attr;
      };
      const aStart = ensure('aStart', 3).array;
      const aEnd = ensure('aEnd', 3).array;
      const aSide = ensure('aSide', 1).array;
      const aCorner = ensure('aCorner', 1).array;
      const aSeg = ensure('aSeg', 1).array;
      const aCol0 = ensure('aCol0', 3).array;
      const aCol1 = ensure('aCol1', 3).array;
      const aDist0 = ensure('aDist0', 1).array;
      const aDist1 = ensure('aDist1', 1).array;

      /* 每段两三角形：
         (A,-)(B,-)(A,+) / (B,-)(B,+)(A,+) */
      let run = 0;
      for (let s = 0; s < nSeg; s++) {
        const i0 = topology === 'segments' ? s * 2 : s;
        const i1 = topology === 'segments' ? s * 2 + 1 : s + 1;
        const x0 = pa[i0 * 3], y0 = pa[i0 * 3 + 1], z0 = pa[i0 * 3 + 2];
        const x1 = pa[i1 * 3], y1 = pa[i1 * 3 + 1], z1 = pa[i1 * 3 + 2];
        const dseg = Math.sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) + (z1 - z0) * (z1 - z0));
        const d0 = topology === 'segments' ? 0 : run;
        const d1 = topology === 'segments' ? dseg : run + dseg;
        if (topology === 'strip') run += dseg;
        let cr0 = 1, cg0 = 1, cb0 = 1, cr1 = 1, cg1 = 1, cb1 = 1;
        if (ca) {
          cr0 = ca[i0 * 3]; cg0 = ca[i0 * 3 + 1]; cb0 = ca[i0 * 3 + 2];
          cr1 = ca[i1 * 3]; cg1 = ca[i1 * 3 + 1]; cb1 = ca[i1 * 3 + 2];
        }
        const b = s * 6;
        for (let q = 0; q < 6; q++) {
          const j = b + q;
          const isEnd = (q === 1 || q === 3 || q === 4) ? 1 : 0;
          aSide[j] = (q === 0 || q === 1 || q === 3) ? -1 : 1;
          aCorner[j] = isEnd;
          aSeg[j] = s;
          aStart[j * 3] = x0; aStart[j * 3 + 1] = y0; aStart[j * 3 + 2] = z0;
          aEnd[j * 3] = x1; aEnd[j * 3 + 1] = y1; aEnd[j * 3 + 2] = z1;
          aCol0[j * 3] = cr0; aCol0[j * 3 + 1] = cg0; aCol0[j * 3 + 2] = cb0;
          aCol1[j * 3] = cr1; aCol1[j * 3 + 1] = cg1; aCol1[j * 3 + 2] = cb1;
          aDist0[j] = d0; aDist1[j] = d1;
        }
      }
      for (const name of ['aStart', 'aEnd', 'aSide', 'aCorner', 'aSeg', 'aCol0', 'aCol1', 'aDist0', 'aDist1']) {
        g.getAttribute(name).needsUpdate = true;
      }

      /* 索引：恒等索引，让渲染器按派生角点数（而非源点 position.count）发起绘制 */
      if (!g.index || g.index.count !== nv) {
        const idx = nv > 65535 ? new Uint32Array(nv) : new Uint16Array(nv);
        for (let k = 0; k < nv; k++) idx[k] = k;
        g.setIndex(new THREE.BufferAttribute(idx, 1));
      }

      this._posVersion = pos.version;
      this._colVersion = col ? col.version : -1;
      if (g.boundingSphere !== null) g.boundingSphere = null;   /* 源点变化后重算剔除球 */
    }
  }

  THREE.FatLine = FatLine;
});
