/* AgentMoons — 交互式地平线掩码可视化 前端逻辑
 *
 * - 左侧 DEM 地图（Canvas 2D，极射投影坐标换算）
 * - 右侧极坐标地平线图（北=0° 顶部，顺时针，Canvas 2D）
 * - 底部时间 slider：从预加载日序列本地插值，无需请求后端
 */
"use strict";

/* ═══════════════════════ DOM 引用 ═══════════════════════ */
const mapCanvas = document.getElementById("mapCanvas");
const polarCanvas = document.getElementById("polarCanvas");
const mapCtx = mapCanvas.getContext("2d");
const polarCtx = polarCanvas.getContext("2d");
const latInput = document.getElementById("latInput");
const lonInput = document.getElementById("lonInput");
const confirmBtn = document.getElementById("confirmBtn");
const datePicker = document.getElementById("datePicker");
const timeSlider = document.getElementById("timeSlider");
const timeLabel = document.getElementById("timeLabel");
const siteLabel = document.getElementById("siteLabel");
const mapSpinner = document.getElementById("mapSpinner");
const polarSpinner = document.getElementById("polarSpinner");
const historyList = document.getElementById("historyList");
const statsBtn = document.getElementById("statsBtn");
const statsOverlay = document.getElementById("statsOverlay");
const statsCloseBtn = document.getElementById("statsCloseBtn");
const statsRefreshBtn = document.getElementById("statsRefreshBtn");
const statsTableBody = document.getElementById("statsTableBody");
const statsTotalVisits = document.getElementById("statsTotalVisits");
const statsTotalIps = document.getElementById("statsTotalIps");
const statsTotalClicks = document.getElementById("statsTotalClicks");

/* ═══════════════════════ 常量 ═══════════════════════ */
const MOON_R = 1737400.0;        // 月球赤道半径（米），与 utils 一致
const MAP_SIZE_PX = 800;         // 地图 Canvas 尺寸（px，与后端 TILE_PX 一致）
const MAP_HALF_M = 10000;        // 半幅（米）→ 20 km × 20 km
const MAP_SCALE_BAR_M = 5000;    // 比例尺长度（米）→ 5 km
const DAILY_N_POINTS = 288;      // 日序列点数（每 5 分钟）
const DAILY_INTERVAL_MIN = 5;

/* ═══════════════════════ 全局状态 ═══════════════════════ */
const state = {
    // 地图元数据（从 /api/dem_tile 响应头获取）
    mapMeta: null,
    // 已加载的光照底图（BitmapImage，刷新观测点标记时重绘用）
    mapTileImg: null,
    // 观测点（当前站点）
    lat: -89.67,
    lon: 129.78,
    // 地图中心：默认 Shackleton 坑底（与默认观测点一致 → 坑底位于地图正中）
    centerLat: -89.67,
    centerLon: 129.78,
    // 历史观测点列表（用户确认过的站点）
    history: [],
    // 地平线掩码数据
    horizon: null,
    // 日序列数据
    daily: { earth: null, relay: null },
    // 加载中的 flag（防抖）
    loading: false,
};

/* ═══════════════════════ 极射投影坐标换算 ═══════════════════════
 * 与 analysis.elevation_video.utils.lonlat_to_polar_stereo_xy 完全一致：
 *   colat = radians(90 + lat), r = 2R·tan(colat/2)
 *   x = r·sin(lon), y = r·cos(lon)
 */
function lonLatToStereo(lonDeg, latDeg) {
    const colat = ((90 + latDeg) * Math.PI) / 180;
    const r = 2 * MOON_R * Math.tan(colat / 2);
    const lonRad = (lonDeg * Math.PI) / 180;
    return { x: r * Math.sin(lonRad), y: r * Math.cos(lonRad) };
}

function stereoToLonLat(x, y) {
    const r = Math.hypot(x, y);
    const colat = 2 * Math.atan2(r, 2 * MOON_R);
    const lat = (colat * 180) / Math.PI - 90;
    const lon = (Math.atan2(x, y) * 180) / Math.PI;
    return { lon, lat };
}

/* 像素 → 极射投影（米）：y 轴翻转（Canvas y 向下） */
function pixelToStereo(px, py, meta) {
    const x = meta.center_x_stereo_m + (px - meta.img_width_px / 2) * meta.pixel_size_m;
    const y = meta.center_y_stereo_m - (py - meta.img_height_px / 2) * meta.pixel_size_m;
    return { x, y };
}

/* 极射投影（米）→ 像素 */
function stereoToPixel(x, y, meta) {
    const px = meta.img_width_px / 2 + (x - meta.center_x_stereo_m) / meta.pixel_size_m;
    const py = meta.img_height_px / 2 - (y - meta.center_y_stereo_m) / meta.pixel_size_m;
    return { px, py };
}

/* 方位角 → Canvas 弧度（北=0° 顶部，顺时针） */
function azToCanvasAngle(azDeg) {
    return -Math.PI / 2 + (azDeg * Math.PI) / 180;
}

/* ═══════════════════════ 加载状态 ═══════════════════════ */
function setSpinner(el, visible) {
    el.classList.toggle("hidden", !visible);
}

function fmtTimeLabel(minutes) {
    const hh = String(Math.floor(minutes / 60)).padStart(2, "0");
    const mm = String(minutes % 60).padStart(2, "0");
    return `${hh}:${mm} UTC`;
}

/* ═══════════════════════ 左侧：光照 DEM 地图 ═══════════════════════ */
async function loadMapTile() {
    setSpinner(mapSpinner, true);
    try {
        // 光照瓦片随当前日期/时刻变化：date + hour 参与 URL（后端按 MJD 算太阳方向）
        const minutes = parseInt(timeSlider.value, 10) || 0;
        const hour = (minutes / 60).toFixed(4);
        const url =
            `/api/dem_tile?center_lat=${state.centerLat}&center_lon=${state.centerLon}` +
            `&width_m=${MAP_HALF_M * 2}&height_m=${MAP_HALF_M * 2}` +
            `&date=${datePicker.value}&hour=${hour}`;
        const resp = await fetch(url);
        if (!resp.ok) throw new Error(`DEM tile HTTP ${resp.status}`);
        const blob = await resp.blob();
        const meta = {
            center_x_stereo_m: parseFloat(resp.headers.get("center_x_stereo_m")),
            center_y_stereo_m: parseFloat(resp.headers.get("center_y_stereo_m")),
            center_lat_deg: parseFloat(resp.headers.get("center_lat_deg")),
            center_lon_deg: parseFloat(resp.headers.get("center_lon_deg")),
            width_m: parseFloat(resp.headers.get("width_m")),
            height_m: parseFloat(resp.headers.get("height_m")),
            pixel_size_m: parseFloat(resp.headers.get("pixel_size_m")),
            img_width_px: parseInt(resp.headers.get("img_width_px"), 10),
            img_height_px: parseInt(resp.headers.get("img_height_px"), 10),
            vmin_m: parseFloat(resp.headers.get("vmin_m")),
            vmax_m: parseFloat(resp.headers.get("vmax_m")),
            sun_elev_deg: parseFloat(resp.headers.get("sun_elev_deg")),
            sun_az_deg: parseFloat(resp.headers.get("sun_az_deg")),
        };
        state.mapMeta = meta;

        const img = await createImageBitmap(blob);
        state.mapTileImg = img; // 缓存底图：点击新观测点时重绘画布，避免旧标记残留
        mapCtx.clearRect(0, 0, MAP_SIZE_PX, MAP_SIZE_PX);
        mapCtx.drawImage(img, 0, 0, MAP_SIZE_PX, MAP_SIZE_PX);
        drawMapOverlays();
        drawSunIndicator();
    } catch (err) {
        mapCtx.fillStyle = "#fff";
        mapCtx.font = "14px sans-serif";
        mapCtx.textAlign = "center";
        mapCtx.fillText("光照瓦片加载失败: " + err.message, MAP_SIZE_PX / 2, MAP_SIZE_PX / 2);
        console.error(err);
    } finally {
        setSpinner(mapSpinner, false);
    }
}

/* 太阳方向指示（左下角比例尺旁）：箭头指向太阳方位 + 高度角标注。
 * 极射投影下"上=北"（Canvas y 向下 → 北在屏幕上方），方位角 0°=北、
 * 90°=东，与地平线图约定一致；屏幕角 = az - 90°（顺时针从北起算）。
 */
function drawSunIndicator() {
    const meta = state.mapMeta;
    if (!meta || !isFinite(meta.sun_elev_deg) || !isFinite(meta.sun_az_deg)) return;
    const cx = 60, cy = MAP_SIZE_PX - 42;
    const r = 30;
    mapCtx.save();
    mapCtx.translate(cx, cy);

    // 半透明底盘圆
    mapCtx.beginPath();
    mapCtx.arc(0, 0, r + 12, 0, 2 * Math.PI);
    mapCtx.fillStyle = "rgba(0,0,0,0.45)";
    mapCtx.fill();
    mapCtx.strokeStyle = "rgba(255,255,255,0.5)";
    mapCtx.lineWidth = 1;
    mapCtx.stroke();

    // 屏幕角：北在上 → 方位角 az（北=0 顺时针）映射为屏幕角 az-90
    const screenAngle = ((meta.sun_az_deg - 90) * Math.PI) / 180;
    mapCtx.beginPath();
    mapCtx.moveTo(0, 0);
    mapCtx.lineTo(r * Math.cos(screenAngle), r * Math.sin(screenAngle));
    mapCtx.strokeStyle = "#ffd24a";
    mapCtx.lineWidth = 3;
    mapCtx.lineCap = "round";
    mapCtx.stroke();

    // 箭头头
    const tipX = r * Math.cos(screenAngle), tipY = r * Math.sin(screenAngle);
    const back = 8, spread = 0.42;
    mapCtx.beginPath();
    mapCtx.moveTo(tipX, tipY);
    mapCtx.lineTo(
        tipX - back * Math.cos(screenAngle - spread),
        tipY - back * Math.sin(screenAngle - spread)
    );
    mapCtx.lineTo(
        tipX - back * Math.cos(screenAngle + spread),
        tipY - back * Math.sin(screenAngle + spread)
    );
    mapCtx.closePath();
    mapCtx.fillStyle = "#ffd24a";
    mapCtx.fill();

    // 中心圆点
    mapCtx.beginPath();
    mapCtx.arc(0, 0, 3.5, 0, 2 * Math.PI);
    mapCtx.fillStyle = "#fff";
    mapCtx.fill();

    // 标注文字（高度角 / 方位角）
    mapCtx.fillStyle = "#ffd24a";
    mapCtx.font = "bold 11px sans-serif";
    mapCtx.textAlign = "center";
    mapCtx.textBaseline = "top";
    mapCtx.fillText("☀ 太阳", 0, r + 16);
    mapCtx.font = "10px sans-serif";
    mapCtx.fillText(
        `仰角 ${meta.sun_elev_deg.toFixed(1)}° · 方位 ${meta.sun_az_deg.toFixed(0)}°`,
        0, r + 30
    );
    mapCtx.restore();
}

/* 整幅地图刷新：先清空画布并重绘底图（缓存的光照瓦片），再叠加所有图层。
 * 关键：每次点击新观测点时调用本函数，确保旧标记被底图覆盖，画布上始终只有当前一个观测点。
 */
function refreshMapCanvas() {
    if (!state.mapMeta) return;
    mapCtx.clearRect(0, 0, MAP_SIZE_PX, MAP_SIZE_PX);
    if (state.mapTileImg) {
        mapCtx.drawImage(state.mapTileImg, 0, 0, MAP_SIZE_PX, MAP_SIZE_PX);
    }
    drawMapOverlays();
    drawSunIndicator();
}

function drawMapOverlays() {
    if (!state.mapMeta) return;
    const meta = state.mapMeta;
    const px = MAP_SIZE_PX, py = MAP_SIZE_PX;

    // 1. 地图中心十字（白色，位于图中心 = Shackleton 坑底）。
    //    用户要求：历史观测点标记不直接标在地图上，改由左侧历史列表呈现。
    mapCtx.strokeStyle = "rgba(255,255,255,0.9)";
    mapCtx.lineWidth = 1.4;
    const len = 16;
    mapCtx.beginPath();
    mapCtx.moveTo(px / 2 - len, py / 2); mapCtx.lineTo(px / 2 + len, py / 2);
    mapCtx.moveTo(px / 2, py / 2 - len); mapCtx.lineTo(px / 2, py / 2 + len);
    mapCtx.stroke();

    // 2. 比例尺（左下角，5 km）
    const scalePx = MAP_SCALE_BAR_M / meta.pixel_size_m;
    const sx = 14, sy = px - 22;
    mapCtx.strokeStyle = "#fff";
    mapCtx.lineWidth = 2.5;
    mapCtx.beginPath();
    mapCtx.moveTo(sx, sy); mapCtx.lineTo(sx + scalePx, sy);
    mapCtx.stroke();
    mapCtx.fillStyle = "#fff";
    mapCtx.font = "10px sans-serif";
    mapCtx.textAlign = "left";
    mapCtx.textBaseline = "top";
    mapCtx.fillText("5 km", sx, sy + 3);

    // 4. 方向指示器（右上角 N↑）
    mapCtx.fillStyle = "#fff";
    mapCtx.font = "bold 13px sans-serif";
    mapCtx.textAlign = "right";
    mapCtx.textBaseline = "top";
    mapCtx.fillText("N↑", px - 10, 8);

    // 5. 当前观测点标记（最后一次确认/点击的经纬度）。
    //    中心十字表示地图中心（Shackleton 坑底），观测点随输入变化。
    const obs = lonLatToStereo(state.lon, state.lat);
    const { px: opx, py: opy } = stereoToPixel(obs.x, obs.y, meta);
    if (isFinite(opx) && isFinite(opy)) {
        // 外圈脉冲环 + 实心黄点（白色描边，深色 DEM 上醒目）
        mapCtx.beginPath();
        mapCtx.arc(opx, opy, 9, 0, 2 * Math.PI);
        mapCtx.strokeStyle = "rgba(255, 210, 90, 0.55)";
        mapCtx.lineWidth = 1.5;
        mapCtx.stroke();
        mapCtx.beginPath();
        mapCtx.arc(opx, opy, 5, 0, 2 * Math.PI);
        mapCtx.fillStyle = "#ffd24a";
        mapCtx.fill();
        mapCtx.strokeStyle = "rgba(255,255,255,0.9)";
        mapCtx.lineWidth = 1.2;
        mapCtx.stroke();

        // 标签：观测点坐标（紧贴标记右上方）
        mapCtx.fillStyle = "#ffd24a";
        mapCtx.font = "bold 11px sans-serif";
        mapCtx.textAlign = "left";
        mapCtx.textBaseline = "bottom";
        mapCtx.fillText(
            `观测点 ${state.lat.toFixed(2)}°, ${state.lon.toFixed(2)}°`,
            opx + 12, opy - 10
        );
    }
}

/* ═══════════════════════ 历史观测点列表 ═══════════════════════
 * 用户要求：观测点不直接标在地图上，改由左侧历史列表呈现。
 * 每条记录：序号 + 经纬度 + 时间；点击可重新定位，✕ 可删除。
 */
function renderHistory() {
    if (!historyList) return;
    historyList.innerHTML = "";
    if (!state.history.length) {
        const li = document.createElement("li");
        li.className = "empty";
        li.textContent = "暂无历史观测点 — 输入经纬度并点击「确定」后加入";
        historyList.appendChild(li);
        return;
    }
    const activeKey = `${state.lat.toFixed(4)}|${state.lon.toFixed(4)}`;
    state.history.forEach((item, i) => {
        const li = document.createElement("li");
        li.dataset.index = String(i);
        if (item.key === activeKey) li.className = "active";

        const idx = document.createElement("span");
        idx.className = "idx";
        idx.textContent = String(i + 1).padStart(2, "0");

        const coord = document.createElement("span");
        coord.className = "coord";
        coord.textContent = `${item.lat.toFixed(2)}°, ${item.lon.toFixed(2)}°`;

        const clear = document.createElement("span");
        clear.className = "clear";
        clear.textContent = "✕";
        clear.title = "从列表移除";
        clear.addEventListener("click", (e) => {
            e.stopPropagation();
            state.history.splice(i, 1);
            renderHistory();
        });

        li.appendChild(idx);
        li.appendChild(coord);
        li.appendChild(clear);
        li.addEventListener("click", () => {
            selectHistoryItem(item);
        });
        historyList.appendChild(li);
    });
}

/* 点击历史项 → 重新定位该观测点（等价于填写输入框 + 确定） */
function selectHistoryItem(item) {
    latInput.value = item.lat.toFixed(4);
    lonInput.value = item.lon.toFixed(4);
    confirmBtn.click();
}

/* 地图点击：反算经纬度 → 填入输入框 → 触发确定 */
mapCanvas.addEventListener("click", (e) => {
    if (!state.mapMeta) return;
    const rect = mapCanvas.getBoundingClientRect();
    const px = (e.clientX - rect.left) * (MAP_SIZE_PX / rect.width);
    const py = (e.clientY - rect.top) * (MAP_SIZE_PX / rect.height);
    const { x, y } = pixelToStereo(px, py, state.mapMeta);
    const { lon, lat } = stereoToLonLat(x, y);
    latInput.value = lat.toFixed(4);
    lonInput.value = lon.toFixed(4);
    confirmBtn.click();
});

/* ═══════════════════════ 右侧：极坐标图 ═══════════════════════ */
function drawPolarGrid(ctx, cx, cy, R, maxElev) {
    ctx.save();
    // 地平线以下区域着色（0° 圆内侧）：提示负高度角目标位于该区
    ctx.beginPath();
    ctx.arc(cx, cy, R, 0, 2 * Math.PI);
    ctx.fillStyle = "rgba(200, 90, 60, 0.05)";
    ctx.fill();

    // 负高度角刻度圈（-5°, -10°, … 每 -5°，最多画到 -maxElev*0.4）
    const negStep = 5;
    const negMax = Math.max(5, Math.floor(maxElev * 0.4));
    for (let e = -negStep; e >= -negMax; e -= negStep) {
        const r = (Math.abs(e) / maxElev) * R;
        if (r <= 1) continue;
        ctx.beginPath();
        ctx.arc(cx, cy, r, 0, 2 * Math.PI);
        ctx.strokeStyle = "rgba(200, 120, 90, 0.45)";
        ctx.lineWidth = 0.6;
        ctx.setLineDash([3, 4]);
        ctx.stroke();
        ctx.setLineDash([]);
        ctx.fillStyle = "#a06a55";
        ctx.font = "9px sans-serif";
        ctx.textAlign = "left";
        ctx.textBaseline = "middle";
        ctx.fillText(`${e}°`, cx + r + 2, cy - 2);
    }

    // 0° 地平线圆（最粗）
    ctx.beginPath();
    ctx.arc(cx, cy, R, 0, 2 * Math.PI);
    ctx.strokeStyle = "#666";
    ctx.lineWidth = 1.6;
    ctx.stroke();
    ctx.fillStyle = "#666";
    ctx.font = "10px sans-serif";
    ctx.textAlign = "left";
    ctx.textBaseline = "middle";
    ctx.fillText("0°", cx + R + 3, cy - 2);

    // 正高度角同心圆（每 15°）
    for (let e = 15; e <= maxElev; e += 15) {
        const r = (e / maxElev) * R;
        ctx.beginPath();
        ctx.arc(cx, cy, r, 0, 2 * Math.PI);
        ctx.strokeStyle = "#ddd";
        ctx.lineWidth = 0.5;
        ctx.stroke();
        ctx.fillStyle = "#999";
        ctx.font = "10px sans-serif";
        ctx.textAlign = "left";
        ctx.textBaseline = "middle";
        ctx.fillText(`${e}°`, cx + r + 2, cy - 2);
    }

    // 方位角辐射线（每 30°）
    for (let az = 0; az < 360; az += 30) {
        const angle = azToCanvasAngle(az);
        ctx.beginPath();
        ctx.moveTo(cx, cy);
        ctx.lineTo(cx + R * Math.cos(angle), cy + R * Math.sin(angle));
        ctx.strokeStyle = "#ddd";
        ctx.lineWidth = 0.5;
        ctx.stroke();
    }
    // 方位角标签 N/E/S/W
    const labels = [["N", 0], ["E", 90], ["S", 180], ["W", 270]];
    labels.forEach(([label, az]) => {
        const angle = azToCanvasAngle(az);
        const x = cx + (R + 16) * Math.cos(angle);
        const y = cy + (R + 16) * Math.sin(angle);
        ctx.fillStyle = "#333";
        ctx.font = "bold 13px sans-serif";
        ctx.textAlign = "center";
        ctx.textBaseline = "middle";
        ctx.fillText(label, x, y);
    });
    ctx.restore();
}

function drawHorizonMask(ctx, cx, cy, R, maxElev, azArr, elevArr) {
    ctx.save();
    ctx.beginPath();
    const n = Math.min(azArr.length, elevArr.length);
    for (let i = 0; i < n; i++) {
        const r = (Math.max(0, elevArr[i]) / maxElev) * R;
        const angle = azToCanvasAngle(azArr[i]);
        const x = cx + r * Math.cos(angle);
        const y = cy + r * Math.sin(angle);
        if (i === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
    }
    ctx.closePath();
    ctx.fillStyle = "rgba(90, 90, 90, 0.35)";
    ctx.fill();
    ctx.strokeStyle = "#444";
    ctx.lineWidth = 1.2;
    ctx.stroke();
    ctx.restore();
}

function drawEarthBand(ctx, cx, cy, R, maxElev, earth) {
    if (!earth) return;
    const { min, max, current, currentAz } = earth;

    // 全范围带（地球高度角不随方位角变化 → 一圈环）。
    // 高度角可为负（地平线以下，环带画在 0° 圆内侧）。
    // 注意：Canvas arc 半径必须 ≥ 0，负高度角取绝对值画在 0° 圆内侧。
    if (min !== undefined && max !== undefined && isFinite(min) && isFinite(max)) {
        const rMin = (Math.abs(min) / maxElev) * R;
        const rMax = (Math.abs(max) / maxElev) * R;
        const lo = Math.min(rMin, rMax);
        const hi = Math.max(rMin, rMax);
        if (hi - lo > 0.5) {
            ctx.beginPath();
            ctx.arc(cx, cy, hi, 0, 2 * Math.PI);
            ctx.arc(cx, cy, lo, 0, 2 * Math.PI, true);
            ctx.fillStyle = "rgba(30, 120, 220, 0.16)";
            ctx.fill();
            ctx.strokeStyle = "rgba(30, 120, 220, 0.5)";
            ctx.lineWidth = 0.8;
            ctx.stroke();
        }
    }

    // 当前时刻位置点：无论是否在地平线上方都绘制。
    // 负高度角 → 半径取绝对值画在 0° 圆内侧，并用虚线连接中心以标识方位。
    if (current !== undefined && isFinite(current)) {
        const below = current < 0;
        const r = (Math.abs(current) / maxElev) * R;
        const angle = azToCanvasAngle(currentAz || 0);
        const x = cx + r * Math.cos(angle);
        const y = cy + r * Math.sin(angle);

        // 方位指示虚线（中心 → 目标，穿过 0° 圆）
        if (below && r > 2) {
            ctx.save();
            ctx.setLineDash([3, 4]);
            ctx.strokeStyle = "rgba(30, 120, 220, 0.45)";
            ctx.lineWidth = 1;
            ctx.beginPath();
            ctx.moveTo(cx, cy);
            ctx.lineTo(x, y);
            ctx.stroke();
            ctx.restore();
        }

        ctx.beginPath();
        ctx.arc(x, y, below ? 5 : 6, 0, 2 * Math.PI);
        ctx.fillStyle = "#1f77b4";
        ctx.fill();
        ctx.strokeStyle = "white";
        ctx.lineWidth = 1.5;
        ctx.stroke();

        const label = below
            ? `Earth ${current.toFixed(1)}° (below horizon)`
            : `Earth ${current.toFixed(1)}°`;
        ctx.fillStyle = below ? "#7a9cc6" : "#1f77b4";
        ctx.font = "11px sans-serif";
        ctx.textAlign = "left";
        ctx.textBaseline = "bottom";
        ctx.fillText(label, x + 10, y - 8);
    }
}

function drawRelayMarker(ctx, cx, cy, R, maxElev, relay) {
    if (!relay) return;
    // 兼容 API 返回的 snake_case 键（elev_deg/az_deg）与旧式短键（elev/az）
    const elev = relay.elev_deg ?? relay.elev;
    const az = relay.az_deg ?? relay.az;
    const visible = relay.visible ?? true;
    if (!isFinite(elev)) return;

    // 通信状态判定：
    //   中断 = 中继低于 0° 地平线，或被该方位角地形遮挡（低于 horizon mask）
    const horizonAtAz = horizonElevAtAzimuth(az);
    const blockedByTerrain = horizonAtAz !== null && elev < horizonAtAz;
    const interrupted = elev < 0 || !visible || blockedByTerrain;

    // 通信状态配色：正常 = 绿色，中断 = 红色
    const COMM_ON = { fill: "#2e9e5b", ring: "rgba(46, 158, 91, 0.5)", label: "#5fbf85" };
    const COMM_OFF = { fill: "#d62728", ring: "rgba(214, 39, 40, 0.5)", label: "#e07070" };
    const color = interrupted ? COMM_OFF : COMM_ON;

    const below = elev < 0;
    const r = (Math.abs(elev) / maxElev) * R;
    const angle = azToCanvasAngle(az);
    const x = cx + r * Math.cos(angle);
    const y = cy + r * Math.sin(angle);

    // 方位指示虚线（中心 → 目标，穿过 0° 圆；中断时红色、正常时绿色）
    if ((interrupted || r > 2) && r > 0.5) {
        ctx.save();
        ctx.setLineDash([3, 4]);
        ctx.strokeStyle = interrupted ? "rgba(214, 39, 40, 0.5)" : "rgba(46, 158, 91, 0.4)";
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(cx, cy);
        ctx.lineTo(x, y);
        ctx.stroke();
        ctx.restore();
    }

    // 卫星图标：三角形（中断 → 半透明）
    ctx.beginPath();
    ctx.moveTo(x, y - (interrupted ? 5.5 : 7));
    ctx.lineTo(x - 5, y + 4);
    ctx.lineTo(x + 5, y + 4);
    ctx.closePath();
    ctx.fillStyle = interrupted ? color.fill + "99" : color.fill;
    ctx.fill();
    ctx.strokeStyle = "white";
    ctx.lineWidth = 1.5;
    ctx.stroke();

    // 标签：正常/中断 + 高度角 + 遮挡角
    let reason = "";
    if (interrupted) {
        if (blockedByTerrain) reason = `terrain ${horizonAtAz.toFixed(1)}° blocks`;
        else if (below) reason = "below horizon";
        else reason = "not visible";
    }
    const label = interrupted
        ? `Queqiao-2 ${elev.toFixed(1)}° · 通信中断 (${reason})`
        : `Queqiao-2 ${elev.toFixed(1)}° · 通信正常`;
    ctx.fillStyle = color.label;
    ctx.font = "11px sans-serif";
    ctx.textAlign = "left";
    ctx.textBaseline = "top";
    ctx.fillText(label, x + 10, y + 4);
}

/* 查询某方位角（0-360°）的地平线遮挡角；无数据返回 null */
function horizonElevAtAzimuth(azDeg) {
    const h = state.horizon;
    if (!h || !h.azimuths_deg || !h.horizon_elev_deg) return null;
    const n = Math.min(h.azimuths_deg.length, h.horizon_elev_deg.length);
    if (!n) return null;
    const az = ((azDeg % 360) + 360) % 360;
    // 线性插值（方位角近似均匀分布，0°↔360° 环绕）
    const idx = Math.min(n - 1, Math.max(0, Math.round((az / 360) * n)));
    if (idx >= n) return h.horizon_elev_deg[n - 1];
    const i0 = (idx + n) % n, i1 = (idx + 1) % n;
    const a0 = h.azimuths_deg[i0], a1 = a0 + (h.azimuths_deg[i1] - a0 + 360) % 360;
    let t = (az - a0) / (a1 - a0 + 1e-9);
    t = Math.max(0, Math.min(1, t));
    return h.horizon_elev_deg[i0] + (h.horizon_elev_deg[i1] - h.horizon_elev_deg[i0]) * t;
}

/* 自适应量程：根据当前数据（地平线、地球、中继）确定极坐标最大高度角，
 * 保证数据铺满画布而非全部挤在圆心。返回 (maxElev, earthObj)。
 */
function computePolarScale(earthPoint, relayPoint) {
    let peak = 15; // 下限：至少 15° 保证网格可读
    if (state.horizon) {
        peak = Math.max(peak, state.horizon.max_elev_deg || 0);
    }
    // 地球当日范围（含负值，取绝对值）
    if (earthPoint) {
        peak = Math.max(peak, Math.abs(earthPoint.elev_deg));
        const eSeries = state.daily.earth || [];
        for (const p of eSeries) {
            peak = Math.max(peak, Math.abs(p.elev_deg));
        }
    }
    // 中继当前点（含负值，取绝对值）
    if (relayPoint) {
        peak = Math.max(peak, Math.abs(relayPoint.elev_deg));
    }
    // 向上取整到 15° 的倍数，留出余量
    const maxElev = Math.min(90, Math.ceil(peak / 15) * 15);
    return maxElev;
}

/* 综合重绘极坐标图 */
function redrawPolar(earthPoint, relayPoint) {
    const ctx = polarCtx;
    const W = polarCanvas.width, H = polarCanvas.height;
    const cx = W / 2, cy = H / 2;
    const pixelRadius = Math.min(cx, cy) - 42;

    // 自适应量程
    const maxElev = computePolarScale(earthPoint, relayPoint);

    ctx.clearRect(0, 0, W, H);
    // 背景
    ctx.fillStyle = "#f4f6f8";
    ctx.fillRect(0, 0, W, H);

    drawPolarGrid(ctx, cx, cy, pixelRadius, maxElev);

    // 地形遮挡（依赖方位角 → 灰色填充曲线）
    if (state.horizon) {
        drawHorizonMask(
            ctx, cx, cy, pixelRadius, maxElev,
            state.horizon.azimuths_deg, state.horizon.horizon_elev_deg
        );
    }

    // 地球带（范围环 = 当日 min~max + 当前时刻点）
    let earth = null;
    if (earthPoint) {
        const eSeries = state.daily.earth || [];
        let eMin = earthPoint.elev_deg;
        let eMax = earthPoint.elev_deg;
        for (const p of eSeries) {
            if (p.elev_deg < eMin) eMin = p.elev_deg;
            if (p.elev_deg > eMax) eMax = p.elev_deg;
        }
        earth = {
            min: eMin,
            max: eMax,
            current: earthPoint.elev_deg,
            currentAz: earthPoint.az_deg,
        };
    }
    drawEarthBand(ctx, cx, cy, pixelRadius, maxElev, earth);

    // 中继卫星标记（始终显示）
    if (relayPoint) {
        drawRelayMarker(ctx, cx, cy, pixelRadius, maxElev, relayPoint);
    }

    // 站点 + 时间标注
    const minutes = parseInt(timeSlider.value, 10);
    const hh = String(Math.floor(minutes / 60)).padStart(2, "0");
    const mm = String(minutes % 60).padStart(2, "0");
    ctx.fillStyle = "#333";
    ctx.font = "11px sans-serif";
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    ctx.fillText(
        `观测点 ${state.lat.toFixed(2)}°, ${state.lon.toFixed(2)}°  ·  ` +
        `${datePicker.value} ${hh}:${mm} UTC  ·  量程 ±${maxElev}°`,
        cx, H - 22
    );
}

/* ═══════════════════════ 数据加载 ═══════════════════════ */
async function loadHorizonMask(lat, lon) {
    setSpinner(polarSpinner, true);
    try {
        const resp = await fetch(`/api/horizon_mask?lat=${lat}&lon=${lon}`);
        if (!resp.ok) {
            const err = await resp.json();
            throw new Error(err.error || `HTTP ${resp.status}`);
        }
        const data = await resp.json();
        state.horizon = {
            azimuths_deg: data.azimuths_deg,
            horizon_elev_deg: data.horizon_elev_deg,
            max_elev_deg: data.max_elev_deg,
            mean_elev_deg: data.mean_elev_deg,
            elev_m: data.elev_m,
        };
        // 重绘（保持当前时刻标记）
        redrawPolar(getCurrentEarth(), getCurrentRelay());
    } catch (err) {
        console.error("horizon_mask:", err);
    } finally {
        setSpinner(polarSpinner, false);
    }
}

async function loadDailyData(lat, lon, date) {
    try {
        const [earthResp, relayResp] = await Promise.all([
            fetch(`/api/earth_daily?lat=${lat}&lon=${lon}&date=${date}`).then((r) => r.json()),
            fetch(`/api/relay_daily?lat=${lat}&lon=${lon}&date=${date}`).then((r) => r.json()),
        ]);
        state.daily.earth = earthResp.points || [];
        state.daily.relay = relayResp.points || [];
        paintSliderTrack();
        updateCommSummary();
        redrawPolar(getCurrentEarth(), getCurrentRelay());
        updateSkyObjects();
    } catch (err) {
        console.error("daily data:", err);
    }
}

/* 更新通信区间汇总文字（底部 commSummary 元素） */
function updateCommSummary() {
    const el = document.getElementById("commSummary");
    if (!el) return;
    el.textContent = fmtCommSummary();
}

/* 从日序列按当前 slider 分钟数取点（288 点对应 0-1435 分钟，每 5 分钟） */
function getPoint(series, minutes) {
    if (!series || !series.length) return null;
    const idx = Math.round(minutes / DAILY_INTERVAL_MIN);
    return series[Math.max(0, Math.min(idx, series.length - 1))] || null;
}

function getCurrentEarth() {
    return getPoint(state.daily.earth, parseInt(timeSlider.value, 10));
}

function getCurrentRelay() {
    return getPoint(state.daily.relay, parseInt(timeSlider.value, 10));
}

/* ───────────── 通信区间统计 ─────────────
 * 通信可用定义（与 drawRelayMarker 一致）：
 *   中继高度角 > 0° 且未被该方位角地形遮挡（elev ≥ horizon_mask）。
 * 返回：{ intervals: [[startMin, endMin], ...], totalMinutes, coverage }。
 */
function computeCommIntervals() {
    const series = state.daily.relay || [];
    if (!series.length) {
        return { intervals: [], totalMinutes: 0, coverage: 0 };
    }
    const okArr = series.map((p) => {
        const elev = p.elev_deg;
        const blocked = horizonElevAtAzimuth(p.az_deg) ?? -Infinity;
        return elev >= 0 && p.visible !== false && elev >= blocked;
    });
    // 相邻同状态合并成区间（用样本中点分钟数近似）
    const intervals = [];
    let start = null;
    for (let i = 0; i < okArr.length; i++) {
        if (okArr[i] && start === null) start = i;
        else if (!okArr[i] && start !== null) {
            intervals.push([start * DAILY_INTERVAL_MIN, i * DAILY_INTERVAL_MIN]);
            start = null;
        }
    }
    if (start !== null) intervals.push([start * DAILY_INTERVAL_MIN, okArr.length * DAILY_INTERVAL_MIN]);
    const totalMinutes = intervals.reduce((s, [a, b]) => s + (b - a), 0);
    const coverage = (totalMinutes / (series.length * DAILY_INTERVAL_MIN)) * 100;
    return { intervals, totalMinutes, coverage };
}

/* 时间（分钟）→ "HH:MM" */
function fmtMin(min) {
    const m = Math.max(0, Math.min(1439, Math.round(min)));
    return `${String(Math.floor(m / 60)).padStart(2, "0")}:${String(m % 60).padStart(2, "0")}`;
}

/* 汇总文字：着陆点 XX 在 XX 日期总计可通信区间 */
function fmtCommSummary() {
    const { intervals, totalMinutes, coverage } = computeCommIntervals();
    if (!intervals.length) {
        return `着陆点 ${state.lat.toFixed(2)}°, ${state.lon.toFixed(2)}° 在 ${datePicker.value} 无通信窗口`;
    }
    const parts = intervals.map(([a, b]) => `${fmtMin(a)}–${fmtMin(b)} UTC`);
    const totalH = (totalMinutes / 60).toFixed(1);
    return (
        `着陆点 ${state.lat.toFixed(2)}°, ${state.lon.toFixed(2)}° 在 ${datePicker.value} ` +
        `总计可通信 ${totalH} h（${coverage.toFixed(0)}%），区间：${parts.join("，")}`
    );
}

/* 进度条轨道着色：以渐变背景（linear-gradient）表示全天通/断状态。
 * 每 5 分钟一个样本 → 288 段色块；绿色 = 可通信，红色 = 中断。 */
function paintSliderTrack() {
    const slider = timeSlider;
    if (!slider) return;
    const series = state.daily.relay || [];
    if (!series.length) {
        slider.style.background = "";
        return;
    }
    const steps = series.length;
    const stops = [];
    for (let i = 0; i < steps; i++) {
        const p = series[i];
        const elev = p.elev_deg;
        const blocked = horizonElevAtAzimuth(p.az_deg) ?? -Infinity;
        const ok = elev >= 0 && p.visible !== false && elev >= blocked;
        const color = ok ? "#2e9e5b" : "#d62728";
        const pct = (i / steps) * 100;
        stops.push(`${color} ${pct}%`);
        stops.push(`${color} ${((i + 1) / steps) * 100}%`);
    }
    slider.style.background = `linear-gradient(to right, ${stops.join(",")})`;
    slider.style.accentColor = "#ffffff";
}

/* ═══════════════════════ 访问统计（IP 访问量 + 按钮点击量） ═══════════════════════
 * 页面加载时上报一次访问（POST /api/visit），关键按钮点击时上报（POST /api/click）；
 * 服务端按 IP + 日期以 CSV 持久化，/api/stats 按日期聚合供统计面板展示。
 */
// 上报访问（页面加载一次；sendBeacon 可确保页面关闭前送达）
function trackVisit() {
    const url = "/api/visit?path=" + encodeURIComponent(location.pathname);
    if (navigator.sendBeacon) {
        navigator.sendBeacon(url);
    } else {
        fetch(url, { method: "POST", keepalive: true }).catch(() => {});
    }
}

// 上报一次按钮点击（静默失败，不影响主流程）
function trackClick(name) {
    const url = "/api/click?button=" + encodeURIComponent(name);
    if (navigator.sendBeacon) {
        navigator.sendBeacon(url);
    } else {
        fetch(url, { method: "POST", keepalive: true }).catch(() => {});
    }
}

/* 转义 HTML，防止按钮名/路径注入 */
function escapeHtml(s) {
    return String(s).replace(/[&<>"']/g, (c) => ({
        "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;",
    }[c]));
}

/* 加载并按日期渲染统计汇总（访问量 / 唯一 IP / 按钮点击量） */
async function loadStats() {
    try {
        const resp = await fetch("/api/stats", { cache: "no-store" });
        if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
        const data = await resp.json();
        statsTotalVisits.textContent = String(data.total_visits ?? 0);
        statsTotalIps.textContent = String(data.total_unique_ips ?? 0);
        statsTotalClicks.textContent = String(data.total_clicks ?? 0);

        if (!data.days || data.days.length === 0) {
            statsTableBody.innerHTML = '<tr class="empty-row"><td colspan="4">暂无数据</td></tr>';
            return;
        }
        statsTableBody.innerHTML = data.days.map((d) => {
            const clicks = Object.entries(d.clicks || {})
                .sort((a, b) => b[1] - a[1])
                .map(([btn, n]) => `${escapeHtml(btn)} ×${n}`)
                .join("、") || "—";
            return `<tr>
                <td>${escapeHtml(d.date)}</td>
                <td>${d.visits}</td>
                <td>${d.unique_ips}</td>
                <td class="clicks-cell">${clicks}</td>
            </tr>`;
        }).join("");
    } catch {
        statsTotalVisits.textContent = "—";
        statsTableBody.innerHTML = '<tr class="empty-row"><td colspan="4">加载失败</td></tr>';
    }
}

function openStats() {
    statsOverlay.classList.add("open");
    loadStats();
}

function closeStats() {
    statsOverlay.classList.remove("open");
}

/* 统计面板交互：打开 / 关闭 / 刷新 / Esc 退出 */
if (statsBtn) {
    statsBtn.addEventListener("click", () => {
        trackClick("打开统计");
        openStats();
    });
}
if (statsCloseBtn) statsCloseBtn.addEventListener("click", closeStats);
if (statsRefreshBtn) statsRefreshBtn.addEventListener("click", loadStats);
if (statsOverlay) {
    statsOverlay.addEventListener("click", (e) => {
        if (e.target === statsOverlay) closeStats();
    });
}
document.addEventListener("keydown", (e) => {
    if (e.key === "Escape" && statsOverlay && statsOverlay.classList.contains("open")) {
        closeStats();
    }
});

/* ═══════════════════════ 事件绑定 ═══════════════════════ */
confirmBtn.addEventListener("click", async () => {
    trackClick("确定");
    if (state.loading) return;
    const lat = parseFloat(latInput.value);
    const lon = parseFloat(lonInput.value);
    if (!isFinite(lat) || !isFinite(lon) || lat < -90 || lat > -80 || lon < -180 || lon > 180) {
        alert("请输入有效经纬度：纬度 -90 ~ -80，经度 -180 ~ 180");
        return;
    }
    state.loading = true;
    confirmBtn.disabled = true;
    state.lat = lat;
    state.lon = lon;
    siteLabel.textContent = `观测点: ${lat.toFixed(2)}°, ${lon.toFixed(2)}°`;

    // 1. 加入历史列表（去重；历史标记不画在地图上）
    const key = `${lat.toFixed(4)}|${lon.toFixed(4)}`;
    const existing = state.history.findIndex((h) => h.key === key);
    if (existing >= 0) {
        state.history.splice(existing, 1); // 移到列表末尾（最近使用）
    }
    state.history.push({
        key,
        lat,
        lon,
        time: datePicker.value + " " + timeLabel.textContent,
    });
    renderHistory();

    // 2. 地图上更新观测点标记：先清空画布并重绘底图，再叠加新标记，
    //    确保每次点击后画布上只有一个观测点（旧标记被底图覆盖，不再残留）
    refreshMapCanvas();

    // 3. 刷新地平线掩码（仅依赖位置）
    await loadHorizonMask(lat, lon);

    // 4. 刷新日序列（位置 + 日期）
    await loadDailyData(lat, lon, datePicker.value);

    // 5. 观测点变化 → 重新加载 3D 地形（若 3D 页签已激活）
    SKY.terrainKey = null;
    if (tabPanes.skyTab && tabPanes.skyTab.classList.contains("active")) {
        SKY.terrainKey = `${lat.toFixed(4)}|${lon.toFixed(4)}`;
        loadTerrainGLTF(lat, lon);
    }

    state.loading = false;
    confirmBtn.disabled = false;
});

/* 防抖：拖动 slider / 切换日期时延迟重载光照瓦片（避免每帧都请求后端） */
let mapReloadTimer = null;
function scheduleMapReload(delayMs = 250) {
    if (mapReloadTimer) clearTimeout(mapReloadTimer);
    mapReloadTimer = setTimeout(() => {
        mapReloadTimer = null;
        loadMapTile();
    }, delayMs);
}

/* 时间 slider：极坐标图/天空为纯前端插值；光照瓦片防抖重载（250ms 后按新时刻刷新太阳方向） */
timeSlider.addEventListener("input", (e) => {
    const minutes = parseInt(e.target.value, 10);
    timeLabel.textContent = fmtTimeLabel(minutes);
    redrawPolar(getCurrentEarth(), getCurrentRelay());
    // 更新当前时刻通信状态指示（进度条已按全天着色，此处刷新当前点状态）
    updateCommStatusBadge();
    // 同步 3D 天空视角中地球 / 鹊桥位置
    updateSkyObjects();
    scheduleMapReload();
});

/* 更新当前时刻通信状态徽标（commBadge） */
function updateCommStatusBadge() {
    const badge = document.getElementById("commBadge");
    if (!badge) return;
    const relay = getCurrentRelay();
    if (!relay) {
        badge.textContent = "通信状态：—";
        badge.className = "comm-badge idle";
        return;
    }
    const elev = relay.elev_deg;
    const blocked = horizonElevAtAzimuth(relay.az_deg) ?? -Infinity;
    const ok = elev >= 0 && relay.visible !== false && elev >= blocked;
    badge.textContent = ok
        ? `当前 ${fmtTimeLabel(parseInt(timeSlider.value, 10))} · 通信正常`
        : `当前 ${fmtTimeLabel(parseInt(timeSlider.value, 10))} · 通信中断`;
    badge.className = ok ? "comm-badge on" : "comm-badge off";
}

/* 日期改变：重新加载日序列（地平线掩码不变，仅依赖位置）；光照瓦片按新日期重载 */
datePicker.addEventListener("change", (e) => {
    loadDailyData(state.lat, state.lon, e.target.value);
    updateCommStatusBadge();
    scheduleMapReload(0);
});

/* ═══════════════════════ 右侧页签（极坐标 / 3D 天空） ═══════════════════════ */
const tabButtons = document.querySelectorAll(".tab-btn");
const tabPanes = {
    polarTab: document.getElementById("polarTab"),
    skyTab: document.getElementById("skyTab"),
};

function switchTab(name) {
    tabButtons.forEach((b) => b.classList.toggle("active", b.dataset.tab === name));
    Object.entries(tabPanes).forEach(([k, el]) => {
        if (el) el.classList.toggle("active", k === name);
    });
    if (name === "skyTab") {
        initSkyScene();
        resizeSky();
        const key = `${state.lat.toFixed(4)}|${state.lon.toFixed(4)}`;
        if (SKY.terrainKey !== key) {
            SKY.terrainKey = key;
            loadTerrainGLTF(state.lat, state.lon);
        }
    }
}
tabButtons.forEach((b) =>
    b.addEventListener("click", () => {
        trackClick(b.textContent.trim() || b.dataset.tab);
        switchTab(b.dataset.tab);
    })
);

/* ═══════════════════════ 3D 第一人称天空视角（three.js） ═══════════════════════
 * PyVista 在服务端把 DEM 地形导出为 GLTF（月球灰度顶点色），
 * 浏览器用 three.js 渲染：星空穹顶 + 蓝色地球球体 + 红色鹊桥三角，
 * 鼠标拖动 = 第一人称 yaw/pitch 环视（看向天空），滚轮 = 眼高缩放。
 *
 * 坐标约定（PyVista export_gltf 节点矩阵实测，glTF 列主序）：
 *   VTK (X=east, Y=north, Z=up) → glTF 世界 (x=north, y=up, z=east)
 * 故 three.js 世界 = X=北, Y=上, Z=东（右手系），相机 up = +Y。
 */
const skyCanvas = document.getElementById("skyCanvas");
const skySpinner = document.getElementById("skySpinner");
const skyHud = document.getElementById("skyHud");
const hudCtx = skyHud ? skyHud.getContext("2d") : null;
const SKY_UP = new THREE.Vector3(0, 1, 0);      // 世界"上" = +Y
const SKY_UP_FALLBACK = new THREE.Vector3(0, 0, 1); // 看向天顶时备用 up（东）

const SKY = {
    renderer: null,
    scene: null,
    camera: null,
    cameraPos: new THREE.Vector3(0, 2, 0), // 观察点（眼高沿 Y 轴）
    terrain: null,
    terrainKey: null,
    terrainHalf: 15000,
    obsY: 0,        // 观测点 DEM 高程（绝对米，VTK Z → glTF Y）
    eyeH: 2,        // 眼高（米）
    earthMesh: null,
    earthGhost: null,
    earthBeam: null,
    relayMesh: null,
    relayBeam: null,
    horizonRing: null,
    ready: false,
    loadSeq: 0,     // 加载序号：丢弃过期响应（快速连续切换观测点时）
    az: 0,          // 视角方位（0=北，顺时针增加）
    el: 55,         // 视角仰角（度，默认看向天空；过高会与 up 退化）
    drag: null,
};

/* 方位角/仰角 → 世界方向（glTF 帧：X=北, Y=上, Z=东） */
function skyDirFromAzEl(azDeg, elDeg) {
    const az = (azDeg * Math.PI) / 180;
    const el = (elDeg * Math.PI) / 180;
    const cosEl = Math.cos(el);
    // 北分量 → X（cos az）；上分量 → Y；东分量 → Z（sin az）
    return new THREE.Vector3(cosEl * Math.cos(az), Math.sin(el), cosEl * Math.sin(az));
}

/* 柔和圆点星纹（canvas 径向渐变 → PointsMaterial map） */
function makeStarTexture() {
    const c = document.createElement("canvas");
    c.width = c.height = 64;
    const g = c.getContext("2d");
    const grd = g.createRadialGradient(32, 32, 0, 32, 32, 32);
    grd.addColorStop(0, "rgba(255,255,255,1)");
    grd.addColorStop(0.35, "rgba(255,255,255,0.85)");
    grd.addColorStop(1, "rgba(255,255,255,0)");
    g.fillStyle = grd;
    g.fillRect(0, 0, 64, 64);
    return new THREE.CanvasTexture(c);
}

/* 星空穹顶：均匀随机方向的白色星点（远球壳，地平线下被地形遮挡） */
function addStars() {
    const N = 2600;
    const R = 55000;
    const pos = new Float32Array(N * 3);
    const col = new Float32Array(N * 3);
    for (let i = 0; i < N; i++) {
        // 单位球均匀采样（拒绝采样）
        let x, y, z;
        do {
            x = Math.random() * 2 - 1;
            y = Math.random() * 2 - 1;
            z = Math.random() * 2 - 1;
        } while (x * x + y * y + z * z > 1 || x * x + y * y + z * z < 1e-4);
        const inv = R / Math.sqrt(x * x + y * y + z * z);
        pos[i * 3] = x * inv;
        pos[i * 3 + 1] = y * inv;
        pos[i * 3 + 2] = z * inv;
        const b = 0.55 + Math.random() * 0.45; // 亮度差异
        col[i * 3] = b;
        col[i * 3 + 1] = b;
        col[i * 3 + 2] = Math.min(1, b + 0.05);
    }
    const geo = new THREE.BufferGeometry();
    geo.setAttribute("position", new THREE.BufferAttribute(pos, 3));
    geo.setAttribute("color", new THREE.BufferAttribute(col, 3));
    const mat = new THREE.PointsMaterial({
        size: 1200,
        map: makeStarTexture(),
        transparent: true,
        vertexColors: true,
        depthWrite: false,
        sizeAttenuation: true,
    });
    SKY.scene.add(new THREE.Points(geo, mat));
}

function initSkyScene() {
    if (SKY.renderer) return;
    SKY.renderer = new THREE.WebGLRenderer({ canvas: skyCanvas, antialias: true });
    SKY.renderer.setClearColor(0x05070c);
    SKY.renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    SKY.scene = new THREE.Scene();
    SKY.camera = new THREE.PerspectiveCamera(75, 1, 0.5, 120000);
    addStars();

    // 蓝色地球球体（半透明 + 忽略深度 → 即使被地形/地平线遮挡也始终可见）
    SKY.earthMesh = new THREE.Mesh(
        new THREE.SphereGeometry(2000, 24, 16),
        new THREE.MeshBasicMaterial({
            color: 0x2b6be4,
            transparent: true,
            opacity: 0.5,
            depthTest: false,
            side: THREE.DoubleSide,
        })
    );
    SKY.earthMesh.renderOrder = 20;
    SKY.scene.add(SKY.earthMesh);

    // 地球线框"幽灵"球：位于地球真实位置，指示其被遮挡/位于地平线下的状态
    SKY.earthGhost = new THREE.Mesh(
        new THREE.SphereGeometry(2000, 20, 14),
        new THREE.MeshBasicMaterial({
            color: 0x6cb6ff,
            wireframe: true,
            transparent: true,
            opacity: 0.35,
            depthTest: false,
        })
    );
    SKY.earthGhost.renderOrder = 21;
    SKY.scene.add(SKY.earthGhost);

    // 地球视线光束（虚线，观察点 → 地球；深度忽略保证始终可见）
    SKY.earthBeam = new THREE.Line(
        new THREE.BufferGeometry(),
        new THREE.LineDashedMaterial({
            color: 0x4a9cff,
            dashSize: 600,
            gapSize: 350,
            transparent: true,
            opacity: 0.7,
            depthTest: false,
        })
    );
    SKY.earthBeam.renderOrder = 19;
    SKY.scene.add(SKY.earthBeam);

    // 红色鹊桥二号三角锥（radialSegments=3 → 三角四面体）
    SKY.relayMesh = new THREE.Mesh(
        new THREE.ConeGeometry(1400, 3200, 3),
        new THREE.MeshBasicMaterial({ color: 0xe23b2e, side: THREE.DoubleSide })
    );
    SKY.relayMesh.renderOrder = 25;
    SKY.scene.add(SKY.relayMesh);

    // 鹊桥视线光束（虚线，观察点 → 鹊桥）
    SKY.relayBeam = new THREE.Line(
        new THREE.BufferGeometry(),
        new THREE.LineDashedMaterial({
            color: 0xe23b2e,
            dashSize: 500,
            gapSize: 300,
            transparent: true,
            opacity: 0.55,
            depthTest: false,
        })
    );
    SKY.relayBeam.renderOrder = 24;
    SKY.scene.add(SKY.relayBeam);

    // 地平线环（0° 参考）：观察点眼高处的水平圆环，半径超出地形范围，
    // 红色/青色高亮；深度忽略 → 始终清晰可见，明确标出天文地平线。
    const HZ_R = 50000;
    const hzPts = [];
    for (let i = 0; i <= 180; i++) {
        const t = (i / 180) * 2 * Math.PI;
        hzPts.push(new THREE.Vector3(HZ_R * Math.cos(t), 0, HZ_R * Math.sin(t)));
    }
    const hzCurve = new THREE.CatmullRomCurve3(hzPts, true);
    SKY.horizonRing = new THREE.Mesh(
        new THREE.TubeGeometry(hzCurve, 256, 14, 6, true),
        new THREE.MeshBasicMaterial({
            color: 0x4fd2ff,
            transparent: true,
            opacity: 0.9,
            depthTest: false,
        })
    );
    SKY.horizonRing.renderOrder = 18;
    SKY.scene.add(SKY.horizonRing);

    resizeSky();
    SKY.renderer.setAnimationLoop(animateSky);
}

function animateSky() {
    if (SKY.renderer && SKY.scene) SKY.renderer.render(SKY.scene, SKY.camera);
    drawHud();
}

/* ── XYZ 轴方向指示 HUD（左下角叠加） ──
 * 将世界轴（X=北, Y=上, Z=东）经相机旋转投影到屏幕，随视角旋转，
 * 直观展示"当前看向哪个方向"；同时显示 az/el 数值。
 */
function drawHud() {
    if (!hudCtx || !SKY.camera || !SKY.ready) return;
    const W = skyHud.width, H = skyHud.height;
    const cx = W / 2, cy = H / 2;
    hudCtx.clearRect(0, 0, W, H);

    // 外圈参考圆
    hudCtx.beginPath();
    hudCtx.arc(cx, cy, 62, 0, 2 * Math.PI);
    hudCtx.strokeStyle = "rgba(255,255,255,0.18)";
    hudCtx.lineWidth = 1;
    hudCtx.stroke();

    // 世界轴 → 相机空间（保证 matrixWorldInverse 最新）
    SKY.camera.updateMatrixWorld();
    const rot = new THREE.Matrix3().setFromMatrix4(SKY.camera.matrixWorldInverse);
    const LEN = 46;
    const axes = [
        { v: new THREE.Vector3(1, 0, 0), color: "#ff6b6b", label: "X·N" },  // X = 北
        { v: new THREE.Vector3(0, 1, 0), color: "#6bff8b", label: "Y·UP" }, // Y = 上
        { v: new THREE.Vector3(0, 0, 1), color: "#6b9bff", label: "Z·E" },  // Z = 东
    ];
    for (const a of axes) {
        const c = a.v.clone().applyMatrix3(rot); // 相机坐标：x=右, y=上, -z=前
        const front = c.z < 0;                    // 朝向镜头或远离
        const k = front ? 1.0 : 0.45;
        const ex = cx + c.x * LEN * k;
        const ey = cy - c.y * LEN * k;
        // 轴线
        hudCtx.beginPath();
        hudCtx.moveTo(cx, cy);
        hudCtx.lineTo(ex, ey);
        hudCtx.strokeStyle = a.color;
        hudCtx.globalAlpha = front ? 0.95 : 0.4;
        hudCtx.lineWidth = front ? 2.2 : 1.4;
        hudCtx.stroke();
        // 箭头端
        hudCtx.beginPath();
        hudCtx.arc(ex, ey, 3.2, 0, 2 * Math.PI);
        hudCtx.fillStyle = a.color;
        hudCtx.fill();
        // 标签
        hudCtx.fillStyle = a.color;
        hudCtx.font = "bold 10px sans-serif";
        hudCtx.textAlign = "center";
        hudCtx.textBaseline = "middle";
        hudCtx.fillText(a.label, ex, ey + (front ? 13 : -13));
    }
    hudCtx.globalAlpha = 1;

    // az / el 数值（旋转程度）
    hudCtx.fillStyle = "rgba(255,255,255,0.85)";
    hudCtx.font = "10px monospace";
    hudCtx.textAlign = "center";
    hudCtx.textBaseline = "top";
    hudCtx.fillText(`az ${SKY.az.toFixed(0)}°`, cx, cy + 40);
    hudCtx.fillText(`el ${SKY.el.toFixed(0)}°`, cx, cy + 52);
}

function resizeSky() {
    if (!SKY.renderer || !SKY.camera) return;
    const w = skyCanvas.clientWidth || 520;
    const h = skyCanvas.clientHeight || 520;
    SKY.renderer.setSize(w, h, false);
    SKY.camera.aspect = w / h;
    SKY.camera.updateProjectionMatrix();
}

/* 按当前 yaw/pitch 看向天空（第一人称） */
function lookFromAzEl() {
    const dir = skyDirFromAzEl(SKY.az, SKY.el);
    SKY.camera.position.copy(SKY.cameraPos);
    // 看向天顶时 up 会与视线退化，回退到"东"方向
    const up = Math.abs(dir.dot(SKY_UP)) > 0.999 ? SKY_UP_FALLBACK : SKY_UP;
    SKY.camera.up.copy(up);
    SKY.camera.lookAt(
        SKY.cameraPos.x + dir.x,
        SKY.cameraPos.y + dir.y,
        SKY.cameraPos.z + dir.z
    );
}

/* 加载观测点周围 3D 地形（GLTF，月球灰度纹理）。
 * 旧地形保留到新地形解析完成后再替换；失败时保留旧地形，
 * 并用 loadSeq 令牌丢弃过期响应（快速连续切换观测点时）。
 */
async function loadTerrainGLTF(lat, lon) {
    const seq = ++SKY.loadSeq;
    setSpinner(skySpinner, true);
    try {
        const resp = await fetch(`/api/terrain_gltf?lat=${lat}&lon=${lon}`);
        if (!resp.ok) throw new Error(`terrain HTTP ${resp.status}`);
        if (seq !== SKY.loadSeq) return; // 已发起更新的加载
        SKY.obsY = parseFloat(resp.headers.get("obs_z_m")) || 0;
        SKY.terrainHalf = parseFloat(resp.headers.get("terrain_half_m")) || 15000;
        const buf = await resp.arrayBuffer();
        const loader = new THREE.GLTFLoader();
        const gltf = await new Promise((res, rej) => loader.parse(buf, "", res, rej));
        if (seq !== SKY.loadSeq) return;

        const group = gltf.scene || gltf.scenes[0];
        group.traverse((o) => {
            if (o.isMesh) {
                o.material = new THREE.MeshBasicMaterial({
                    vertexColors: true,
                    side: THREE.DoubleSide,
                });
                o.castShadow = false;
                o.receiveShadow = false;
            }
        });

        // 新地形解析完成 → 先加入场景，再移除旧地形（避免空窗）
        SKY.scene.add(group);
        if (SKY.terrain) {
            SKY.scene.remove(SKY.terrain);
            SKY.terrain.traverse((o) => { if (o.geometry) o.geometry.dispose(); });
        }
        SKY.terrain = group;
        SKY.terrainKey = `${lat.toFixed(4)}|${lon.toFixed(4)}`;

        // 重置相机：眼高 2 m 于观测点正上方（Y-up 帧：高程沿 Y 轴），看向天空
        SKY.eyeH = 2;
        SKY.cameraPos.set(0, SKY.obsY + SKY.eyeH, 0);
        // 地平线环保持在观察者眼高水平面（天文地平线，0° 仰角参考）
        if (SKY.horizonRing) SKY.horizonRing.position.y = SKY.cameraPos.y;
        SKY.az = 0;
        SKY.el = 55;
        SKY.ready = true; // 先置 ready，updateSkyObjects 才能定位地球/鹊桥
        lookFromAzEl();
        updateSkyObjects();
        updateCommStatusBadge();
    } catch (err) {
        // 失败：保留旧地形，仅提示（不抛异常，避免白屏）
        console.error("terrain_gltf:", lat, lon, err);
    } finally {
        setSpinner(skySpinner, false);
    }
}

/* 依据当前时间 slider 的位置更新地球 / 鹊桥天空位置 */
function updateSkyObjects() {
    if (!SKY.ready || !SKY.scene) return;
    const earth = getCurrentEarth();
    const relay = getCurrentRelay();
    if (SKY.earthMesh) {
        if (earth) {
            const dir = skyDirFromAzEl(earth.az_deg || 0, earth.elev_deg || 0);
            const earthPos = SKY.cameraPos.clone().addScaledVector(dir, 15000);
            SKY.earthMesh.position.copy(earthPos);
            SKY.earthMesh.visible = true;
            // 线框幽灵球：与地球同位置，深度忽略 → 被地形遮挡时仍可见
            if (SKY.earthGhost) {
                SKY.earthGhost.position.copy(earthPos);
                SKY.earthGhost.visible = true;
            }
            // 地球视线光束：观察点 → 地球
            if (SKY.earthBeam) {
                const bGeo = SKY.earthBeam.geometry;
                bGeo.setFromPoints([SKY.cameraPos, earthPos]);
                bGeo.computeBoundingSphere();
                SKY.earthBeam.computeLineDistances();
                SKY.earthBeam.visible = true;
            }
        } else {
            SKY.earthMesh.visible = false;
            if (SKY.earthGhost) SKY.earthGhost.visible = false;
            if (SKY.earthBeam) SKY.earthBeam.visible = false;
        }
    }
    if (SKY.relayMesh) {
        if (relay) {
            const dir = skyDirFromAzEl(relay.az_deg || 0, relay.elev_deg || 0);
            const relayPos = SKY.cameraPos.clone().addScaledVector(dir, 20000);
            SKY.relayMesh.position.copy(relayPos);
            // ConeGeometry 局部轴为 +Y，需从局部 +Y 旋转到世界方向 dir
            const CONE_AXIS = new THREE.Vector3(0, 1, 0);
            SKY.relayMesh.quaternion.setFromUnitVectors(CONE_AXIS, dir);
            SKY.relayMesh.visible = true;
            // 鹊桥视线光束：观察点 → 鹊桥
            if (SKY.relayBeam) {
                const bGeo = SKY.relayBeam.geometry;
                bGeo.setFromPoints([SKY.cameraPos, relayPos]);
                bGeo.computeBoundingSphere();
                SKY.relayBeam.computeLineDistances();
                SKY.relayBeam.visible = true;
            }
        } else {
            SKY.relayMesh.visible = false;
            if (SKY.relayBeam) SKY.relayBeam.visible = false;
        }
    }
}

/* ── 鼠标拖动（第一人称环视）+ 滚轮（眼高缩放） ── */
skyCanvas.addEventListener("pointerdown", (e) => {
    SKY.drag = { x: e.clientX, y: e.clientY };
    skyCanvas.classList.add("dragging");
    skyCanvas.setPointerCapture(e.pointerId);
});
skyCanvas.addEventListener("pointermove", (e) => {
    if (!SKY.drag || !SKY.ready) return;
    const dx = e.clientX - SKY.drag.x;
    const dy = e.clientY - SKY.drag.y;
    SKY.drag = { x: e.clientX, y: e.clientY };
    // 第一人称：右拖 → 视角右转（az 增加），上拖 → 视角上仰（el 增加）
    SKY.az = ((SKY.az + dx * 0.25) % 360 + 360) % 360;
    SKY.el = Math.max(-85, Math.min(88, SKY.el - dy * 0.25));
    lookFromAzEl();
});
skyCanvas.addEventListener("pointerup", () => {
    SKY.drag = null;
    skyCanvas.classList.remove("dragging");
});
skyCanvas.addEventListener("pointercancel", () => {
    SKY.drag = null;
    skyCanvas.classList.remove("dragging");
});
skyCanvas.addEventListener(
    "wheel",
    (e) => {
        e.preventDefault();
        SKY.eyeH = Math.max(1, Math.min(200, SKY.eyeH * (1 + Math.sign(e.deltaY) * 0.1)));
        SKY.cameraPos.y = SKY.obsY + SKY.eyeH; // Y-up 帧：眼高沿 Y 轴
        SKY.camera.position.y = SKY.cameraPos.y;
        // 地平线环始终跟随观察者眼高水平面
        if (SKY.horizonRing) SKY.horizonRing.position.y = SKY.cameraPos.y;
        lookFromAzEl();
    },
    { passive: false }
);

window.addEventListener("resize", resizeSky);

/* ═══════════════════════ 3D 天空全屏查看 ═══════════════════════
 * 原生 Fullscreen API：全屏元素 = #skyOverlay（画布 + HUD + 按钮整体），
 * 元素保留在原 DOM 位置，由浏览器全屏层接管；拖拽/滚轮监听仍挂在
 * #skyCanvas 上，全屏后交互无需重新绑定。兼容 Safari 前缀 API。
 */
const skyOverlay = document.getElementById("skyOverlay");
const skyFsBtn = document.getElementById("skyFsBtn");

/* 当前是否处于 3D 全屏（兼容 webkit 前缀） */
function isSkyFullscreen() {
    return (
        document.fullscreenElement === skyOverlay ||
        document.webkitFullscreenElement === skyOverlay
    );
}

function enterSkyFullscreen() {
    if (!skyOverlay) return;
    const fn = skyOverlay.requestFullscreen || skyOverlay.webkitRequestFullscreen;
    if (fn) fn.call(skyOverlay);
}

function exitSkyFullscreen() {
    const fn = document.exitFullscreen || document.webkitExitFullscreen;
    if (fn) fn.call(document);
}

function toggleSkyFullscreen() {
    if (isSkyFullscreen()) exitSkyFullscreen();
    else enterSkyFullscreen();
}

/* 全屏进出：更新按钮文案 + 重算画布尺寸（渲染缓冲 & 相机宽高比）。
 * 用 rAF 延后到浏览器应用完全屏布局后再测量，避免读到旧尺寸。 */
function onSkyFullscreenChange() {
    const fs = isSkyFullscreen();
    if (skyFsBtn) {
        skyFsBtn.textContent = fs ? "⛶ 退出全屏" : "⛶ 全屏查看";
        skyFsBtn.title = fs ? "退出全屏（Esc）" : "全屏查看 3D 天空（Esc 退出）";
    }
    requestAnimationFrame(resizeSky);
}

if (skyFsBtn) {
    skyFsBtn.addEventListener("click", () => {
        trackClick("3D全屏");
        toggleSkyFullscreen();
    });
}
document.addEventListener("fullscreenchange", onSkyFullscreenChange);
document.addEventListener("webkitfullscreenchange", onSkyFullscreenChange);

/* ═══════════════════════ 初始化 ═══════════════════════ */
async function init() {
    // 0. 上报本次页面访问（IP + 日期，CSV 持久化）
    trackVisit();

    timeLabel.textContent = fmtTimeLabel(parseInt(timeSlider.value, 10));

    // 1. 默认观测点（Shackleton 坑底）预置为历史首项
    state.history.push({
        key: `${state.lat.toFixed(4)}|${state.lon.toFixed(4)}`,
        lat: state.lat,
        lon: state.lon,
        time: datePicker.value + " " + timeLabel.textContent,
    });
    renderHistory();

    // 2. 加载 DEM 瓦片（中心 = Shackleton 坑底）
    await loadMapTile();

    // 3. 加载地平线掩码（默认 Shackleton 坑底）
    await loadHorizonMask(state.lat, state.lon);

    // 4. 加载日序列（默认日期）
    await loadDailyData(state.lat, state.lon, datePicker.value);

    // 初始完整重绘 + 通信状态
    redrawPolar(getCurrentEarth(), getCurrentRelay());
    updateCommStatusBadge();
}

init();
