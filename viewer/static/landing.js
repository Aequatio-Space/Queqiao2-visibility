/* 首页访问统计上报（与 viewer/static/app.js 的 trackVisit/trackClick 同口径）。
 *
 * 页面加载时上报一次访问（POST /api/visit），三张功能区卡片点击时上报
 * （POST /api/click）；sendBeacon 保证跳转前送达，静默失败不影响导航。
 */
(function () {
    "use strict";

    function track(url) {
        if (navigator.sendBeacon) {
            navigator.sendBeacon(url);
        } else {
            fetch(url, { method: "POST", keepalive: true }).catch(function () {});
        }
    }

    track("/api/visit?path=" + encodeURIComponent(location.pathname));

    document.querySelectorAll("[data-track]").forEach(function (el) {
        el.addEventListener("click", function () {
            var name = el.getAttribute("data-track") || "unknown";
            track("/api/click?button=" + encodeURIComponent(name));
        });
    });
})();
