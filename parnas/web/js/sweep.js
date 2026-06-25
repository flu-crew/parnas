/**
 * SweepChart — generalized X/Y curve for parameter sweeps (n-diversity or threshold-reps).
 * Renders on-screen and exports publication-ready PNG, SVG, and EPS.
 */
import { downloadBlob, downloadText, canvasToEps } from "./export.js";

export class SweepChart {
  constructor(canvas) {
    this._canvas  = canvas;
    this._ctx     = canvas.getContext("2d");
    this._points  = []; // [{x, y}]
    this._currentX  = null;
    this._elbowX    = null;
    this._elbowAxis = "x";
    this._title     = "";
    this._xTitle   = "x";
    this._yTitle   = "y";
    this._xFmt     = v => String(v);
    this._yFmt     = v => v.toFixed(1);
    this._rect     = null;
  }

  /**
   * Load a complete series at once (used after a full sweep completes).
   * @param {Array<{x:number,y:number}>} points
   * @param {{title?:string, xTitle?:string, yTitle?:string,
   *           xFmt?:(v:number)=>string, yFmt?:(v:number)=>string}} opts
   */
  setSeries(points, { title = "", xTitle = "x", yTitle = "y",
                      xFmt = null, yFmt = null } = {}) {
    this._points  = [...points].sort((a, b) => a.x - b.x);
    this._title   = title;
    this._xTitle  = xTitle;
    this._yTitle  = yTitle;
    if (xFmt) this._xFmt = xFmt;
    if (yFmt) this._yFmt = yFmt;
    this._currentX  = null;
    this._elbowX    = null;
    this._elbowAxis = "x";
    this._rect      = null;
    this._canvas.style.display = this._points.length >= 2 ? "block" : "none";
    if (this._points.length >= 2) this.draw();
  }

  /**
   * Append a single point (used by the live +/− stepper in n-mode).
   */
  addPoint(x, y) {
    if (y == null) return;
    const idx = this._points.findIndex(p => p.x === x);
    if (idx >= 0) this._points[idx].y = y;
    else this._points.push({ x, y });
    this._points.sort((a, b) => a.x - b.x);
    this._currentX = x;
    this._canvas.style.display = this._points.length >= 2 ? "block" : "none";
    if (this._points.length >= 2) this.draw();
  }

  setCurrentX(x) {
    this._currentX = x;
    if (this._points.length >= 2) this.draw();
  }

  invalidateRect() { this._rect = null; }

  setElbow(x, { axis = "x" } = {}) {
    this._elbowX    = x;
    this._elbowAxis = axis;
    if (this._points.length >= 2) this.draw();
  }

  reset() {
    this._points    = [];
    this._currentX  = null;
    this._elbowX    = null;
    this._elbowAxis = "x";
    this._rect      = null;
    this._ctx.clearRect(0, 0, this._canvas.width, this._canvas.height);
    this._canvas.style.display = "none";
  }

  // ── Layout ────────────────────────────────────────────────────────────────

  _layout(W, H) {
    const PAD = { l: 58, r: 18, t: this._title ? 36 : 14, b: 48 };
    const pts  = this._points;
    if (!pts.length) return null;

    const xs = pts.map(p => p.x);
    const ys = pts.map(p => p.y);
    const minX = Math.min(...xs), maxX = Math.max(...xs);
    const minY = Math.min(...ys), maxY = Math.max(...ys);
    const rangeX = maxX - minX || 1;
    const rangeY = maxY - minY || 1;

    const plotW = W - PAD.l - PAD.r;
    const plotH = H - PAD.t - PAD.b;

    const toX = x => PAD.l + (x - minX) / rangeX * plotW;
    const toY = y => PAD.t + (1 - (y - minY) / rangeY) * plotH;

    // x ticks — up to 8, prefer integers
    const nXticks = Math.min(8, xs.length);
    const xStep   = Math.ceil(rangeX / (nXticks - 1));
    const xTicks  = [];
    for (let v = minX; v <= maxX + 1e-9; v += xStep) xTicks.push(v);
    if (!xTicks.includes(maxX)) xTicks.push(maxX);

    // y ticks — 5 evenly spaced
    const yTicks = [];
    for (let i = 0; i <= 4; i++) yTicks.push(minY + rangeY * i / 4);

    return { PAD, minX, maxX, minY, maxY, rangeX, rangeY, plotW, plotH, toX, toY, xTicks, yTicks };
  }

  // ── Core render ───────────────────────────────────────────────────────────

  _renderFigure(ctx, W, H, { dark = true, background = null } = {}) {
    const lay = this._layout(W, H);
    if (!lay) return;
    const { PAD, toX, toY, xTicks, yTicks, minY, maxY } = lay;
    const pts = this._points;

    const bg        = background || (dark ? "#1a1a2e" : "#ffffff");
    const lineColor = dark ? "rgba(0,229,200,0.9)"  : "rgba(0,140,120,0.95)";
    const gridColor = dark ? "rgba(255,255,255,0.1)" : "rgba(0,0,0,0.08)";
    const textColor = dark ? "rgba(255,255,255,0.6)" : "rgba(0,0,0,0.55)";
    const titleClr  = dark ? "rgba(255,255,255,0.9)" : "rgba(0,0,0,0.85)";
    const markerClr = dark ? "rgba(255,255,255,0.5)" : "rgba(0,0,0,0.4)";
    const elbowClr  = dark ? "rgba(255,180,0,0.95)"  : "rgba(200,120,0,0.95)";
    const axisClr   = dark ? "rgba(255,255,255,0.2)" : "rgba(0,0,0,0.15)";

    // Background
    ctx.fillStyle = bg;
    ctx.fillRect(0, 0, W, H);

    // Title
    if (this._title) {
      ctx.fillStyle  = titleClr;
      ctx.font       = `bold ${Math.round(W * 0.018)}px sans-serif`;
      ctx.textAlign  = "center";
      ctx.textBaseline = "top";
      ctx.fillText(this._title, W / 2, 10);
    }

    // Axis lines
    ctx.strokeStyle = axisClr;
    ctx.lineWidth   = 1;
    ctx.beginPath();
    ctx.moveTo(PAD.l, PAD.t); ctx.lineTo(PAD.l, H - PAD.b);
    ctx.moveTo(PAD.l, H - PAD.b); ctx.lineTo(W - PAD.r, H - PAD.b);
    ctx.stroke();

    const tickFontSz = Math.max(9, Math.round(W * 0.011));

    // Y gridlines + labels
    ctx.fillStyle   = textColor;
    ctx.font        = `${tickFontSz}px sans-serif`;
    ctx.textBaseline = "middle";
    for (const yv of yTicks) {
      const py = toY(yv);
      ctx.strokeStyle = gridColor;
      ctx.lineWidth   = 1;
      ctx.beginPath(); ctx.moveTo(PAD.l, py); ctx.lineTo(W - PAD.r, py); ctx.stroke();
      ctx.textAlign = "right";
      ctx.fillText(this._yFmt(yv), PAD.l - 5, py);
    }

    // X tick labels
    ctx.textBaseline = "top";
    for (const xv of xTicks) {
      const px = toX(xv);
      ctx.textAlign = "center";
      ctx.fillStyle = textColor;
      ctx.fillText(this._xFmt(xv), px, H - PAD.b + 5);
    }

    // Axis titles
    const axTitleSz = Math.max(11, Math.round(W * 0.013));
    ctx.fillStyle = titleClr;
    ctx.font      = `${axTitleSz}px sans-serif`;
    ctx.textAlign = "center";
    ctx.textBaseline = "bottom";
    ctx.fillText(this._xTitle, PAD.l + lay.plotW / 2, H - 4);

    // Y axis title (rotated)
    ctx.save();
    ctx.translate(13, PAD.t + lay.plotH / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign    = "center";
    ctx.textBaseline = "top";
    ctx.fillText(this._yTitle, 0, 0);
    ctx.restore();

    // Curve
    ctx.strokeStyle = lineColor;
    ctx.lineWidth   = Math.max(1.5, W * 0.0015);
    ctx.lineJoin    = "round";
    ctx.beginPath();
    pts.forEach((p, i) => {
      const px = toX(p.x), py = toY(p.y);
      if (i === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
    });
    ctx.stroke();

    // Dots
    const dotR = Math.max(2.5, W * 0.003);
    for (const p of pts) {
      ctx.beginPath();
      ctx.arc(toX(p.x), toY(p.y), dotR, 0, Math.PI * 2);
      ctx.fillStyle = lineColor;
      ctx.fill();
    }

    // Elbow marker
    if (this._elbowX != null) {
      ctx.strokeStyle = elbowClr;
      ctx.lineWidth   = Math.max(1.5, W * 0.002);
      ctx.fillStyle   = elbowClr;
      ctx.font        = `bold ${tickFontSz + 1}px sans-serif`;
      if (this._elbowAxis === "y") {
        const ey = toY(this._elbowX);
        ctx.beginPath(); ctx.moveTo(PAD.l, ey); ctx.lineTo(W - PAD.r, ey); ctx.stroke();
        ctx.textAlign    = "left";
        ctx.textBaseline = "bottom";
        ctx.fillText(`elbow: ${this._yFmt(this._elbowX)}`, PAD.l + 4, ey - 3);
      } else {
        const ex = toX(this._elbowX);
        ctx.beginPath(); ctx.moveTo(ex, PAD.t); ctx.lineTo(ex, H - PAD.b); ctx.stroke();
        ctx.textAlign    = "center";
        ctx.textBaseline = "top";
        ctx.fillText(`elbow: ${this._xFmt(this._elbowX)}`, ex, PAD.t + 4);
      }
    }

    // Current-x marker (dashed)
    if (this._currentX != null && this._currentX !== this._elbowX) {
      const cx = toX(this._currentX);
      ctx.strokeStyle = markerClr;
      ctx.lineWidth   = 1;
      ctx.setLineDash([3, 3]);
      ctx.beginPath(); ctx.moveTo(cx, PAD.t); ctx.lineTo(cx, H - PAD.b); ctx.stroke();
      ctx.setLineDash([]);
    }
  }

  // ── On-screen draw ────────────────────────────────────────────────────────

  draw() {
    const canvas = this._canvas;
    if (!this._rect) this._rect = canvas.getBoundingClientRect();
    const rect = this._rect;
    const dpr  = window.devicePixelRatio || 1;
    canvas.width  = Math.round(rect.width  * dpr);
    canvas.height = Math.round(rect.height * dpr);
    const ctx = this._ctx;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    const W = rect.width, H = rect.height;
    const dark = document.documentElement.getAttribute("data-theme") !== "light";
    const bg   = getComputedStyle(document.documentElement)
                   .getPropertyValue("--pico-card-background-color").trim()
                 || (dark ? "#1e1e2e" : "#ffffff");
    this._renderFigure(ctx, W, H, { dark, background: bg });
  }

  // ── SVG export ────────────────────────────────────────────────────────────

  toSVG(W, H, { dark = false } = {}) {
    const lay = this._layout(W, H);
    if (!lay) return "";
    const { PAD, toX, toY, xTicks, yTicks } = lay;
    const pts = this._points;

    const bg        = dark ? "#1a1a2e" : "#ffffff";
    const lineColor = dark ? "rgba(0,229,200,0.9)"  : "rgba(0,140,120,0.95)";
    const gridColor = dark ? "rgba(255,255,255,0.1)" : "rgba(0,0,0,0.08)";
    const textColor = dark ? "rgba(255,255,255,0.6)" : "rgba(0,0,0,0.55)";
    const titleClr  = dark ? "rgba(255,255,255,0.9)" : "rgba(0,0,0,0.85)";
    const elbowClr  = dark ? "rgba(255,180,0,0.95)"  : "rgba(200,120,0,0.95)";
    const axisClr   = dark ? "rgba(255,255,255,0.2)" : "rgba(0,0,0,0.15)";

    const e = (s) => String(s).replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;");
    const tickFontSz = Math.max(9, Math.round(W * 0.011));
    const axTitleSz  = Math.max(11, Math.round(W * 0.013));
    const dotR       = Math.max(2.5, W * 0.003);

    let s = `<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H}" viewBox="0 0 ${W} ${H}">\n`;
    s += `<rect width="${W}" height="${H}" fill="${e(bg)}"/>\n`;

    // Title
    if (this._title) {
      const tsz = Math.round(W * 0.018);
      s += `<text x="${W/2}" y="${10+tsz}" text-anchor="middle" font-size="${tsz}" font-weight="bold" fill="${e(titleClr)}" font-family="sans-serif">${e(this._title)}</text>\n`;
    }

    // Axis lines
    s += `<line x1="${PAD.l}" y1="${PAD.t}" x2="${PAD.l}" y2="${H-PAD.b}" stroke="${e(axisClr)}" stroke-width="1"/>\n`;
    s += `<line x1="${PAD.l}" y1="${H-PAD.b}" x2="${W-PAD.r}" y2="${H-PAD.b}" stroke="${e(axisClr)}" stroke-width="1"/>\n`;

    // Y gridlines + labels
    for (const yv of yTicks) {
      const py = toY(yv);
      s += `<line x1="${PAD.l}" y1="${py.toFixed(1)}" x2="${W-PAD.r}" y2="${py.toFixed(1)}" stroke="${e(gridColor)}" stroke-width="1"/>\n`;
      s += `<text x="${PAD.l-5}" y="${py.toFixed(1)}" text-anchor="end" dominant-baseline="middle" font-size="${tickFontSz}" fill="${e(textColor)}" font-family="sans-serif">${e(this._yFmt(yv))}</text>\n`;
    }

    // X tick labels
    for (const xv of xTicks) {
      const px = toX(xv);
      s += `<text x="${px.toFixed(1)}" y="${H-PAD.b+5}" text-anchor="middle" dominant-baseline="hanging" font-size="${tickFontSz}" fill="${e(textColor)}" font-family="sans-serif">${e(this._xFmt(xv))}</text>\n`;
    }

    // Axis titles
    s += `<text x="${PAD.l+lay.plotW/2}" y="${H-4}" text-anchor="middle" dominant-baseline="auto" font-size="${axTitleSz}" fill="${e(titleClr)}" font-family="sans-serif">${e(this._xTitle)}</text>\n`;
    s += `<text x="0" y="0" text-anchor="middle" font-size="${axTitleSz}" fill="${e(titleClr)}" font-family="sans-serif" transform="translate(13,${PAD.t+lay.plotH/2}) rotate(-90)">${e(this._yTitle)}</text>\n`;

    // Curve
    const ptStr = pts.map(p => `${toX(p.x).toFixed(1)},${toY(p.y).toFixed(1)}`).join(" ");
    s += `<polyline points="${ptStr}" fill="none" stroke="${e(lineColor)}" stroke-width="${Math.max(1.5,W*0.0015).toFixed(1)}" stroke-linejoin="round"/>\n`;

    // Dots
    for (const p of pts) {
      s += `<circle cx="${toX(p.x).toFixed(1)}" cy="${toY(p.y).toFixed(1)}" r="${dotR.toFixed(1)}" fill="${e(lineColor)}"/>\n`;
    }

    // Elbow
    if (this._elbowX != null) {
      const sw = Math.max(1.5, W * 0.002).toFixed(1);
      if (this._elbowAxis === "y") {
        const ey = toY(this._elbowX).toFixed(1);
        s += `<line x1="${PAD.l}" y1="${ey}" x2="${W-PAD.r}" y2="${ey}" stroke="${e(elbowClr)}" stroke-width="${sw}"/>\n`;
        s += `<text x="${PAD.l+4}" y="${toY(this._elbowX)-3}" text-anchor="start" dominant-baseline="auto" font-size="${tickFontSz+1}" font-weight="bold" fill="${e(elbowClr)}" font-family="sans-serif">elbow: ${e(this._yFmt(this._elbowX))}</text>\n`;
      } else {
        const ex = toX(this._elbowX).toFixed(1);
        s += `<line x1="${ex}" y1="${PAD.t}" x2="${ex}" y2="${H-PAD.b}" stroke="${e(elbowClr)}" stroke-width="${sw}"/>\n`;
        s += `<text x="${ex}" y="${PAD.t+4}" text-anchor="middle" dominant-baseline="hanging" font-size="${tickFontSz+1}" font-weight="bold" fill="${e(elbowClr)}" font-family="sans-serif">elbow: ${e(this._xFmt(this._elbowX))}</text>\n`;
      }
    }

    // Current-x dashed
    if (this._currentX != null && this._currentX !== this._elbowX) {
      const cx = toX(this._currentX).toFixed(1);
      s += `<line x1="${cx}" y1="${PAD.t}" x2="${cx}" y2="${H-PAD.b}" stroke="rgba(128,128,128,0.5)" stroke-width="1" stroke-dasharray="3 3"/>\n`;
    }

    s += `</svg>`;
    return s;
  }

  // ── Export ────────────────────────────────────────────────────────────────

  exportFigure(fmt, { dark = false } = {}) {
    const W = 1200, H = 750;
    if (fmt === "svg") {
      downloadText(this.toSVG(W, H, { dark }), "parnas-sweep.svg", "image/svg+xml");
      return;
    }
    // PNG or EPS — render to offscreen canvas @2×
    const off = document.createElement("canvas");
    off.width  = W * 2;
    off.height = H * 2;
    const ctx = off.getContext("2d");
    ctx.setTransform(2, 0, 0, 2, 0, 0);
    this._renderFigure(ctx, W, H, { dark, background: dark ? "#1a1a2e" : "#ffffff" });

    if (fmt === "eps") {
      const eps = canvasToEps(off, W, H);
      downloadText(eps, "parnas-sweep.eps", "application/postscript");
    } else {
      off.toBlob(blob => downloadBlob(blob, "parnas-sweep.png"), "image/png");
    }
  }
}
