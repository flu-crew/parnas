/**
 * SweepChart — diversity-vs-n sparkline for the parameter sweep bar.
 */
export class SweepChart {
  constructor(canvas) {
    this._canvas = canvas;
    this._ctx    = canvas.getContext("2d");
    this._points = []; // [{n, diversity}]
    this._currentN = null;
    this._rect   = null; // cached bounding rect; invalidated on resize
  }

  addPoint(n, diversity) {
    if (diversity == null) return;
    // Replace existing point for same n
    const idx = this._points.findIndex(p => p.n === n);
    if (idx >= 0) this._points[idx].diversity = diversity;
    else this._points.push({ n, diversity });
    this._points.sort((a, b) => a.n - b.n);
    this._currentN = n;
    this._canvas.style.display = this._points.length >= 2 ? "block" : "none";
    this.draw();
  }

  /** Call when the canvas is resized so the cached rect is re-measured on next draw. */
  invalidateRect() { this._rect = null; }

  reset() {
    this._points = [];
    this._currentN = null;
    this._rect = null;
    const ctx = this._ctx;
    ctx.clearRect(0, 0, this._canvas.width, this._canvas.height);
    this._canvas.style.display = "none";
  }

  draw() {
    const pts = this._points;
    if (pts.length < 2) return;

    const canvas = this._canvas;
    const ctx    = this._ctx;
    // Size backing buffer to displayed CSS box at device pixel ratio (prevents clipping/stretch)
    // Cache the rect — re-measure only when explicitly invalidated (resize)
    if (!this._rect) this._rect = canvas.getBoundingClientRect();
    const rect = this._rect;
    const dpr  = window.devicePixelRatio || 1;
    canvas.width  = Math.round(rect.width  * dpr);
    canvas.height = Math.round(rect.height * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    const W = rect.width;
    const H = rect.height;
    const PAD = { l: 34, r: 8, t: 6, b: 20 };

    ctx.clearRect(0, 0, W, H);

    const minN   = pts[0].n;
    const maxN   = pts[pts.length - 1].n;
    const minDiv = Math.min(...pts.map(p => p.diversity));
    const maxDiv = Math.max(...pts.map(p => p.diversity));
    const rangeDiv = maxDiv - minDiv || 1;

    const toX = n   => PAD.l + (n - minN) / Math.max(maxN - minN, 1) * (W - PAD.l - PAD.r);
    const toY = div => PAD.t + (1 - (div - minDiv) / rangeDiv) * (H - PAD.t - PAD.b);

    // Detect theme
    const dark = document.documentElement.getAttribute("data-theme") !== "light";
    const lineColor  = dark ? "rgba(0,229,200,0.85)" : "rgba(0,140,120,0.9)";
    const gridColor  = dark ? "rgba(255,255,255,0.08)" : "rgba(0,0,0,0.06)";
    const textColor  = dark ? "rgba(255,255,255,0.45)" : "rgba(0,0,0,0.40)";
    const markerColor = dark ? "#ffffff" : "#000000";

    // Grid lines at min/max diversity
    ctx.strokeStyle = gridColor;
    ctx.lineWidth   = 1;
    ctx.beginPath();
    ctx.moveTo(PAD.l, toY(minDiv)); ctx.lineTo(W - PAD.r, toY(minDiv));
    ctx.moveTo(PAD.l, toY(maxDiv)); ctx.lineTo(W - PAD.r, toY(maxDiv));
    ctx.stroke();

    // Axis labels
    ctx.fillStyle = textColor;
    ctx.font      = "9px sans-serif";
    ctx.textAlign = "right";
    ctx.fillText(maxDiv.toFixed(1) + "%", PAD.l - 2, toY(maxDiv) + 3);
    ctx.fillText(minDiv.toFixed(1) + "%", PAD.l - 2, toY(minDiv) + 3);
    ctx.textAlign = "center";
    ctx.fillText(minN, toX(minN), H - 4);
    ctx.fillText(maxN, toX(maxN), H - 4);

    // Line
    ctx.strokeStyle = lineColor;
    ctx.lineWidth   = 1.5;
    ctx.lineJoin    = "round";
    ctx.beginPath();
    pts.forEach((p, i) => {
      if (i === 0) ctx.moveTo(toX(p.n), toY(p.diversity));
      else         ctx.lineTo(toX(p.n), toY(p.diversity));
    });
    ctx.stroke();

    // Dots
    for (const p of pts) {
      ctx.beginPath();
      ctx.arc(toX(p.n), toY(p.diversity), 2.5, 0, Math.PI * 2);
      ctx.fillStyle = lineColor;
      ctx.fill();
    }

    // Current-n marker
    if (this._currentN != null) {
      const cx = toX(this._currentN);
      ctx.strokeStyle = markerColor;
      ctx.lineWidth   = 1;
      ctx.setLineDash([2, 2]);
      ctx.beginPath();
      ctx.moveTo(cx, PAD.t);
      ctx.lineTo(cx, H - PAD.b);
      ctx.stroke();
      ctx.setLineDash([]);
    }
  }
}
