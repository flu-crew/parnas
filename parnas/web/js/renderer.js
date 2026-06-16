/**
 * paper.js canvas tree renderer.
 * Handles rectangular and radial layouts, pan/zoom, hit-testing,
 * annotation overlays, metadata strips, and clade collapse.
 *
 * Usage:
 *   const r = new TreeRenderer(canvasEl);
 *   r.render(root, coords, annotations, opts);
 *   r.setAnnotations(annotations);   // re-style without re-layout
 *   r.resetView();
 *   r.zoom(factor);
 *   r.destroy();
 */

import { leaves, postorder, allNodes } from "./treemodel.js";
import { polarToCartesian } from "./layout.js";

// ── Constants ──────────────────────────────────────────────────────────────
const LABEL_PAD   = 6;   // px between node dot and label start
const NODE_R      = 2.5; // default leaf dot radius
const REP_R       = 6;   // representative node radius
const STRIP_W     = 12;  // metadata strip column width (px)
const STRIP_GAP   = 3;   // gap between strips
const PAD_LEFT    = 40;
const PAD_TOP     = 20;
const PAD_BOTTOM  = 20;
const MAX_LABEL_W = 240; // cap on space reserved for labels

export class TreeRenderer {
  /**
   * @param {HTMLCanvasElement} canvas
   */
  constructor(canvas) {
    this._canvas = canvas;
    const p = window.paper;
    p.setup(canvas);
    this._p = p;

    // Layer stack (back to front)
    this._layers = {
      bg:      new p.Layer(),
      branch:  new p.Layer(),
      node:    new p.Layer(),
      label:   new p.Layer(),
      strip:   new p.Layer(),
      overlay: new p.Layer(), // clade triangles, user shapes, legend
    };

    this._root   = null;
    this._coords = null;  // Map<id, {x,y}> normalised [0,1]
    this._ann    = null;
    this._opts   = {};

    // id → paper.Item for click hit-testing
    this._hitNodes   = new Map();
    this._hitBranches = new Map();

    // Callbacks set by app
    this.onNodeClick   = null; // (nodeId, node, event) => void
    this.onBranchClick = null; // (nodeId, node, event) => void

    this._setupInteraction();
  }

  // ── Public API ─────────────────────────────────────────────────────────

  /** Full re-draw with new layout. */
  render(root, coords, annotations, opts = {}) {
    this._root   = root;
    this._coords = coords;
    this._ann    = annotations;
    this._opts   = opts;
    this._redraw();
  }

  /** Re-apply annotation styles; skips layout re-computation. */
  setAnnotations(annotations) {
    this._ann = annotations;
    this._redraw();
  }

  /** Notify renderer that canvas dimensions changed. Call before render(). */
  resize() {
    const c = this._canvas;
    this._p.view.viewSize = new this._p.Size(c.width, c.height);
    this._p.view.update();
  }

  resetView() {
    this._p.view.matrix = new this._p.Matrix();
    this._p.view.update();
  }

  zoom(factor) {
    const c = this._p.view.center;
    this._p.view.scale(factor, c);
    this._p.view.update();
  }

  fitLabels() {
    // Scale so the rightmost label edge touches canvas right edge
    this.resetView();
  }

  destroy() {
    this._p.project.clear();
    this._p.view.remove();
  }

  /** Export current view as PNG blob (offscreen scaled). */
  exportPng(widthPx, heightPx, bgColor) {
    return new Promise(resolve => {
      const off = document.createElement("canvas");
      off.width  = widthPx;
      off.height = heightPx;
      const ctx = off.getContext("2d");
      ctx.fillStyle = bgColor;
      ctx.fillRect(0, 0, widthPx, heightPx);
      ctx.drawImage(
        this._canvas, 0, 0, this._canvas.width, this._canvas.height,
        0, 0, widthPx, heightPx
      );
      off.toBlob(resolve, "image/png");
    });
  }

  /** Export as SVG string via paper.js. */
  exportSvg() {
    return this._p.project.exportSVG({ asString: true });
  }

  /** Re-render with a specific dark/light theme, export SVG, then restore. */
  exportSvgWithTheme(dark) {
    const origDark = this._opts.dark;
    if (origDark === dark) return this.exportSvg();
    this._opts = { ...this._opts, dark };
    this._redraw();
    const svg = this._p.project.exportSVG({ asString: true });
    this._opts = { ...this._opts, dark: origDark };
    this._redraw();
    return svg;
  }

  /**
   * Re-render with export theme, copy canvas pixels to an offscreen canvas
   * at target dimensions, then restore original theme.
   * Returns a Promise<Blob>.
   */
  exportPngWithTheme(dark, bgColor, targetW, targetH, fmt = "png") {
    const origDark = this._opts.dark;
    if (origDark !== dark) {
      this._opts = { ...this._opts, dark };
      this._redraw();
      this._p.view.update();
    }

    const src = this._canvas;
    const off = document.createElement("canvas");
    off.width  = targetW;
    off.height = targetH;
    const ctx = off.getContext("2d");
    ctx.fillStyle = bgColor;
    ctx.fillRect(0, 0, targetW, targetH);
    ctx.drawImage(src, 0, 0, src.width, src.height, 0, 0, targetW, targetH);

    const mime = fmt === "jpg" ? "image/jpeg" : "image/png";
    return new Promise(resolve => {
      off.toBlob(blob => {
        if (origDark !== dark) {
          this._opts = { ...this._opts, dark: origDark };
          this._redraw();
        }
        resolve(blob);
      }, mime, 0.95);
    });
  }

  // ── Internal ────────────────────────────────────────────────────────────

  _clearLayers() {
    for (const layer of Object.values(this._layers)) {
      layer.activate();
      layer.removeChildren();
    }
    this._hitNodes.clear();
    this._hitBranches.clear();
  }

  _pad() {
    const w = this._canvas.width  || this._canvas.offsetWidth  || 800;
    const h = this._canvas.height || this._canvas.offsetHeight || 600;
    const isRadial = this._opts.radial;

    if (isRadial) {
      return { w, h, left: 0, right: 0, top: 0, bottom: 0 };
    }

    const leafList = this._root ? leaves(this._root) : [];
    const maxLen   = Math.max(...leafList.map(l => (l.name || "").length), 0);
    const labelW   = Math.min(maxLen * 6, MAX_LABEL_W);
    const nStrips  = this._ann?.metadataStrips?.length || 0;
    const stripW   = nStrips * (STRIP_W + STRIP_GAP);

    return {
      w, h,
      left:   PAD_LEFT,
      right:  labelW + LABEL_PAD + stripW + 20,
      top:    PAD_TOP,
      bottom: PAD_BOTTOM,
    };
  }

  /** Normalised [0,1] → canvas px. */
  _toPx(nx, ny, pad) {
    if (this._opts.radial) {
      const cx = pad.w / 2;
      const cy = pad.h / 2;
      const r  = Math.min(pad.w, pad.h) / 2 * 0.82 * ny;
      return polarToCartesian(nx, r, cx, cy);
    }
    return {
      x: pad.left + nx * (pad.w - pad.left - pad.right),
      y: pad.top  + ny * (pad.h - pad.top  - pad.bottom),
    };
  }

  _effectiveColors() {
    // Build Map<nodeId, hex|null> by propagating leaf colours upward.
    const ann    = this._ann    || {};
    const opts   = this._opts   || {};
    const nc     = ann.nodeColors  || new Map();
    const repSet = opts.repsSet    || new Set();
    const repCol = opts.repColors  || {};

    const colorOf = new Map();

    // Leaves
    for (const leaf of leaves(this._root)) {
      const c = nc.get(leaf.name) || repCol[leaf.name] || null;
      colorOf.set(leaf.id, c);
    }

    // Internal nodes (postorder)
    postorder(this._root, node => {
      if (node.children.length === 0) return;
      const childColors = node.children.map(c => colorOf.get(c.id));
      const nonNull = childColors.filter(Boolean);
      const unique  = new Set(nonNull);
      colorOf.set(node.id,
        unique.size === 1 && nonNull.length === childColors.length
          ? [...unique][0] : null
      );
    });

    return colorOf;
  }

  _redraw() {
    if (!this._root || !this._coords) return;

    this._clearLayers();

    const p       = this._p;
    const pad     = this._pad();
    const ann     = this._ann    || {};
    const opts    = this._opts   || {};
    const isRadial = opts.radial || false;
    const isDark  = opts.dark !== false;
    const fontSize = opts.fontSize || 11;
    const coords  = this._coords;
    const root    = this._root;

    const branchDef = isDark ? "rgba(255,255,255,0.18)" : "rgba(0,0,0,0.22)";
    const textDef   = isDark ? "rgba(255,255,255,0.50)" : "rgba(0,0,0,0.55)";
    const repsSet   = opts.repsSet || new Set();

    const effectiveColor = this._effectiveColors();

    // Collapsed clade node ids
    const collapsedSet = new Set(
      (ann.cladeGroups || []).filter(g => g.collapsed).map(g => g.nodeId)
    );

    // ── Recursive draw ────────────────────────────────────────
    const drawSubtree = (node, parentPx, parentNc, isSkipped) => {
      const nc = coords.get(node.id);
      if (!nc) return;
      const px = this._toPx(nc.x, nc.y, pad);

      // Branch to parent
      if (parentPx) {
        this._layers.branch.activate();

        let path;
        if (!isRadial) {
          // Elbow: vertical stem then horizontal arm
          path = new p.Path([
            new p.Point(parentPx.x, parentPx.y),
            new p.Point(parentPx.x, px.y),
            new p.Point(px.x, px.y),
          ]);
        } else {
          // Radial: arc at parent radius then spoke outward
          const pR  = Math.min(pad.w, pad.h) / 2 * 0.82 * parentNc.y;
          const arc = polarToCartesian(nc.x, pR, pad.w / 2, pad.h / 2);
          path = new p.Path([
            new p.Point(parentPx.x, parentPx.y),
            new p.Point(arc.x, arc.y),
            new p.Point(px.x, px.y),
          ]);
        }

        const branchColor = ann.branchColors?.get(node.id)
          || effectiveColor.get(node.id)
          || branchDef;
        path.strokeColor = branchColor;
        path.strokeWidth = ann.branchWidths?.get(node.id) || 1.5;
        path.fillColor   = null;
        this._hitBranches.set(node.id, path);
      }

      const collapsed = collapsedSet.has(node.id) && node !== root;

      if (collapsed) {
        // Draw triangle for collapsed clade
        this._layers.overlay.activate();
        const leafList = leaves(node);
        if (leafList.length > 1) {
          const firstNc = coords.get(leafList[0].id);
          const lastNc  = coords.get(leafList[leafList.length - 1].id);
          if (firstNc && lastNc) {
            const tip = px;
            const top = this._toPx(1.0, firstNc.y, pad);
            const bot = this._toPx(1.0, lastNc.y,  pad);
            const tri = new p.Path([
              new p.Point(tip.x, tip.y),
              new p.Point(top.x, top.y),
              new p.Point(bot.x, bot.y),
            ]);
            tri.closed      = true;
            const triColor  = new p.Color(
              ann.branchColors?.get(node.id) || effectiveColor.get(node.id) || branchDef
            );
            tri.fillColor   = triColor;
            tri.fillColor.alpha = 0.30;
            tri.strokeColor = triColor;
            tri.strokeWidth = 1;

            // Clade label
            const cg = (ann.cladeGroups || []).find(g => g.nodeId === node.id);
            const label = cg?.label || `[${leafList.length}]`;
            this._layers.label.activate();
            const tx = new p.PointText(new p.Point(top.x + LABEL_PAD, (top.y + bot.y) / 2 + 4));
            tx.content   = label;
            tx.fontSize  = fontSize;
            tx.fillColor = triColor;
          }
        }
        return; // Don't recurse into collapsed clade
      }

      // Recurse children
      for (const child of node.children) {
        drawSubtree(child, px, nc, false);
      }

      // Node dot
      const isLeaf  = node.children.length === 0;
      const color   = effectiveColor.get(node.id);
      const isRep   = node.name && repsSet.has(node.name);
      const dotR    = isRep ? REP_R : (isLeaf && color ? NODE_R : (isLeaf ? 1.5 : 0));

      if (dotR > 0) {
        this._layers.node.activate();
        const dot = new p.Path.Circle(new p.Point(px.x, px.y), dotR);
        dot.fillColor   = color || textDef;
        dot.strokeColor = isRep ? (isDark ? "#ffffff" : "#000000") : null;
        dot.strokeWidth = isRep ? 1.5 : 0;
        this._hitNodes.set(node.id, dot);
      }

      // Leaf label
      if (isLeaf && node.name && !isRadial) {
        this._layers.label.activate();
        const displayName = ann.labelOverrides?.get(node.id) || node.name;
        const tx = new p.PointText(new p.Point(px.x + LABEL_PAD + dotR, px.y + 4));
        tx.content    = displayName;
        tx.fontSize   = fontSize;
        tx.fontFamily = isRep
          ? '"JetBrains Mono", monospace'
          : '"Outfit", sans-serif';
        tx.fontWeight  = isRep ? "bold" : "normal";
        tx.fillColor   = color || textDef;
      }

      // Radial leaf label
      if (isLeaf && node.name && isRadial) {
        this._layers.label.activate();
        const displayName = ann.labelOverrides?.get(node.id) || node.name;
        const angleDeg = nc.x * 180 / Math.PI;
        const flip     = angleDeg > 90 && angleDeg < 270;
        const tx = new p.PointText(new p.Point(px.x, px.y));
        tx.content   = displayName;
        tx.fontSize  = fontSize;
        tx.fillColor = color || textDef;
        tx.rotate(flip ? angleDeg + 180 : angleDeg);
        tx.translate(
          flip
            ? new p.Point(-tx.bounds.width - LABEL_PAD, 4)
            : new p.Point(LABEL_PAD, 4)
        );
      }

      // Internal node support value label (optional)
      if (!isLeaf && node.support !== null && node.support !== undefined && opts.showSupport) {
        this._layers.label.activate();
        const tx = new p.PointText(new p.Point(px.x + 2, px.y - 3));
        tx.content   = node.support.toFixed(0);
        tx.fontSize  = Math.max(fontSize - 2, 7);
        tx.fillColor = textDef;
      }

      // Clade group bracket/label (non-collapsed)
      const cg = (ann.cladeGroups || []).find(g => g.nodeId === node.id && !g.collapsed);
      if (cg?.label) {
        this._layers.overlay.activate();
        const tx = new p.PointText(new p.Point(px.x, px.y - 9));
        tx.content    = cg.label;
        tx.fontSize   = fontSize + 1;
        tx.fontWeight = "bold";
        tx.fillColor  = cg.color || effectiveColor.get(node.id) || textDef;
      }
    };

    drawSubtree(root, null, null, false);

    // ── Metadata strips ────────────────────────────────────────
    if (!isRadial && ann.metadataStrips?.length) {
      this._layers.strip.activate();
      const leafList  = leaves(root);
      const nLeaves   = leafList.length;
      const rowH      = Math.max(1.5, (pad.h - pad.top - pad.bottom) / nLeaves);
      // Start strips just after the label area
      const stripStartX = pad.w - pad.right + LABEL_PAD + 10;

      ann.metadataStrips.forEach((strip, si) => {
        const sx = stripStartX + si * (STRIP_W + STRIP_GAP);

        leafList.forEach(leaf => {
          const entry = strip.data?.get(leaf.name);
          if (!entry) return;
          const nc = coords.get(leaf.id);
          if (!nc) return;
          const lpx = this._toPx(nc.x, nc.y, pad);
          const rect = new p.Path.Rectangle(
            new p.Rectangle(sx, lpx.y - rowH / 2, STRIP_W, rowH)
          );
          rect.fillColor   = entry.color;
          rect.strokeColor = null;
        });

        // Column header (rotated)
        const hdr = new p.PointText(new p.Point(sx + STRIP_W / 2, pad.top - 6));
        hdr.content      = strip.name.slice(0, 10);
        hdr.fontSize     = 8;
        hdr.fillColor    = textDef;
        hdr.justification = "center";
      });
    }

    // ── Legend ─────────────────────────────────────────────────
    if (ann.legend?.length) {
      this._layers.overlay.activate();
      // Position at top-right, inside the label margin
      const legendW = 160;
      const lx = pad.w - Math.max(pad.right, legendW + 8);
      let   ly = pad.top + 10;
      ann.legend.forEach(entry => {
        const dot = new p.Path.Circle(new p.Point(lx + 5, ly + 5), 5);
        dot.fillColor = entry.color;
        const tx = new p.PointText(new p.Point(lx + 14, ly + 9));
        tx.content   = entry.label;
        tx.fontSize  = 10;
        tx.fillColor = entry.color;
        ly += 16;
      });
    }

    // ── User shapes ────────────────────────────────────────────
    this._layers.overlay.activate();
    for (const shape of (ann.shapes || [])) {
      if (shape.type === "rect") {
        const r = new p.Path.Rectangle(
          new p.Rectangle(shape.x, shape.y, shape.w || 60, shape.h || 20)
        );
        r.strokeColor = shape.color || "#00e5c8";
        r.strokeWidth = 1.5;
        r.fillColor   = null;
      } else if (shape.type === "text") {
        const tx = new p.PointText(new p.Point(shape.x, shape.y));
        tx.content   = shape.text || "";
        tx.fontSize  = 12;
        tx.fillColor = shape.color || "#00e5c8";
      } else if (shape.type === "line") {
        const ln = new p.Path.Line(
          new p.Point(shape.x, shape.y),
          new p.Point(shape.x2 || shape.x + 40, shape.y2 || shape.y)
        );
        ln.strokeColor = shape.color || "#00e5c8";
        ln.strokeWidth = 1.5;
      }
    }

    p.view.update();
  }

  // ── Interaction ────────────────────────────────────────────────────────

  _setupInteraction() {
    const p      = this._p;
    const canvas = this._canvas;
    let panning  = false;
    let lastPt   = null;

    canvas.addEventListener("mousedown", e => {
      if (e.button !== 0) return;
      panning = true;
      lastPt  = { x: e.clientX, y: e.clientY };
    });

    window.addEventListener("mousemove", e => {
      if (!panning || !lastPt) return;
      const dx = e.clientX - lastPt.x;
      const dy = e.clientY - lastPt.y;
      lastPt = { x: e.clientX, y: e.clientY };
      p.view.translate(new p.Point(dx, dy));
      p.view.update();
    });

    window.addEventListener("mouseup", e => {
      // Short drag = click; long drag = pan
      if (panning && lastPt) {
        const dx = Math.abs(e.clientX - lastPt.x);
        const dy = Math.abs(e.clientY - lastPt.y);
        if (dx + dy < 4) this._handleClick(e);
      }
      panning = false;
      lastPt  = null;
    });

    canvas.addEventListener("wheel", e => {
      if (!e.ctrlKey) return;
      e.preventDefault();
      const factor = e.deltaY < 0 ? 1.12 : 1 / 1.12;
      const mouse  = p.view.viewToProject(new p.Point(e.offsetX, e.offsetY));
      p.view.scale(factor, mouse);
      p.view.update();
    }, { passive: false });
  }

  _handleClick(e) {
    const p     = this._p;
    const rect  = this._canvas.getBoundingClientRect();
    const px    = e.clientX - rect.left;
    const py    = e.clientY - rect.top;
    const pt    = p.view.viewToProject(new p.Point(px, py));

    // Hit-test node dots first (smaller, higher priority)
    for (const [id, item] of this._hitNodes) {
      if (item.contains && item.contains(pt)) {
        const node = this._findNodeById(id);
        if (this.onNodeClick) this.onNodeClick(id, node, e);
        return;
      }
    }

    // Hit-test branches (wider tolerance via bounds)
    for (const [id, path] of this._hitBranches) {
      if (path.hitTest && path.hitTest(pt, { stroke: true, tolerance: 5 })) {
        const node = this._findNodeById(id);
        if (this.onBranchClick) this.onBranchClick(id, node, e);
        return;
      }
    }
  }

  _findNodeById(id) {
    if (!this._root) return null;
    let found = null;
    const search = node => {
      if (found) return;
      if (node.id === id) { found = node; return; }
      for (const c of node.children) search(c);
    };
    search(this._root);
    return found;
  }
}
