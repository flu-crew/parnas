/**
 * Annotation panel — Visual Options drawer (PearTree-inspired).
 * Wires HTML panel controls to annotation state mutations.
 *
 * Usage:
 *   const panel = new AnnotationPanel(containerEl, onAnnotationsChange);
 *   panel.setTree(root, nameIndex);
 *   panel.setAnnotations(annotations);
 *   // When user changes something, onAnnotationsChange(newAnnotations) fires.
 */

import {
  setBranchColor, setBranchWidth, setNodeColor,
  upsertCladeGroup, removeCladeGroup, setLabelOverride,
  upsertMetadataStrip, removeMetadataStrip, addShape, setLegend,
  parseMetadataCsv,
} from "./annotations.js";
import { mrca, subtreeLeafNames } from "./treemodel.js";

export class AnnotationPanel {
  /**
   * @param {HTMLElement} el   - container element for the panel
   * @param {function}    onChange  - (newAnnotations) => void
   */
  constructor(el, onChange) {
    this._el       = el;
    this._onChange = onChange;
    this._ann      = null;
    this._root     = null;
    this._nameIndex = null;

    // Currently selected node id (set by renderer click)
    this._selectedNodeId = null;
    this._selectedNode   = null;

    this._build();
  }

  /** Called when tree is loaded. */
  setTree(root, nameIndex) {
    this._root      = root;
    this._nameIndex = nameIndex;
  }

  /** Called when annotation state changes externally. */
  setAnnotations(ann) {
    this._ann = ann;
    this._syncUI();
  }

  /** Called by renderer when user clicks a node. */
  selectNode(nodeId, node) {
    this._selectedNodeId = nodeId;
    this._selectedNode   = node;
    this._syncNodeSection();
    this._openSection("sel-section");
  }

  /** Called by renderer when user clicks a branch. */
  selectBranch(nodeId, node) {
    this._selectedNodeId = nodeId;
    this._selectedNode   = node;
    this._syncBranchSection();
    this._openSection("branch-section");
  }

  // ── Build HTML ────────────────────────────────────────────────────────

  _build() {
    this._el.innerHTML = `
<div class="vopt-panel">

  <!-- Header -->
  <div class="vopt-header">
    <span>Visual Options</span>
    <button class="vopt-close" id="vopt-close-btn" aria-label="Close">✕</button>
  </div>

  <!-- ── Selection: node/clade ────────────────────────────── -->
  <details class="vopt-section" id="sel-section" open>
    <summary>Selection</summary>
    <div class="vopt-body">
      <p class="vopt-hint" id="sel-hint">Click a node or branch on the tree.</p>
      <div id="sel-controls" style="display:none">
        <label class="vopt-label">Node</label>
        <span class="vopt-mono" id="sel-name">—</span>

        <label class="vopt-label" style="margin-top:0.6rem">Node colour</label>
        <div class="vopt-row">
          <input type="color" id="sel-node-color" class="vopt-color-swatch"/>
          <button class="vopt-btn-sm" id="sel-node-color-clear">Clear</button>
        </div>

        <label class="vopt-label" style="margin-top:0.5rem">Label override</label>
        <div class="vopt-row">
          <input type="text" id="sel-label-input" class="vopt-input" placeholder="Custom label"/>
          <button class="vopt-btn-sm" id="sel-label-clear">Clear</button>
        </div>

        <div class="vopt-divider"></div>
        <label class="vopt-label">Clade group</label>
        <div class="vopt-row">
          <input type="text" id="clade-label-input" class="vopt-input" placeholder="Group name"/>
          <input type="color" id="clade-color" class="vopt-color-swatch" value="#00e5c8"/>
        </div>
        <div class="vopt-row" style="margin-top:0.3rem">
          <button class="vopt-btn-sm primary" id="clade-label-set">Set label</button>
          <button class="vopt-btn-sm" id="clade-collapse-btn">Collapse</button>
          <button class="vopt-btn-sm danger" id="clade-remove-btn">Remove</button>
        </div>
      </div>
    </div>
  </details>

  <!-- ── Branch styling ───────────────────────────────────── -->
  <details class="vopt-section" id="branch-section">
    <summary>Branch / Edge</summary>
    <div class="vopt-body">
      <p class="vopt-hint" id="branch-hint">Click a branch on the tree.</p>
      <div id="branch-controls" style="display:none">
        <label class="vopt-label">Colour</label>
        <div class="vopt-row">
          <input type="color" id="branch-color" class="vopt-color-swatch"/>
          <button class="vopt-btn-sm" id="branch-color-clear">Clear</button>
        </div>
        <label class="vopt-label" style="margin-top:0.5rem">Width (px)</label>
        <input type="range" id="branch-width" min="0.5" max="6" step="0.5" value="1.5"
               style="width:100%"/>
        <div style="text-align:right;font-size:0.7rem;opacity:0.6" id="branch-width-val">1.5 px</div>
        <div class="vopt-row" style="margin-top:0.4rem">
          <button class="vopt-btn-sm primary" id="branch-apply">Apply to branch</button>
          <button class="vopt-btn-sm" id="branch-apply-clade">Apply to clade</button>
        </div>
      </div>
    </div>
  </details>

  <!-- ── Display options ──────────────────────────────────── -->
  <details class="vopt-section" id="display-section">
    <summary>Display</summary>
    <div class="vopt-body">
      <label class="vopt-switch-row">
        <input type="checkbox" id="opt-show-support" role="switch"/>
        <span>Support values</span>
      </label>
      <label class="vopt-switch-row">
        <input type="checkbox" id="opt-radial" role="switch"/>
        <span>Radial layout</span>
      </label>
      <label class="vopt-label" style="margin-top:0.5rem">Font size</label>
      <input type="range" id="opt-font-size" min="6" max="16" step="1" value="11"
             style="width:100%"/>
      <div style="text-align:right;font-size:0.7rem;opacity:0.6" id="opt-font-size-val">11 px</div>
    </div>
  </details>

  <!-- ── Metadata strips ──────────────────────────────────── -->
  <details class="vopt-section" id="meta-section">
    <summary>Metadata Strips</summary>
    <div class="vopt-body">
      <label class="vopt-label">Import CSV/TSV</label>
      <p class="vopt-hint">Columns: taxon, value[, color]</p>
      <div class="vopt-row">
        <button class="vopt-btn-sm primary" id="meta-import-btn">Choose file…</button>
        <span class="vopt-mono" id="meta-file-name" style="opacity:0.5">No file</span>
      </div>
      <input type="file" id="meta-file-input" accept=".csv,.tsv,.txt" style="display:none"/>
      <div id="meta-strip-list" style="margin-top:0.5rem"></div>
    </div>
  </details>

  <!-- ── Legend ───────────────────────────────────────────── -->
  <details class="vopt-section" id="legend-section">
    <summary>Legend</summary>
    <div class="vopt-body">
      <div id="legend-list"></div>
      <button class="vopt-btn-sm" id="legend-clear-btn" style="margin-top:0.4rem">Clear legend</button>
    </div>
  </details>

  <!-- ── Add shapes ───────────────────────────────────────── -->
  <details class="vopt-section" id="shapes-section">
    <summary>Annotations / Shapes</summary>
    <div class="vopt-body">
      <div class="vopt-row">
        <button class="vopt-btn-sm" id="shape-add-text">+ Text</button>
        <button class="vopt-btn-sm" id="shape-add-rect">+ Box</button>
        <button class="vopt-btn-sm" id="shape-add-line">+ Line</button>
      </div>
      <div id="shapes-list" style="margin-top:0.5rem"></div>
    </div>
  </details>

</div>`;

    this._wire();
  }

  _wire() {
    const $ = id => document.getElementById(id);

    $("vopt-close-btn").addEventListener("click", () => {
      this._el.style.display = "none";
    });

    // ── Node colour ──────────────────────────────────────────
    $("sel-node-color").addEventListener("input", e => {
      if (!this._selectedNodeId) return;
      this._emit(setNodeColor(this._ann, this._selectedNode?.name || this._selectedNodeId, e.target.value));
    });
    $("sel-node-color-clear").addEventListener("click", () => {
      if (!this._selectedNode) return;
      this._emit(setNodeColor(this._ann, this._selectedNode?.name || this._selectedNodeId, null));
    });

    // ── Label override ───────────────────────────────────────
    $("sel-label-input").addEventListener("change", e => {
      if (!this._selectedNodeId) return;
      this._emit(setLabelOverride(this._ann, this._selectedNodeId, e.target.value || null));
    });
    $("sel-label-clear").addEventListener("click", () => {
      if (!this._selectedNodeId) return;
      $("sel-label-input").value = "";
      this._emit(setLabelOverride(this._ann, this._selectedNodeId, null));
    });

    // ── Clade group ──────────────────────────────────────────
    $("clade-label-set").addEventListener("click", () => {
      if (!this._selectedNodeId) return;
      const label = $("clade-label-input").value.trim();
      const color = $("clade-color").value;
      this._emit(upsertCladeGroup(this._ann, {
        nodeId: this._selectedNodeId, label, color, collapsed: false,
      }));
    });

    $("clade-collapse-btn").addEventListener("click", () => {
      if (!this._selectedNodeId) return;
      const existing = (this._ann.cladeGroups || []).find(g => g.nodeId === this._selectedNodeId);
      const label = $("clade-label-input").value.trim() || existing?.label || "";
      const color = $("clade-color").value;
      const wasCollapsed = existing?.collapsed || false;
      this._emit(upsertCladeGroup(this._ann, {
        nodeId: this._selectedNodeId, label, color, collapsed: !wasCollapsed,
      }));
      $("clade-collapse-btn").textContent = wasCollapsed ? "Collapse" : "Expand";
    });

    $("clade-remove-btn").addEventListener("click", () => {
      if (!this._selectedNodeId) return;
      this._emit(removeCladeGroup(this._ann, this._selectedNodeId));
    });

    // ── Branch colour / width ────────────────────────────────
    $("branch-width").addEventListener("input", e => {
      $("branch-width-val").textContent = e.target.value + " px";
    });

    $("branch-apply").addEventListener("click", () => {
      if (!this._selectedNodeId) return;
      let ann = this._ann;
      const color = $("branch-color").value;
      const width = parseFloat($("branch-width").value);
      ann = setBranchColor(ann, this._selectedNodeId, color);
      ann = setBranchWidth(ann, this._selectedNodeId, width);
      this._emit(ann);
    });

    $("branch-apply-clade").addEventListener("click", () => {
      if (!this._selectedNode || !this._root) return;
      const color = $("branch-color").value;
      const width = parseFloat($("branch-width").value);
      let ann = this._ann;
      const applyToSubtree = node => {
        ann = setBranchColor(ann, node.id, color);
        ann = setBranchWidth(ann, node.id, width);
        for (const c of node.children) applyToSubtree(c);
      };
      applyToSubtree(this._selectedNode);
      this._emit(ann);
    });

    $("branch-color-clear").addEventListener("click", () => {
      if (!this._selectedNodeId) return;
      this._emit(setBranchColor(this._ann, this._selectedNodeId, null));
    });

    // ── Display options ──────────────────────────────────────
    $("opt-show-support").addEventListener("change", e => {
      this._emit(this._ann, { showSupport: e.target.checked });
    });
    $("opt-radial").addEventListener("change", e => {
      this._emit(this._ann, { radial: e.target.checked });
    });
    $("opt-font-size").addEventListener("input", e => {
      $("opt-font-size-val").textContent = e.target.value + " px";
      this._emit(this._ann, { fontSize: parseInt(e.target.value) });
    });

    // ── Metadata CSV import ──────────────────────────────────
    $("meta-import-btn").addEventListener("click", () => {
      $("meta-file-input").click();
    });
    $("meta-file-input").addEventListener("change", async e => {
      const file = e.target.files[0];
      if (!file) return;
      $("meta-file-name").textContent = file.name;
      const text  = await file.text();
      const sep   = file.name.endsWith(".tsv") ? "\t" : ",";
      // Use tab-separated if TSV
      const strip = parseMetadataCsv(file.name.replace(/\.[^.]+$/, ""), text, sep);
      let ann = upsertMetadataStrip(this._ann, strip);
      // Auto-update legend from strip
      ann = setLegend(ann, [
        ...(ann.legend || []),
        ...strip.legend.filter(e => !(ann.legend || []).find(l => l.label === e.label)),
      ]);
      this._emit(ann);
    });

    // ── Legend clear ─────────────────────────────────────────
    $("legend-clear-btn").addEventListener("click", () => {
      this._emit(setLegend(this._ann, []));
    });

    // ── Add shapes ───────────────────────────────────────────
    $("shape-add-text").addEventListener("click", () => {
      const text = prompt("Label text:");
      if (!text) return;
      this._emit(addShape(this._ann, { type: "text", x: 60, y: 60, text, color: "#00e5c8" }));
    });
    $("shape-add-rect").addEventListener("click", () => {
      this._emit(addShape(this._ann, { type: "rect", x: 60, y: 60, w: 80, h: 40, color: "#00e5c8" }));
    });
    $("shape-add-line").addEventListener("click", () => {
      this._emit(addShape(this._ann, { type: "line", x: 60, y: 60, x2: 140, y2: 60, color: "#00e5c8" }));
    });
  }

  // ── Sync UI state ────────────────────────────────────────────────────

  _syncUI() {
    this._syncNodeSection();
    this._syncBranchSection();
    this._syncMetaList();
    this._syncLegend();
    this._syncShapesList();
  }

  _syncNodeSection() {
    const node = this._selectedNode;
    document.getElementById("sel-hint").style.display    = node ? "none" : "block";
    document.getElementById("sel-controls").style.display = node ? "block" : "none";
    if (!node) return;
    document.getElementById("sel-name").textContent =
      (node.name || `node:${node.id}`) + (node.children.length ? ` [${node.children.length} children]` : "");
    const existingColor = this._ann?.nodeColors?.get(node.name || node.id);
    if (existingColor) document.getElementById("sel-node-color").value = existingColor;
    const existingLabel = this._ann?.labelOverrides?.get(node.id);
    document.getElementById("sel-label-input").value = existingLabel || "";
    const cg = (this._ann?.cladeGroups || []).find(g => g.nodeId === node.id);
    document.getElementById("clade-label-input").value = cg?.label || "";
    if (cg?.color) document.getElementById("clade-color").value = cg.color;
    document.getElementById("clade-collapse-btn").textContent =
      cg?.collapsed ? "Expand" : "Collapse";
  }

  _syncBranchSection() {
    const node = this._selectedNode;
    document.getElementById("branch-hint").style.display     = node ? "none" : "block";
    document.getElementById("branch-controls").style.display = node ? "block" : "none";
    if (!node) return;
    const bc = this._ann?.branchColors?.get(node.id);
    if (bc) document.getElementById("branch-color").value = bc;
    const bw = this._ann?.branchWidths?.get(node.id);
    if (bw) {
      document.getElementById("branch-width").value     = bw;
      document.getElementById("branch-width-val").textContent = bw + " px";
    }
  }

  _syncMetaList() {
    const el = document.getElementById("meta-strip-list");
    if (!el) return;
    const strips = this._ann?.metadataStrips || [];
    el.innerHTML = strips.map((s, i) =>
      `<div class="vopt-row" style="margin-bottom:0.2rem">
         <span class="vopt-mono" style="flex:1">${s.name}</span>
         <button class="vopt-btn-sm danger" data-strip-idx="${i}">✕</button>
       </div>`
    ).join("");
    el.querySelectorAll("[data-strip-idx]").forEach(btn => {
      btn.addEventListener("click", () => {
        const idx = parseInt(btn.dataset.stripIdx);
        const name = (this._ann?.metadataStrips || [])[idx]?.name;
        if (name) this._emit(removeMetadataStrip(this._ann, name));
      });
    });
  }

  _syncLegend() {
    const el = document.getElementById("legend-list");
    if (!el) return;
    const legend = this._ann?.legend || [];
    el.innerHTML = legend.map(e =>
      `<div class="vopt-row" style="margin-bottom:0.2rem;align-items:center;gap:0.4rem">
         <span style="display:inline-block;width:10px;height:10px;border-radius:50%;
               background:${e.color};flex-shrink:0"></span>
         <span style="font-size:0.75rem;flex:1">${e.label}</span>
       </div>`
    ).join("");
  }

  _syncShapesList() {
    const el = document.getElementById("shapes-list");
    if (!el) return;
    const shapes = this._ann?.shapes || [];
    el.innerHTML = shapes.map(s =>
      `<div class="vopt-row" style="margin-bottom:0.2rem">
         <span class="vopt-mono" style="flex:1">${s.type}${s.text ? ": " + s.text : ""}</span>
       </div>`
    ).join("");
  }

  _openSection(id) {
    const el = document.getElementById(id);
    if (el) el.open = true;
  }

  /**
   * Emit annotation change.
   * @param {object} ann      - new annotation state
   * @param {object} optsPatch - optional opts patch (radial, fontSize, etc.)
   */
  _emit(ann, optsPatch = null) {
    this._ann = ann;
    this._syncUI();
    this._onChange(ann, optsPatch);
  }
}
