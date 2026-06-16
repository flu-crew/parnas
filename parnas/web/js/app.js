/**
 * PARNAS web UI — app orchestrator.
 * Replaces the inline phylotree.js rendering in index.html.
 *
 * Responsibilities:
 *  - Wire up form submission → /api/run
 *  - Parse Newick response → tree model
 *  - Compute layout
 *  - Drive TreeRenderer (paper.js canvas)
 *  - Drive AnnotationPanel (Visual Options drawer)
 *  - Drive export dialog
 */

import { parseNewick }                        from "./newick.js";
import { leaves, indexByName }                from "./treemodel.js";
import { rectangularLayout }                  from "./layout.js";
import { TreeRenderer }                       from "./renderer.js";
import {
  emptyAnnotations, setClusterColors, setLegend, deserialise,
}                                             from "./annotations.js";
import {
  buildAnnotatedNexus, buildSessionJson,
  downloadBlob, downloadText,
  svgToRasterBlob, canvasToEps,
}                                             from "./export.js";
import { SweepChart }                        from "./sweep.js";

// ── State ──────────────────────────────────────────────────────────────────
let treeFile    = null;
let weightsFile = null;
let alnFile     = null;

let lastNewick  = null;   // raw newick string from /api/run
let lastData    = null;   // full /api/run JSON response
let treeRoot    = null;   // parsed tree root
let treeCoords  = null;   // Map<id, {x,y}>

let annotations = emptyAnnotations();
let renderOpts  = { dark: true, radial: false, fontSize: 11, showSupport: false };

/** @type {TreeRenderer|null} */
let renderer    = null;
/** @type {ResizeObserver|null} */
let resizeObs   = null;
/** @type {SweepChart|null} */
let sweepChart  = null;
let sweepAbort  = null;
let currentSweepN = 1;

// ── DOM refs ───────────────────────────────────────────────────────────────
const treeCard     = document.getElementById("tree-card");
const stateOverlay = document.getElementById("state-overlay");
const runBtn       = document.getElementById("run-btn");
const zoomControls = document.getElementById("zoom-controls");
const exportBtn    = document.getElementById("export-btn");
const exportDialog = document.getElementById("export-dialog");

// ── Initialise ─────────────────────────────────────────────────────────────

function init() {
  // Zoom buttons
  document.getElementById("zoom-in-btn")?.addEventListener("click",    () => renderer?.zoom(1.25));
  document.getElementById("zoom-out-btn")?.addEventListener("click",   () => renderer?.zoom(1 / 1.25));
  document.getElementById("zoom-reset-btn")?.addEventListener("click", () => renderer?.resetView());

  // Export button opens dialog
  exportBtn?.addEventListener("click", openExportDialog);

  // Theme toggle (re-render on switch)
  document.getElementById("theme-toggle")?.addEventListener("click", () => {
    setTimeout(() => {
      renderOpts.dark = document.documentElement.getAttribute("data-theme") !== "light";
      if (renderer && treeRoot) renderer.setAnnotations(annotations);
    }, 30);
  });

  // Sweep chart
  const sweepCanvas = document.getElementById("sweep-chart");
  if (sweepCanvas) sweepChart = new SweepChart(sweepCanvas);

  wireRunButton();
  wireExportDialog();
  wireFilePickers();
  wireDependentFields();
  wireSweepBar();
  restoreTheme();
}

// ── Run button ──────────────────────────────────────────────────────────────

let currentAbort = null;

function wireRunButton() {
  runBtn.addEventListener("click", async () => {
    if (currentAbort) { currentAbort.abort(); return; }

    if (!treeFile) {
      showOverlay("error", "⚠️", "Please upload a tree file first.");
      highlightField(document.getElementById("drop-zone"));
      return;
    }

    const coverCb    = document.getElementById("cover-cb");
    const evaluateCb = document.getElementById("evaluate-cb");
    const nInput     = document.getElementById("n-input");
    const cover      = coverCb.checked;
    const evaluate   = evaluateCb.checked;
    const n          = parseInt(nInput.value, 10);

    if (!cover && !evaluate && (!n || n < 1)) {
      showOverlay("error", "⚠️", "Please enter a valid number of representatives (≥ 1).");
      highlightField(nInput);
      return;
    }
    if (evaluate && !document.getElementById("prior-input").value.trim()) {
      showOverlay("error", "⚠️", "Evaluate mode requires a prior regex (--prior).");
      highlightField(document.getElementById("prior-input"));
      return;
    }

    currentAbort = new AbortController();
    runBtn.setAttribute("aria-busy", "true");
    runBtn.textContent = "Cancel";
    showOverlay("loading", "⏳", "Running PARNAS analysis…");

    const fd = new FormData();
    fd.append("tree",            treeFile);
    fd.append("n",               n || 1);
    fd.append("cover",           cover    ? "true" : "false");
    fd.append("evaluate",        evaluate ? "true" : "false");
    fd.append("binary",          document.getElementById("binary-cb").checked ? "true" : "false");
    fd.append("prior",           document.getElementById("prior-input").value.trim());
    fd.append("radius",          document.getElementById("radius-input").value.trim());
    fd.append("threshold",       document.getElementById("threshold-input").value.trim());
    fd.append("exclude_rep",     document.getElementById("exc-rep").value.trim());
    fd.append("exclude_obj",     document.getElementById("exc-obj").value.trim());
    fd.append("exclude_fully",   document.getElementById("exc-full").value.trim());
    fd.append("constrain_fully", document.getElementById("constrain").value.trim());
    fd.append("aln_type", document.querySelector("input[name='aln-type']:checked")?.value || "nt");
    if (weightsFile) fd.append("weights",   weightsFile);
    if (alnFile && parseFloat(document.getElementById("threshold-input").value) > 0)
      fd.append("alignment", alnFile);

    try {
      const resp = await fetch("/api/run", {
        method: "POST", body: fd, signal: currentAbort.signal,
      });
      const data = await resp.json();
      if (!resp.ok) throw new Error(data.error || "Server error");

      lastData   = data;
      lastNewick = data.tree;

      // Parse + layout + render
      annotations = emptyAnnotations();
      if (data.colors) annotations = setClusterColors(annotations, data.colors);

      // Build legend from PARNAS cluster colours
      const repColors = data.colors || {};
      const repsAll   = [
        ...(data.representatives      || []),
        ...(data.best_representatives || []),
        ...(data.prior_centers        || []),
      ];
      const legendEntries = repsAll
        .filter(r => repColors[r])
        .map(r => ({ label: r, color: repColors[r] }));
      if (legendEntries.length) annotations = setLegend(annotations, legendEntries);

      renderOpts.dark = document.documentElement.getAttribute("data-theme") !== "light";

      loadAndRender(data);
      showResults(data);

      zoomControls.style.display = "flex";
      exportBtn.style.display    = "flex";

      // Show sweep bar only for sample mode (not cover, not evaluate)
      const sweepBar = document.getElementById("sweep-bar");
      if (sweepBar) {
        const showSweep = data.mode === "sample" && !cover;
        sweepBar.style.display = showSweep ? "flex" : "none";
        if (showSweep) {
          currentSweepN = n || 1;
          document.getElementById("sweep-n-val").textContent = currentSweepN;
          const divText = data.diversity != null ? data.diversity.toFixed(2) + "% diversity" : "";
          document.getElementById("sweep-diversity-display").textContent = divText;
          if (sweepChart) {
            sweepChart.reset();
            sweepChart.addPoint(currentSweepN, data.diversity);
          }
        }
      }

    } catch (err) {
      if (err.name === "AbortError") {
        showOverlay("", "🛑", "Analysis cancelled.");
      } else {
        showOverlay("error", "✕", "Error: " + err.message);
        highlightField(fieldForError(err.message));
        console.error(err);
      }
    } finally {
      currentAbort = null;
      runBtn.removeAttribute("aria-busy");
      runBtn.textContent = "Run Analysis";
    }
  });
}

// ── Tree load & render ──────────────────────────────────────────────────────

function loadAndRender(data) {
  try {
    treeRoot = parseNewick(data.tree);
  } catch (err) {
    showOverlay("error", "✕", "Tree parse error: " + err.message);
    console.error(err);
    return;
  }


  computeLayout();
  mountRenderer();
  renderOpts.repsSet   = new Set([
    ...(data.representatives      || []),
    ...(data.best_representatives || []),
    ...(data.prior_centers        || []),
  ]);
  renderOpts.repColors = data.colors || {};

  const _nL = leaves(treeRoot).length;
  const _canvas = document.getElementById("tree-canvas");
  if (_canvas) applyCanvasHeight(_canvas, _nL);
  renderer.render(treeRoot, treeCoords, annotations, {
    ...renderOpts,
    fontSize: computeFontSize(_nL),
  });

  hideOverlay();
}

function recomputeLayout() {
  if (!treeRoot) return;
  computeLayout();
  const _nL2 = leaves(treeRoot).length;
  const _c2 = document.getElementById("tree-canvas");
  if (_c2) applyCanvasHeight(_c2, _nL2);
  renderer?.render(treeRoot, treeCoords, annotations, {
    ...renderOpts,
    fontSize: computeFontSize(_nL2),
  });
}

function computeLayout() {
  treeCoords = renderOpts.radial
    ? radialLayout(treeRoot)
    : rectangularLayout(treeRoot);
}

function computeFontSize(nLeaves) {
  if (nLeaves <= 50)  return 12;
  if (nLeaves <= 150) return 10;
  if (nLeaves <= 400) return 9;
  if (nLeaves <= 1500) return 8;
  return 7;
}

const PX_PER_LEAF = 16;

function canvasHeight(nLeaves) {
  return Math.max(treeCard.clientHeight || 600, nLeaves * PX_PER_LEAF);
}

function applyCanvasHeight(canvas, nLeaves) {
  const h = canvasHeight(nLeaves);
  canvas.width  = treeCard.clientWidth || canvas.width;
  canvas.height = h;
  canvas.style.height = h + "px";
  if (renderer) renderer.resize();
}

function mountRenderer() {
  // Clear tree card but keep the state overlay
  treeCard.innerHTML = "";
  treeCard.appendChild(stateOverlay);

  // Create canvas
  const canvas = document.createElement("canvas");
  canvas.id     = "tree-canvas";
  canvas.style.cssText = "display:block;width:100%;";
  canvas.width  = treeCard.clientWidth  || 800;
  canvas.height = treeCard.clientHeight || 600;
  treeCard.appendChild(canvas);

  // Resize observer
  if (window.ResizeObserver) {
    if (resizeObs) resizeObs.disconnect();
    resizeObs = new ResizeObserver(() => {
      const nL = treeRoot ? leaves(treeRoot).length : 0;
      canvas.width  = treeCard.clientWidth;
      if (nL) applyCanvasHeight(canvas, nL); else canvas.height = treeCard.clientHeight || 600;
      if (renderer && treeRoot) {
        renderer.render(treeRoot, treeCoords, annotations, {
          ...renderOpts,
          fontSize: computeFontSize(nL),
        });
      }
    });
    resizeObs.observe(treeCard);
  }

  if (renderer) renderer.destroy();
  renderer = new TreeRenderer(canvas);

}

// ── Results panel ────────────────────────────────────────────────────────────

function showResults(data) {
  const resultsPanel = document.getElementById("results-panel");
  resultsPanel.style.display = "block";

  const sampleDiv = document.getElementById("sample-results");
  const evalDiv   = document.getElementById("eval-results");

  if (data.mode === "sample") {
    sampleDiv.style.display = "block";
    evalDiv.style.display   = "none";

    const pill = document.getElementById("diversity-pill");
    if (data.diversity != null) {
      document.getElementById("diversity-val").textContent = data.diversity.toFixed(2) + "%";
      pill.style.display = "inline-flex";
    } else {
      pill.style.display = "none";
    }

    const reps = data.representatives || [];
    document.getElementById("rep-label").textContent = `Representatives (${reps.length})`;
    const list = document.getElementById("rep-list");
    list.innerHTML = "";
    reps.forEach(rep => {
      const li = document.createElement("li");
      li.className = "rep-item";
      const color = data.colors?.[rep] || "#888";
      li.innerHTML = `<span class="rep-dot" style="background:${color};color:${color}"></span>
                      <span>${escapeHtml(rep)}</span>`;
      list.appendChild(li);
    });

  } else {
    sampleDiv.style.display = "none";
    evalDiv.style.display   = "block";

    const statsEl = document.getElementById("eval-stats");
    statsEl.innerHTML = "";

    if (data.mode === "evaluate_cover") {
      statsEl.appendChild(makeStat("Taxa covered within radius", data.coverage_pct + "%"));
    } else {
      if (data.prior_diversity != null)
        statsEl.appendChild(makeStat("Prior diversity",   data.prior_diversity + "%"));
      if (data.best_diversity != null)
        statsEl.appendChild(makeStat("Best possible",     data.best_diversity  + "%"));
      statsEl.appendChild(makeStat("Better than random", data.percentile + "% of sets"));
    }

    const list = document.getElementById("eval-rep-list");
    list.innerHTML = "";
    const reps = data.best_representatives || data.prior_centers || [];
    reps.forEach(rep => {
      const li = document.createElement("li");
      li.className = "rep-item";
      const color = data.colors?.[rep] || "#888";
      li.innerHTML = `<span class="rep-dot" style="background:${color};color:${color}"></span>
                      <span>${escapeHtml(rep)}</span>`;
      list.appendChild(li);
    });
  }
}

function makeStat(label, value) {
  const div = document.createElement("div");
  div.className = "stat-item";
  div.innerHTML = `<span class="stat-label">${escapeHtml(label)}</span>
                   <span class="stat-value">${escapeHtml(String(value))}</span>`;
  return div;
}

// ── Sweep ────────────────────────────────────────────────────────────────────

function wireSweepBar() {
  document.getElementById("sweep-dec")?.addEventListener("click", () => sweepTo(currentSweepN - 1));
  document.getElementById("sweep-inc")?.addEventListener("click", () => sweepTo(currentSweepN + 1));
}

async function sweepTo(newN) {
  if (!treeFile || !treeRoot) return;
  if (newN < 1) return;
  // Cover mode: n is irrelevant
  if (document.getElementById("cover-cb")?.checked) return;

  // Cancel prior in-flight sweep
  if (sweepAbort) sweepAbort.abort();
  sweepAbort = new AbortController();

  currentSweepN = newN;
  document.getElementById("sweep-n-val").textContent = newN;
  document.getElementById("sweep-diversity-display").textContent = "…";

  const fd = new FormData();
  fd.append("tree",            treeFile);
  fd.append("n",               newN);
  fd.append("cover",           "false");
  fd.append("evaluate",        "false");
  fd.append("binary",          document.getElementById("binary-cb").checked ? "true" : "false");
  fd.append("prior",           document.getElementById("prior-input").value.trim());
  fd.append("radius",          document.getElementById("radius-input").value.trim());
  fd.append("threshold",       document.getElementById("threshold-input").value.trim());
  fd.append("exclude_rep",     document.getElementById("exc-rep").value.trim());
  fd.append("exclude_obj",     document.getElementById("exc-obj").value.trim());
  fd.append("exclude_fully",   document.getElementById("exc-full").value.trim());
  fd.append("constrain_fully", document.getElementById("constrain").value.trim());
  if (weightsFile) fd.append("weights", weightsFile);

  try {
    const resp = await fetch("/api/run", {
      method: "POST", body: fd, signal: sweepAbort.signal,
    });
    const data = await resp.json();
    if (!resp.ok) throw new Error(data.error || "Server error");

    // Update colors without re-layout
    annotations = emptyAnnotations();
    if (data.colors) annotations = setClusterColors(annotations, data.colors);
    const legendEntries = (data.representatives || [])
      .filter(r => data.colors?.[r])
      .map(r => ({ label: r, color: data.colors[r] }));
    if (legendEntries.length) annotations = setLegend(annotations, legendEntries);

    renderOpts.repsSet   = new Set(data.representatives || []);
    renderOpts.repColors = data.colors || {};
    if (renderer) renderer.setAnnotations(annotations);

    // Update sweep bar display
    const divText = data.diversity != null ? data.diversity.toFixed(2) + "% diversity" : "";
    document.getElementById("sweep-diversity-display").textContent = divText;
    if (sweepChart) sweepChart.addPoint(newN, data.diversity);

    // Sync panel if open

  } catch (err) {
    if (err.name !== "AbortError") {
      document.getElementById("sweep-diversity-display").textContent = "Error";
      console.error("Sweep failed:", err);
    }
  }
}

// ── Export ───────────────────────────────────────────────────────────────────

let exportFormat = "svg";
let exportTheme  = "dark";

function openExportDialog() {
  exportTheme = document.documentElement.getAttribute("data-theme") || "dark";
  document.querySelectorAll("#exp-theme-group .exp-fmt-btn").forEach(b => {
    b.classList.toggle("active", b.dataset.exptheme === exportTheme);
  });
  // Pre-fill height from actual canvas so export isn't clipped
  const cv = document.getElementById("tree-canvas");
  if (cv) {
    document.getElementById("exp-width").value  = cv.width  || 1200;
    document.getElementById("exp-height").value = cv.height || 800;
  }
  exportDialog.showModal();
}

function wireExportDialog() {
  document.getElementById("exp-cancel-btn")?.addEventListener("click", () => exportDialog.close());
  exportDialog.addEventListener("click", e => { if (e.target === exportDialog) exportDialog.close(); });

  document.querySelectorAll("[data-fmt]").forEach(btn => {
    btn.addEventListener("click", () => {
      document.querySelectorAll("[data-fmt]").forEach(b => b.classList.remove("active"));
      btn.classList.add("active");
      exportFormat = btn.dataset.fmt;
      const isSvg = exportFormat === "svg" || exportFormat === "nex";
      document.getElementById("exp-dpi-field").style.opacity = isSvg ? "0.4" : "1";
      document.getElementById("exp-dpi").disabled = isSvg;
      // Lock height field for raster: always use canvas height
      const hField = document.getElementById("exp-height");
      const cvEl   = document.getElementById("tree-canvas");
      if (!isSvg && cvEl) {
        hField.value    = cvEl.height;
        hField.readOnly = true;
        hField.style.opacity = "0.5";
      } else {
        hField.readOnly = false;
        hField.style.opacity = "1";
      }
    });
  });

  document.querySelectorAll("[data-exptheme]").forEach(btn => {
    btn.addEventListener("click", () => {
      document.querySelectorAll("[data-exptheme]").forEach(b => b.classList.remove("active"));
      btn.classList.add("active");
      exportTheme = btn.dataset.exptheme;
    });
  });

  document.getElementById("exp-download-btn")?.addEventListener("click", async () => {
    if (!treeRoot) return;
    const cv     = document.getElementById("tree-canvas");
    const width  = Math.max(100, parseInt(document.getElementById("exp-width").value)  || 1200);
    // Always use actual canvas pixel height for raster exports to avoid clipping
    const isRaster = exportFormat === "png" || exportFormat === "jpg" || exportFormat === "eps";
    const height = isRaster
      ? (cv?.height || Math.max(100, parseInt(document.getElementById("exp-height").value) || 800))
      : Math.max(100, parseInt(document.getElementById("exp-height").value) || 800);
    const dpi    = Math.max(72,  parseInt(document.getElementById("exp-dpi").value)    || 96);
    const bgColor = exportTheme === "light" ? "#f5f7fa" : "#0d1117";
    const dlBtn = document.getElementById("exp-download-btn");
    dlBtn.setAttribute("aria-busy", "true");
    dlBtn.textContent = "Exporting…";

    try {
      if (exportFormat === "nex") {
        const nexus = buildAnnotatedNexus(treeRoot, annotations, lastData || {});
        downloadText(nexus, "parnas-tree.nex", "text/plain");

      } else if (exportFormat === "svg") {
        const svgStr = renderer.exportSvgWithTheme(exportTheme !== "light");
        downloadBlob(new Blob([svgStr], { type: "image/svg+xml" }), "parnas-tree.svg");

      } else if (exportFormat === "png" || exportFormat === "jpg") {
        const canvasW = Math.round(width  * dpi / 96);
        const canvasH = Math.round(height * dpi / 96);
        const blob = await renderer.exportPngWithTheme(exportTheme !== "light", bgColor, canvasW, canvasH, exportFormat);
        downloadBlob(blob, `parnas-tree.${exportFormat}`);

      } else if (exportFormat === "eps") {
        const canvasW = Math.round(width  * dpi / 96);
        const canvasH = Math.round(height * dpi / 96);
        const pngBlob = await renderer.exportPngWithTheme(exportTheme !== "light", bgColor, canvasW, canvasH, "png");
        const imgBitmap = await createImageBitmap(pngBlob);
        const off = document.createElement("canvas");
        off.width = canvasW; off.height = canvasH;
        const ctx = off.getContext("2d");
        ctx.fillStyle = bgColor;
        ctx.fillRect(0, 0, canvasW, canvasH);
        ctx.drawImage(imgBitmap, 0, 0);
        const epsStr = canvasToEps(off, width, height);
        downloadText(epsStr, "parnas-tree.eps", "application/postscript");
      }

      exportDialog.close();
    } catch (err) {
      console.error("Export failed:", err);
      dlBtn.textContent = "Error — retry";
      return;
    } finally {
      dlBtn.removeAttribute("aria-busy");
      dlBtn.textContent = "Download";
    }
  });

  // Session save
  document.getElementById("exp-save-session-btn")?.addEventListener("click", () => {
    if (!treeRoot) return;
    const json = buildSessionJson(lastNewick, annotations, lastData || {});
    downloadBlob(new Blob([json], { type: "application/json" }), "parnas-session.json");
  });

  // Session load
  const loadSessionInput = document.getElementById("load-session-input");
  document.getElementById("load-session-btn")?.addEventListener("click", () => {
    loadSessionInput?.click();
  });
  loadSessionInput?.addEventListener("change", () => {
    const f = loadSessionInput.files[0];
    if (!f) return;
    const reader = new FileReader();
    reader.onload = e => {
      try {
        loadSession(JSON.parse(e.target.result));
      } catch (err) {
        showOverlay("error", "✕", "Session load failed: " + err.message);
        console.error(err);
      }
    };
    reader.readAsText(f);
    loadSessionInput.value = "";
  });
}

function loadSession(json) {
  if (!json?.version || !json?.newick) throw new Error("Invalid session file.");

  lastNewick = json.newick;
  lastData   = json.appState || {};

  try {
    treeRoot = parseNewick(json.newick);
  } catch (err) {
    throw new Error("Tree parse error: " + err.message);
  }

  annotations = json.annotations ? deserialise(json.annotations) : emptyAnnotations();


  renderOpts = {
    ...renderOpts,
    dark:     document.documentElement.getAttribute("data-theme") !== "light",
    repsSet:  new Set([
      ...(lastData.representatives      || []),
      ...(lastData.best_representatives || []),
      ...(lastData.prior_centers        || []),
    ]),
    repColors: lastData.colors || {},
  };

  computeLayout();
  mountRenderer();
  const _nL3 = leaves(treeRoot).length;
  const _c3 = document.getElementById("tree-canvas");
  if (_c3) applyCanvasHeight(_c3, _nL3);
  renderer.render(treeRoot, treeCoords, annotations, {
    ...renderOpts,
    fontSize: computeFontSize(_nL3),
  });

  hideOverlay();
  zoomControls.style.display = "flex";
  exportBtn.style.display    = "flex";

  if (lastData.mode) showResults(lastData);
}

// ── File pickers & dependent fields ──────────────────────────────────────────

function wireFilePickers() {
  const dropZone  = document.getElementById("drop-zone");
  const treeInput = document.getElementById("tree-file");
  const dzPrimary = document.getElementById("dz-primary");
  const dzSecondary = document.getElementById("dz-secondary");

  dropZone.addEventListener("click",   () => treeInput.click());
  dropZone.addEventListener("keydown", e => { if (e.key==="Enter"||e.key===" ") treeInput.click(); });
  dropZone.addEventListener("dragover",  e => { e.preventDefault(); dropZone.classList.add("drag-active"); });
  dropZone.addEventListener("dragleave", ()=> dropZone.classList.remove("drag-active"));
  dropZone.addEventListener("drop", e => {
    e.preventDefault(); dropZone.classList.remove("drag-active");
    if (e.dataTransfer.files[0]) setTreeFile(e.dataTransfer.files[0]);
  });
  treeInput.addEventListener("change", () => { if (treeInput.files[0]) setTreeFile(treeInput.files[0]); });

  function setTreeFile(f) {
    treeFile = f;
    dzPrimary.textContent   = f.name;
    dzSecondary.textContent = (f.size / 1024).toFixed(1) + " KB · click to change";
  }

  bindFilePicker("weights-btn", "weights-file", "weights-name", f => { weightsFile = f; });
  bindFilePicker("aln-btn",     "aln-file",     "aln-name",     f => { alnFile = f; });
}

function bindFilePicker(btnId, inputId, nameId, onPick) {
  const btn    = document.getElementById(btnId);
  const input  = document.getElementById(inputId);
  const nameEl = document.getElementById(nameId);
  if (!btn || !input) return;
  btn.addEventListener("click", () => input.click());
  input.addEventListener("change", () => {
    const f = input.files[0];
    if (!f) return;
    nameEl.textContent = f.name;
    nameEl.classList.remove("empty");
    onPick(f);
  });
}

function wireDependentFields() {
  const radiusInput    = document.getElementById("radius-input");
  const thresholdInput = document.getElementById("threshold-input");
  const coverCb        = document.getElementById("cover-cb");
  const binaryCb       = document.getElementById("binary-cb");
  const binaryRow      = document.getElementById("binary-row");
  const coverRow       = document.getElementById("cover-row");
  const nInput         = document.getElementById("n-input");
  const nRow           = document.getElementById("n-row");
  const alnSection     = document.getElementById("alignment-section");

  function update() {
    const radiusOk = parseFloat(radiusInput.value) > 0;
    const coverOk  = radiusOk || parseFloat(thresholdInput.value) > 0;
    binaryCb.disabled = !radiusOk;
    binaryRow.classList.toggle("dimmed", !radiusOk);
    if (!radiusOk) binaryCb.checked = false;
    coverCb.disabled = !coverOk;
    coverRow.classList.toggle("dimmed", !coverOk);
    if (!coverOk) coverCb.checked = false;
    const covered = coverCb.checked;
    nInput.disabled    = covered;
    nRow.style.opacity = covered ? "0.45" : "1";
    alnSection.style.display = parseFloat(thresholdInput.value) > 0 ? "block" : "none";
  }

  radiusInput?.addEventListener("input",    update);
  thresholdInput?.addEventListener("input", update);
  coverCb?.addEventListener("change",       update);
  update();
}

// ── Theme ────────────────────────────────────────────────────────────────────

function restoreTheme() {
  const saved = localStorage.getItem("parnas-theme") || "dark";
  applyTheme(saved);
  document.getElementById("theme-toggle")?.addEventListener("click", () => {
    const next = document.documentElement.getAttribute("data-theme") === "dark" ? "light" : "dark";
    applyTheme(next);
  });
}

function applyTheme(theme) {
  document.documentElement.setAttribute("data-theme", theme);
  const iconDark  = document.getElementById("icon-dark");
  const iconLight = document.getElementById("icon-light");
  if (iconDark)  iconDark.style.display  = theme === "dark"  ? "block" : "none";
  if (iconLight) iconLight.style.display = theme === "light" ? "block" : "none";
  localStorage.setItem("parnas-theme", theme);
  renderOpts.dark = theme !== "light";
}

// ── Overlay helpers ────────────────────────────────────────────────────────

function showOverlay(type, icon, msg) {
  const el = document.getElementById("state-overlay");
  el.className = "state-overlay" + (type === "error" ? " error" : "");
  el.innerHTML = `<span class="state-icon">${icon}</span><span>${msg}</span>`;
  el.style.display = "flex";
}

function hideOverlay() {
  const el = document.getElementById("state-overlay");
  if (el) el.style.display = "none";
}

// ── Field error helpers ────────────────────────────────────────────────────

function clearHighlights() {
  document.querySelectorAll(".field-error,.drop-zone-error")
    .forEach(el => el.classList.remove("field-error", "drop-zone-error"));
}

function highlightField(el) {
  if (!el) return;
  clearHighlights();
  const parent = el.closest("details.field-group");
  if (parent) parent.open = true;
  el.classList.add(el.id === "drop-zone" ? "drop-zone-error" : "field-error");
  el.scrollIntoView({ behavior: "smooth", block: "center" });
  ["input", "change", "click"].forEach(ev =>
    el.addEventListener(ev, clearHighlights, { once: true }));
}

function fieldForError(msg) {
  if (!msg) return null;
  const m = msg.toLowerCase();
  if (m.includes("tree file") || m.includes("parse tree")) return document.getElementById("drop-zone");
  if (m.includes("n must be"))  return document.getElementById("n-input");
  if (m.includes("radius"))     return document.getElementById("radius-input");
  if (m.includes("threshold"))  return document.getElementById("threshold-input");
  if (m.includes("alignment") || m.includes("treetime")) return document.getElementById("aln-btn");
  if (m.includes("weights"))    return document.getElementById("weights-btn");
  if (m.includes("prior") || m.includes("evaluate requires")) return document.getElementById("prior-input");
  return null;
}

// ── Misc helpers ─────────────────────────────────────────────────────────────

function escapeHtml(s) {
  return String(s).replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;");
}

// ── Boot ──────────────────────────────────────────────────────────────────────

document.addEventListener("DOMContentLoaded", init);
