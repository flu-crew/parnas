/**
 * PARNAS web UI — app orchestrator.
 *
 * Responsibilities:
 *  - Wire up form submission → /api/run
 *  - Parse Newick response → tree model
 *  - Compute layout
 *  - Drive TreeRenderer (paper.js canvas)
 *  - Drive export dialog
 */

import { parseNewick }                        from "./newick.js";
import { leaves }                             from "./treemodel.js";
import { rectangularLayout, cladogramLayout, radialLayout } from "./layout.js";
import { TreeRenderer }                       from "./renderer.js";
import {
  emptyAnnotations, setClusterColors, setLegend, setLegendCaption, deserialise,
}                                             from "./annotations.js";
import {
  buildAnnotatedNexus, buildSessionJson,
  downloadBlob, downloadText,
  canvasToEps,
}                                             from "./export.js";
import { SweepChart }                        from "./sweep.js";

// ── Utilities ──────────────────────────────────────────────────────────────
function debounce(fn, ms) {
  let t;
  return (...args) => { clearTimeout(t); t = setTimeout(() => fn(...args), ms); };
}

// ── State ──────────────────────────────────────────────────────────────────
let treeFile    = null;
let weightsFile = null;
let alnFile     = null;

let lastNewick  = null;   // raw newick string from /api/run
let lastData    = null;   // full /api/run JSON response
let treeRoot    = null;   // parsed tree root
let treeCoords  = null;   // Map<id, {x,y}>

let annotations      = emptyAnnotations();
let annotationsPrior = emptyAnnotations();
let annotationsBest  = emptyAnnotations();
let renderOpts  = { dark: true, radial: false, layout: "phylogram", fontSize: 11, showSupport: false };

/** @type {TreeRenderer|null} */
let renderer    = null;
/** @type {TreeRenderer|null} — evaluate mode left panel (prior reps) */
let rendererPrior = null;
/** @type {TreeRenderer|null} — evaluate mode right panel (best reps) */
let rendererBest  = null;
/** @type {ResizeObserver|null} */
let resizeObs   = null;
/** @type {SweepChart|null} */
let sweepChart  = null;
let sweepAbort  = null;
let currentSweepN = 1;
let lastNTaxa   = null; // from most recent /api/run response
let lastReps    = [];   // representatives from most recent analysis or sweep
/** Session-only memo cache: param-key → /api/run response. Not serialized in export. */
const sweepCache = new Map();
const SWEEP_CACHE_MAX = 20;
/** LRU-insert into sweepCache; evicts oldest entry when over cap. */
function _sweepCacheSet(key, value) {
  if (sweepCache.has(key)) sweepCache.delete(key); // move to end (LRU)
  sweepCache.set(key, value);
  if (sweepCache.size > SWEEP_CACHE_MAX) sweepCache.delete(sweepCache.keys().next().value);
}

// ── DOM refs ───────────────────────────────────────────────────────────────
const treeCard     = document.getElementById("tree-card");
const stateOverlay = document.getElementById("state-overlay");
const runBtn       = document.getElementById("run-btn");
const exportBtn     = document.getElementById("export-btn");
const layoutToggle  = document.getElementById("layout-toggle");
const exportDialog = document.getElementById("export-dialog");

// ── Initialise ─────────────────────────────────────────────────────────────

function init() {
  // Generic section collapse buttons (data-target → id of collapsible div)
  document.querySelectorAll(".collapse-btn[data-target]").forEach(btn => {
    btn.addEventListener("click", () => {
      const el = document.getElementById(btn.dataset.target);
      if (!el) return;
      const collapsed = el.classList.toggle("collapsed");
      btn.textContent = collapsed ? "▸" : "▾";
    });
  });

  document.getElementById("rep-export-btn")?.addEventListener("click", () => {
    if (!lastReps.length) return;
    downloadText(lastReps.join("\n") + "\n", "parnas-representatives.txt", "text/plain");
  });

  // Sidebar tab switching
  function activateTab(name) {
    document.querySelectorAll(".sb-tab").forEach(b => b.classList.toggle("active", b.dataset.tab === name));
    document.querySelectorAll(".sb-panel").forEach(p => p.classList.toggle("active", p.id === "tab-" + name));
  }
  document.querySelectorAll(".sb-tab").forEach(btn => {
    btn.addEventListener("click", () => activateTab(btn.dataset.tab));
  });
  window._activateSidebarTab = activateTab;

  // Sidebar toggle (toolbar button - reopen)
  document.getElementById("sidebar-toggle")?.addEventListener("click", () => {
    const body = document.querySelector(".app-body");
    const btn  = document.getElementById("sidebar-toggle");
    const collapsed = body.classList.toggle("sidebar-collapsed");
    btn.textContent = collapsed ? "▶" : "◀";
    setTimeout(() => renderer?.resize(), 220);
  });

  // Export button opens dialog
  exportBtn?.addEventListener("click", () => {
    if (rendererPrior && rendererBest) exportEvaluateImage();
    else openExportDialog();
  });

  // Layout toggle (Phylogram ⇄ Cladogram)
  document.getElementById("layout-toggle")?.addEventListener("click", () => {
    renderOpts.layout = renderOpts.layout === "cladogram" ? "phylogram" : "cladogram";
    const btn = document.getElementById("layout-toggle");
    if (btn) btn.textContent = renderOpts.layout === "cladogram" ? "Phylogram" : "Cladogram";
    rerenderActive();
  });

  // Reset tree view button
  document.getElementById("reset-view-btn")?.addEventListener("click", () => {
    renderer?.resetView();
    rendererPrior?.resetView();
    rendererBest?.resetView();
  });

  // Sweep chart
  const sweepCanvas = document.getElementById("sweep-chart");
  if (sweepCanvas) {
    sweepChart = new SweepChart(sweepCanvas);

    // Hover tooltip on sweep canvas
    const sweepTooltip = document.getElementById("sweep-tooltip");
    sweepCanvas.addEventListener("mousemove", e => {
      if (!sweepChart || !sweepTooltip) return;
      const br = sweepCanvas.getBoundingClientRect();
      const cssX = e.clientX - br.left;
      const cssY = e.clientY - br.top;
      const hit = sweepChart.hitTest(cssX, cssY);
      if (hit) {
        sweepTooltip.textContent =
          `${sweepChart._xTitle}: ${sweepChart._xFmt(hit.x)}  ·  ${sweepChart._yTitle}: ${sweepChart._yFmt(hit.y)}`;
        sweepTooltip.style.left = hit.px + "px";
        sweepTooltip.style.top  = hit.py + "px";
        sweepTooltip.style.display = "block";
      } else {
        sweepTooltip.style.display = "none";
      }
    });
    sweepCanvas.addEventListener("mouseleave", () => {
      if (sweepTooltip) sweepTooltip.style.display = "none";
    });
  }

  wireRunButton();
  wireExportDialog();
  wireFilePickers();
  wireDependentFields();
  wireSweepBar();
  wireElbowBar();
  wireSweepToggles();
  wireSweepExport();
  restoreTheme();
}

// ── Run button ──────────────────────────────────────────────────────────────

let currentAbort = null;

function wireRunButton() {
  runBtn.addEventListener("click", async () => {
    if (currentAbort) { currentAbort.abort(); return; }

    // Delegate to sweep modes when active
    if (document.getElementById("sweep-n-cb")?.checked) {
      await runNSweep();
      return;
    }
    if (document.getElementById("sweep-thr-cb")?.checked) {
      await runThresholdSweep();
      return;
    }

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

    // For sample (non-cover, non-evaluate) runs, check shared sweep cache first
    const runCacheKey = (!cover && !evaluate) ? JSON.stringify({
      n:               n || 1,
      binary:          document.getElementById("binary-cb").checked,
      prior:           document.getElementById("prior-input").value.trim(),
      radius:          document.getElementById("radius-input").value.trim(),
      threshold:       document.getElementById("threshold-input").value.trim(),
      exclude_rep:     document.getElementById("exc-rep").value.trim(),
      exclude_obj:     document.getElementById("exc-obj").value.trim(),
      exclude_fully:   document.getElementById("exc-full").value.trim(),
      constrain_fully: document.getElementById("constrain").value.trim(),
      weightsKey:      weightsFile ? `${weightsFile.name}:${weightsFile.size}` : null,
    }) : null;

    if (runCacheKey && sweepCache.has(runCacheKey)) {
      const cached = sweepCache.get(runCacheKey);
      lastData   = cached;
      lastNewick = cached.tree;
      annotations = emptyAnnotations();
      if (cached.colors) annotations = setClusterColors(annotations, cached.colors);
      const repColors = cached.colors || {};
      const repsAll = [
        ...(cached.representatives      || []),
        ...(cached.best_representatives || []),
        ...(cached.prior_centers        || []),
      ];
      const legendEntries = repsAll.filter(r => repColors[r]).map(r => ({ label: r, color: repColors[r] }));
      if (legendEntries.length) annotations = setLegend(annotations, legendEntries);
      renderOpts.dark = document.documentElement.getAttribute("data-theme") !== "light";
      loadAndRender(cached);
      showResults(cached);
      exportBtn.style.display = "flex";
      if (layoutToggle) layoutToggle.style.display = "flex";
      document.getElementById("reset-view-btn").style.display = "flex";
      const sweepSection = document.getElementById("sweep-section");
      if (sweepSection) sweepSection.style.display = "";
      currentSweepN = n || 1;
      document.getElementById("sweep-n-val").textContent = currentSweepN;
      const divText = cached.diversity != null ? cached.diversity.toFixed(2) + "% diversity" : "";
      document.getElementById("sweep-diversity-display").textContent = divText;
      if (sweepChart) { sweepChart.reset(); sweepChart.addPoint(currentSweepN, cached.diversity); }
      if (cached.n_taxa != null) lastNTaxa = cached.n_taxa;
      currentAbort = null;
      runBtn.removeAttribute("aria-busy");
      runBtn.textContent = "Run Analysis";
      return;
    }

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

      // Cache sample results so sweep and repeat runs skip the server
      if (runCacheKey) _sweepCacheSet(runCacheKey, data);

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

      exportBtn.style.display = "flex";
      if (layoutToggle) layoutToggle.style.display = "flex";
      document.getElementById("reset-view-btn").style.display = "flex";

      // Sweep submenu: visible only for non-cover sample mode
      const showSweep = data.mode === "sample" && !cover;
      const sweepSection = document.getElementById("sweep-section");
      if (sweepSection) sweepSection.style.display = showSweep ? "" : "none";
      if (showSweep) {
        currentSweepN = n || 1;
        document.getElementById("sweep-n-val").textContent = currentSweepN;
        const divText = data.diversity != null ? data.diversity.toFixed(2) + "% diversity" : "";
        document.getElementById("sweep-diversity-display").textContent = divText;
        if (sweepChart) {
          sweepChart.reset();
          sweepChart.addPoint(currentSweepN, data.diversity);
        }
        if (data.n_taxa != null) lastNTaxa = data.n_taxa;
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
  const _nL = leaves(treeRoot).length;
  const _fs  = computeFontSize(_nL);

  if (data.mode === "evaluate") {
    // ── Dual split view ───────────────────────────────────────
    mountDualRenderer();

    const priorReps = data.prior_centers        || [];
    const bestReps  = data.best_representatives || [];

    annotationsPrior = emptyAnnotations();
    if (data.colors_prior) annotationsPrior = setClusterColors(annotationsPrior, data.colors_prior);
    const legendPrior = priorReps
      .filter(r => data.colors_prior?.[r])
      .map(r => ({ label: r, color: data.colors_prior[r] }));
    if (legendPrior.length) annotationsPrior = setLegend(annotationsPrior, legendPrior);
    if (data.prior_diversity != null)
      annotationsPrior = setLegendCaption(annotationsPrior, "Prior diversity: " + data.prior_diversity + "%");

    annotationsBest = emptyAnnotations();
    if (data.colors_best) annotationsBest = setClusterColors(annotationsBest, data.colors_best);
    const legendBest = bestReps
      .filter(r => data.colors_best?.[r])
      .map(r => ({ label: r, color: data.colors_best[r] }));
    if (legendBest.length) annotationsBest = setLegend(annotationsBest, legendBest);
    if (data.best_diversity != null)
      annotationsBest = setLegendCaption(annotationsBest, "Best possible: " + data.best_diversity + "%");

    renderOpts.repsSet   = new Set([...priorReps, ...bestReps]);
    renderOpts.repColors = {};

    const halfW = Math.floor(treeCard.clientWidth / 2) || 400;
    const hh    = canvasHeight(_nL);

    const cvPrior = document.getElementById("tree-canvas-prior");
    const cvBest  = document.getElementById("tree-canvas-best");
    if (cvPrior) { cvPrior.width = halfW; cvPrior.height = hh; cvPrior.style.height = hh + "px"; rendererPrior.resize(); }
    if (cvBest)  { cvBest.width  = halfW; cvBest.height  = hh; cvBest.style.height  = hh + "px"; rendererBest.resize();  }

    rendererPrior.render(treeRoot, treeCoords, annotationsPrior, { ...renderOpts, fontSize: _fs });
    rendererBest.render(treeRoot, treeCoords, annotationsBest,   { ...renderOpts, fontSize: _fs });

  } else {
    // ── Single canvas ─────────────────────────────────────────
    if (rendererPrior) { rendererPrior.destroy(); rendererPrior = null; }
    if (rendererBest)  { rendererBest.destroy();  rendererBest  = null; }
    mountRenderer();

    renderOpts.repsSet   = new Set([
      ...(data.representatives      || []),
      ...(data.best_representatives || []),
      ...(data.prior_centers        || []),
    ]);
    renderOpts.repColors = data.colors || {};

    const _canvas = document.getElementById("tree-canvas");
    if (_canvas) applyCanvasHeight(_canvas, _nL);
    renderer.render(treeRoot, treeCoords, annotations, {
      ...renderOpts,
      fontSize: _fs,
    });
  }

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
  if (renderOpts.radial) {
    treeCoords = radialLayout(treeRoot);
  } else if (renderOpts.layout === "cladogram") {
    treeCoords = cladogramLayout(treeRoot);
  } else {
    treeCoords = rectangularLayout(treeRoot);
  }
}

/** Re-render whichever view is currently mounted. */
function rerenderActive() {
  if (!treeRoot) return;
  computeLayout();
  const _nL = leaves(treeRoot).length;
  const _fs  = computeFontSize(_nL);
  if (rendererPrior && rendererBest) {
    rendererPrior.render(treeRoot, treeCoords, annotationsPrior, { ...renderOpts, fontSize: _fs });
    rendererBest.render(treeRoot,  treeCoords, annotationsBest,  { ...renderOpts, fontSize: _fs });
  } else if (renderer) {
    const _c = document.getElementById("tree-canvas");
    if (_c) applyCanvasHeight(_c, _nL);
    renderer.render(treeRoot, treeCoords, annotations, { ...renderOpts, fontSize: _fs });
  }
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
    resizeObs = new ResizeObserver(debounce(() => {
      const nL = treeRoot ? leaves(treeRoot).length : 0;
      canvas.width  = treeCard.clientWidth;
      if (nL) applyCanvasHeight(canvas, nL); else canvas.height = treeCard.clientHeight || 600;
      if (renderer && treeRoot) {
        renderer.render(treeRoot, treeCoords, annotations, {
          ...renderOpts,
          fontSize: computeFontSize(nL),
        });
      }
    }, 140));
    resizeObs.observe(treeCard);
  }

  if (renderer) renderer.destroy();
  renderer = new TreeRenderer(canvas);

}

function mountDualRenderer() {
  treeCard.innerHTML = "";
  treeCard.appendChild(stateOverlay);

  // Destroy old renderers
  if (renderer)      { renderer.destroy();      renderer      = null; }
  if (rendererPrior) { rendererPrior.destroy();  rendererPrior = null; }
  if (rendererBest)  { rendererBest.destroy();   rendererBest  = null; }

  const halfW = Math.floor(treeCard.clientWidth / 2) || 400;
  const h     = treeCard.clientHeight || 600;

  const splitView = document.createElement("div");
  splitView.id = "tree-split-view";

  const makePanel = (id, label) => {
    const panel = document.createElement("div");
    panel.className = "tree-panel";
    const lbl = document.createElement("div");
    lbl.className = "tree-panel-label";
    lbl.textContent = label;
    const cv = document.createElement("canvas");
    cv.id = id;
    cv.style.cssText = "display:block;width:100%;";
    cv.width  = halfW;
    cv.height = h;
    panel.appendChild(lbl);
    panel.appendChild(cv);
    return { panel, cv };
  };

  const { panel: panelPrior, cv: cvPrior } = makePanel("tree-canvas-prior", "Prior Representatives");
  const { panel: panelBest,  cv: cvBest  } = makePanel("tree-canvas-best",  "Best Coverage");

  splitView.appendChild(panelPrior);
  splitView.appendChild(panelBest);
  treeCard.appendChild(splitView);

  rendererPrior = new TreeRenderer(cvPrior);
  rendererBest  = new TreeRenderer(cvBest);

  if (window.ResizeObserver) {
    if (resizeObs) resizeObs.disconnect();
    resizeObs = new ResizeObserver(debounce(() => {
      if (!treeRoot) return;
      const nL  = leaves(treeRoot).length;
      const w2  = Math.floor(treeCard.clientWidth / 2);
      const hh  = canvasHeight(nL);
      const fs  = computeFontSize(nL);
      [cvPrior, cvBest].forEach(c => { c.width = w2; c.height = hh; c.style.height = hh + "px"; });
      rendererPrior.resize();
      rendererBest.resize();
      rendererPrior.render(treeRoot, treeCoords, annotationsPrior, { ...renderOpts, fontSize: fs });
      rendererBest.render(treeRoot, treeCoords, annotationsBest,   { ...renderOpts, fontSize: fs });
    }, 140));
    resizeObs.observe(treeCard);
  }
}

// ── Results panel ────────────────────────────────────────────────────────────

function showResults(data) {
  // Reveal the Results tab and switch to it
  const tabBtn = document.getElementById("tab-btn-results");
  if (tabBtn) tabBtn.style.display = "";
  if (window._activateSidebarTab) window._activateSidebarTab("results");
  // Ensure sidebar is open
  const body = document.querySelector(".app-body");
  if (body?.classList.contains("sidebar-collapsed")) {
    body.classList.remove("sidebar-collapsed");
    setTimeout(() => renderer?.resize(), 220);
  }

  const sampleDiv = document.getElementById("sample-results");
  const evalDiv   = document.getElementById("eval-results");

  if (data.mode === "sample") {
    sampleDiv.style.display = "block";
    evalDiv.style.display   = "none";
    renderSampleReps(data);

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

    // Prior representatives (evaluate mode only)
    const priorRepsData = data.prior_centers || [];
    const priorLabel = document.getElementById("eval-prior-label");
    const priorList  = document.getElementById("eval-prior-list");
    const showPrior  = data.mode !== "evaluate_cover" && priorRepsData.length > 0;
    if (priorLabel) priorLabel.style.display = showPrior ? "" : "none";
    if (priorList) {
      priorList.style.display = showPrior ? "" : "none";
      priorList.innerHTML = "";
      if (showPrior) {
        priorRepsData.forEach(rep => {
          const li = document.createElement("li");
          li.className = "rep-item";
          const color = data.colors_prior?.[rep] || data.colors?.[rep] || "#888";
          li.innerHTML = `<span class="rep-dot" style="background:${color};color:${color}"></span>
                          <span>${escapeHtml(rep)}</span>`;
          priorList.appendChild(li);
        });
      }
    }

    const list = document.getElementById("eval-rep-list");
    list.innerHTML = "";
    const bestReps = data.best_representatives || data.prior_centers || [];
    bestReps.forEach(rep => {
      const li = document.createElement("li");
      li.className = "rep-item";
      const color = data.colors_best?.[rep] || data.colors?.[rep] || "#888";
      li.innerHTML = `<span class="rep-dot" style="background:${color};color:${color}"></span>
                      <span>${escapeHtml(rep)}</span>`;
      list.appendChild(li);
    });
  }
}

function renderSampleReps(data) {
  const pill = document.getElementById("diversity-pill");
  if (data.diversity != null) {
    document.getElementById("diversity-val").textContent = data.diversity.toFixed(2) + "%";
    pill.style.display = "inline-flex";
  } else {
    pill.style.display = "none";
  }

  const reps = data.representatives || [];
  lastReps = reps;
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

function _applySweepResult(newN, data) {
  annotations = emptyAnnotations();
  if (data.colors) annotations = setClusterColors(annotations, data.colors);
  const legendEntries = (data.representatives || [])
    .filter(r => data.colors?.[r])
    .map(r => ({ label: r, color: data.colors[r] }));
  if (legendEntries.length) annotations = setLegend(annotations, legendEntries);

  renderOpts.repsSet   = new Set(data.representatives || []);
  renderOpts.repColors = data.colors || {};
  if (renderer) renderer.setAnnotations(annotations);

  const divText = data.diversity != null ? data.diversity.toFixed(2) + "% diversity" : "";
  document.getElementById("sweep-diversity-display").textContent = divText;
  if (sweepChart) sweepChart.addPoint(newN, data.diversity);
  renderSampleReps(data);
}

/** Build the FormData + cacheKey for a sample-mode request at a given n. */
function _buildSampleRequest(n) {
  const params = {
    n,
    binary:          document.getElementById("binary-cb").checked,
    prior:           document.getElementById("prior-input").value.trim(),
    radius:          document.getElementById("radius-input").value.trim(),
    threshold:       document.getElementById("threshold-input").value.trim(),
    exclude_rep:     document.getElementById("exc-rep").value.trim(),
    exclude_obj:     document.getElementById("exc-obj").value.trim(),
    exclude_fully:   document.getElementById("exc-full").value.trim(),
    constrain_fully: document.getElementById("constrain").value.trim(),
    weightsKey:      weightsFile ? `${weightsFile.name}:${weightsFile.size}` : null,
  };
  const cacheKey = JSON.stringify(params);
  const fd = new FormData();
  fd.append("tree",            treeFile);
  fd.append("n",               n);
  fd.append("cover",           "false");
  fd.append("evaluate",        "false");
  fd.append("binary",          params.binary ? "true" : "false");
  fd.append("prior",           params.prior);
  fd.append("radius",          params.radius);
  fd.append("threshold",       params.threshold);
  fd.append("exclude_rep",     params.exclude_rep);
  fd.append("exclude_obj",     params.exclude_obj);
  fd.append("exclude_fully",   params.exclude_fully);
  fd.append("constrain_fully", params.constrain_fully);
  if (weightsFile) fd.append("weights", weightsFile);
  return { fd, cacheKey };
}

async function sweepTo(newN) {
  if (!treeFile || !treeRoot) return;
  if (newN < 1) return;
  // Cover mode: n is irrelevant
  if (document.getElementById("cover-cb")?.checked) return;

  currentSweepN = newN;
  document.getElementById("sweep-n-val").textContent = newN;

  const { fd, cacheKey } = _buildSampleRequest(newN);

  if (sweepCache.has(cacheKey)) {
    _applySweepResult(newN, sweepCache.get(cacheKey));
    return;
  }

  // Cancel prior in-flight sweep
  if (sweepAbort) sweepAbort.abort();
  sweepAbort = new AbortController();

  document.getElementById("sweep-diversity-display").textContent = "…";

  try {
    const resp = await fetch("/api/run", {
      method: "POST", body: fd, signal: sweepAbort.signal,
    });
    const data = await resp.json();
    if (!resp.ok) throw new Error(data.error || "Server error");

    _sweepCacheSet(cacheKey, data);
    _applySweepResult(newN, data);

  } catch (err) {
    if (err.name !== "AbortError") {
      document.getElementById("sweep-diversity-display").textContent = "Error";
      console.error("Sweep failed:", err);
    }
  }
}

// ── Elbow sweep ──────────────────────────────────────────────────────────────

/**
 * Kneedle / max-triangle elbow detection on a concave-increasing diversity curve.
 * points: [{n, diversity}] sorted ascending by n.
 * Returns the n at the elbow, or null when there are too few points.
 */
/**
 * Generic elbow detector: max perpendicular distance from chord first→last.
 * Works for {x,y} points regardless of curve direction.
 * @param {Array<{x:number,y:number}>} points
 * @returns {number|null} x-value at the elbow
 */
function findElbow(points) {
  if (!points || points.length < 3) return points?.length ? points[points.length - 1].x : null;
  const first = points[0], last = points[points.length - 1];
  const dx = last.x - first.x, dy = last.y - first.y;
  const len = Math.sqrt(dx * dx + dy * dy);
  if (len === 0) return first.x;
  const a = dy / len, b = -dx / len;
  const c = -(a * first.x + b * first.y);
  let maxDist = -Infinity, elbowX = null;
  for (const p of points) {
    const dist = Math.abs(a * p.x + b * p.y + c);
    if (dist > maxDist) { maxDist = dist; elbowX = p.x; }
  }
  return elbowX;
}
// Legacy shim kept for backward compat with sweepTo path
function findElbowN(legacyPoints) {
  if (!legacyPoints) return null;
  return findElbow(legacyPoints.map(p => ({ x: p.n, y: p.diversity })));
}

// ── Sweep mode toggles ────────────────────────────────────────────────────

function wireSweepToggles() {
  const sweepNCb  = document.getElementById("sweep-n-cb");
  const sweepThrCb = document.getElementById("sweep-thr-cb");
  const sweepNRange  = document.getElementById("sweep-n-range");
  const sweepThrRange = document.getElementById("sweep-thr-range");

  function onSweepNChange() {
    if (sweepNCb?.checked && sweepThrCb) sweepThrCb.checked = false;
    if (sweepNRange) sweepNRange.style.display = sweepNCb?.checked ? "" : "none";
    if (sweepThrRange) sweepThrRange.style.display = "none";
    // Trigger alignment-section visibility update
    wireDependentFieldsUpdate();
  }
  function onSweepThrChange() {
    if (sweepThrCb?.checked && sweepNCb) sweepNCb.checked = false;
    if (sweepThrRange) sweepThrRange.style.display = sweepThrCb?.checked ? "" : "none";
    if (sweepNRange) sweepNRange.style.display = "none";
    // Show alignment section since threshold sweep always needs it
    const alnSection = document.getElementById("alignment-section");
    if (alnSection && sweepThrCb?.checked) alnSection.style.display = "block";
    wireDependentFieldsUpdate();
  }

  sweepNCb?.addEventListener("change",  onSweepNChange);
  sweepThrCb?.addEventListener("change", onSweepThrChange);
}

// External hook so wireSweepToggles can call wireDependentFields update logic
let wireDependentFieldsUpdate = () => {};

function wireElbowBar() {
  // No-op: elbow sweep now triggered via Run Analysis button
}

// ── n sweep (run via Run Analysis when sweep-n-cb is checked) ─────────────

async function runNSweep() {
  if (!treeFile) {
    showOverlay("error", "⚠️", "Upload a tree file first.");
    return;
  }

  let nMin = parseInt(document.getElementById("elbow-min")?.value, 10) || 2;
  let nMax = parseInt(document.getElementById("elbow-max")?.value, 10) || 20;
  nMin = Math.max(2, nMin);
  if (lastNTaxa != null) nMax = Math.min(nMax, lastNTaxa - 1);
  if (nMin >= nMax) {
    showOverlay("error", "⚠️", "From n must be less than To n (and ≥ 2).");
    return;
  }

  showOverlay("loading", "⏳", "Running n sweep…");
  runBtn.setAttribute("aria-busy", "true");
  runBtn.textContent = "Cancel";
  currentAbort = new AbortController();

  try {
    // One request at nMax — server returns diversity_curve[k-2] for k=2..nMax
    const { fd, cacheKey } = _buildSampleRequest(nMax);
    let data;
    if (sweepCache.has(cacheKey)) {
      data = sweepCache.get(cacheKey);
    } else {
      const resp = await fetch("/api/run", { method: "POST", body: fd, signal: currentAbort.signal });
      data = await resp.json();
      if (!resp.ok) throw new Error(data.error || "Server error");
      _sweepCacheSet(cacheKey, data);
    }
    if (data.n_taxa != null) lastNTaxa = data.n_taxa;
    if (!treeRoot && data.tree) treeRoot = parseNewick(data.tree);
    if (!renderer && treeRoot) loadAndRender(data);

    const curve = data.diversity_curve || [];
    const points = [];
    for (let k = nMin; k <= nMax; k++) {
      const val = curve[k - 2];
      if (val != null) points.push({ x: k, y: val });
    }
    if (points.length === 0) throw new Error("No diversity data in that n range.");

    const elbowX = findElbow(points);

    // Show figure — activate tab FIRST so canvas has real dimensions when drawn
    const figSection = document.getElementById("sweep-figure-section");
    if (figSection) figSection.style.display = "";
    const sweepSection = document.getElementById("sweep-section");
    if (sweepSection) sweepSection.style.display = "";
    window._activateSidebarTab?.("results");
    if (sweepChart) {
      sweepChart.invalidateRect();
      sweepChart.setSeries(points, {
        title:  "Diversity vs. number of representatives",
        xTitle: "Number of representatives (n)",
        yTitle: "Diversity covered (%)",
        xFmt:   v => String(Math.round(v)),
        yFmt:   v => v.toFixed(1) + "%",
      });
      sweepChart.setElbow(elbowX);
    }
    const caption = document.getElementById("sweep-caption");
    if (caption) caption.textContent = elbowX != null ? `Elbow at n = ${elbowX}` : "";

    // Auto-select elbow n and re-render tree
    if (elbowX != null) {
      const nInput = document.getElementById("n-input");
      if (nInput) nInput.value = elbowX;
      await sweepTo(elbowX);
    }
    hideOverlay();
  } catch (err) {
    if (err.name !== "AbortError") showOverlay("error", "⚠️", err.message);
  } finally {
    currentAbort = null;
    runBtn.removeAttribute("aria-busy");
    runBtn.textContent = "Run Analysis";
  }
}

// ── Threshold sweep (run via Run Analysis when sweep-thr-cb is checked) ───

async function runThresholdSweep() {
  if (!treeFile) {
    showOverlay("error", "⚠️", "Upload a tree file first.");
    return;
  }
  if (!alnFile) {
    showOverlay("error", "⚠️", "Threshold sweep requires an alignment FASTA file.");
    return;
  }

  const thrMin = parseFloat(document.getElementById("thr-min")?.value);
  const thrMax = parseFloat(document.getElementById("thr-max")?.value);
  if (isNaN(thrMin) || isNaN(thrMax) || thrMin >= thrMax) {
    showOverlay("error", "⚠️", "From % must be less than To % for the threshold sweep.");
    return;
  }
  if (thrMin <= 0 || thrMax >= 100) {
    showOverlay("error", "⚠️", "Threshold values must be between 0 and 100 (exclusive).");
    return;
  }

  showOverlay("loading", "⏳", "Running threshold sweep…");
  runBtn.setAttribute("aria-busy", "true");
  runBtn.textContent = "Cancel";
  currentAbort = new AbortController();

  try {
    const fd = new FormData();
    fd.append("tree",            treeFile);
    fd.append("sweep",           "threshold");
    fd.append("threshold_min",   thrMin);
    fd.append("threshold_max",   thrMax);
    fd.append("threshold_steps", "20");
    fd.append("aln_type", document.querySelector("input[name='aln-type']:checked")?.value || "nt");
    fd.append("alignment", alnFile);
    fd.append("prior",           document.getElementById("prior-input")?.value.trim() || "");
    fd.append("exclude_rep",     document.getElementById("exc-rep")?.value.trim()     || "");
    fd.append("exclude_obj",     document.getElementById("exc-obj")?.value.trim()     || "");
    fd.append("exclude_fully",   document.getElementById("exc-full")?.value.trim()    || "");
    fd.append("constrain_fully", document.getElementById("constrain")?.value.trim()   || "");
    if (weightsFile) fd.append("weights", weightsFile);

    const resp = await fetch("/api/run", { method: "POST", body: fd, signal: currentAbort.signal });
    const data = await resp.json();
    if (!resp.ok) throw new Error(data.error || "Server error");
    if (data.n_taxa != null) lastNTaxa = data.n_taxa;

    const curve = data.threshold_curve || [];
    if (curve.length === 0) throw new Error("No threshold curve data returned.");
    // Detect elbow on natural orientation (threshold→reps), then display swapped (reps→threshold)
    const elbowThreshold = findElbow(curve.map(pt => ({ x: pt.threshold, y: pt.reps })));
    const points = curve.map(pt => ({ x: pt.reps, y: pt.threshold }));

    // Show figure — activate tab FIRST so canvas has real dimensions when drawn
    const figSection = document.getElementById("sweep-figure-section");
    if (figSection) figSection.style.display = "";
    window._activateSidebarTab?.("results");
    if (sweepChart) {
      sweepChart.invalidateRect();
      sweepChart.setSeries(points, {
        title:  "Similarity threshold vs. representatives needed",
        xTitle: "Representatives to cover all taxa",
        yTitle: "Similarity threshold (%)",
        xFmt:   v => String(Math.round(v)),
        yFmt:   v => v.toFixed(1) + "%",
      });
      sweepChart.setElbow(elbowThreshold, { axis: "y" });
    }
    const caption = document.getElementById("sweep-caption");
    if (caption) caption.textContent = elbowThreshold != null ? `Elbow at threshold = ${elbowThreshold.toFixed(1)}%` : "";

    // Auto-select: set threshold, uncheck sweep, check cover, run cover
    if (elbowThreshold != null) {
      const thrInput = document.getElementById("threshold-input");
      if (thrInput) thrInput.value = elbowThreshold.toFixed(1);
      const sweepThrCb = document.getElementById("sweep-thr-cb");
      if (sweepThrCb) { sweepThrCb.checked = false; wireDependentFieldsUpdate(); }
      const coverCb = document.getElementById("cover-cb");
      if (coverCb) coverCb.checked = true;
      // Re-trigger cover run at elbow threshold (will re-enter wireRunButton handler)
      // Defer so this sweep's finally block cleans up first
      setTimeout(() => runBtn.click(), 0);
    }

    hideOverlay();
  } catch (err) {
    if (err.name !== "AbortError") showOverlay("error", "⚠️", err.message);
  } finally {
    currentAbort = null;
    runBtn.removeAttribute("aria-busy");
    runBtn.textContent = "Run Analysis";
  }
}

// ── Sweep figure export ────────────────────────────────────────────────────

function wireSweepExport() {
  const fmtSel   = document.getElementById("sweep-export-fmt");
  const themeSel = document.getElementById("sweep-export-theme");

  function updateThemeVisibility() {
    if (themeSel) themeSel.style.display = fmtSel?.value === "csv" ? "none" : "";
  }
  fmtSel?.addEventListener("change", updateThemeVisibility);
  updateThemeVisibility();

  document.getElementById("sweep-export-btn")?.addEventListener("click", () => {
    if (!sweepChart) return;
    const fmt      = fmtSel?.value || "png";
    const themeVal = themeSel?.value || "light";
    sweepChart.exportFigure(fmt, { dark: themeVal !== "light" });
  });
}

// ── Export ───────────────────────────────────────────────────────────────────

let exportFormat = "svg";
let exportTheme  = "dark";

function exportEvaluateImage() {
  const cvPrior = document.getElementById("tree-canvas-prior");
  const cvBest  = document.getElementById("tree-canvas-best");
  if (!cvPrior || !cvBest) return;

  const dark   = document.documentElement.getAttribute("data-theme") !== "light";
  const bgColor = dark ? "#0d1117" : "#f5f7fa";
  const textColor = dark ? "#e6edf3" : "#24292f";
  const GAP     = 12;
  const LABEL_H = 22;
  const FONT    = "bold 11px 'Outfit', sans-serif";

  const w = cvPrior.width + GAP + cvBest.width;
  const h = LABEL_H + Math.max(cvPrior.height, cvBest.height);

  const off = document.createElement("canvas");
  off.width = w; off.height = h;
  const ctx = off.getContext("2d");

  ctx.fillStyle = bgColor;
  ctx.fillRect(0, 0, w, h);

  ctx.fillStyle = textColor;
  ctx.font = FONT;
  ctx.textAlign = "center";
  ctx.fillText("Prior Representatives", cvPrior.width / 2, 14);
  ctx.fillText("Best Coverage", cvPrior.width + GAP + cvBest.width / 2, 14);

  ctx.drawImage(cvPrior, 0, LABEL_H);
  ctx.drawImage(cvBest,  cvPrior.width + GAP, LABEL_H);

  off.toBlob(blob => {
    if (blob) downloadBlob(blob, "parnas-evaluate.png");
  }, "image/png");
}

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
  exportDialog.show();
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

  sweepCache.clear(); // session-only cache; does not persist across loads
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
  exportBtn.style.display = "flex";
  if (layoutToggle) layoutToggle.style.display = "flex";
  document.getElementById("reset-view-btn").style.display = "flex";

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
  const sweepNRow      = document.getElementById("sweep-n-row");
  const alnSection     = document.getElementById("alignment-section");

  function update() {
    const sweepN   = document.getElementById("sweep-n-cb")?.checked;
    const sweepThr = document.getElementById("sweep-thr-cb")?.checked;

    const radiusOk = parseFloat(radiusInput.value) > 0;
    const thrVal   = parseFloat(thresholdInput.value);
    const coverOk  = radiusOk || thrVal > 0 || sweepThr;
    binaryCb.disabled = !radiusOk;
    binaryRow.classList.toggle("dimmed", !radiusOk);
    if (!radiusOk) binaryCb.checked = false;
    coverCb.disabled = !coverOk;
    coverRow.classList.toggle("dimmed", !coverOk);
    if (!coverOk) coverCb.checked = false;

    // Grey out n when covered OR n-sweep active
    const nLocked = coverCb.checked || sweepN;
    nInput.disabled    = nLocked;
    nRow.style.opacity = nLocked ? "0.45" : "1";
    if (sweepNRow) sweepNRow.style.opacity = coverCb.checked ? "0.45" : "1";

    // Grey out threshold input when threshold-sweep active
    if (thresholdInput) {
      thresholdInput.disabled    = !!sweepThr;
      thresholdInput.style.opacity = sweepThr ? "0.45" : "1";
    }

    // Show alignment section if threshold has a value OR threshold-sweep is active
    alnSection.style.display = (thrVal > 0 || sweepThr) ? "block" : "none";
  }

  // Expose update so wireSweepToggles can call it
  wireDependentFieldsUpdate = update;

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
    if (renderer && treeRoot) renderer.setAnnotations(annotations);
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
