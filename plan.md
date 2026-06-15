# PARNAS UI — paper.js Annotation Layer Plan

## Status: Phase 1 complete (unverified in browser)

---

## What was built

- Removed phylotree.js (D3/SVG); added paper.js CDN
- Flask now serves `/<path:filename>` for all static files under `web/`
- 8 ES modules in `parnas/web/js/`:
  - `newick.js` — Newick/NEXUS parser + serialiser + `toNexus()`
  - `treemodel.js` — traversal, MRCA, depth, parent map
  - `layout.js` — Reingold-Tilford rectangular + radial polar layouts
  - `renderer.js` — paper.js canvas renderer, pan/zoom, hit-testing, clade collapse triangles, metadata strips, legend, user shapes
  - `annotations.js` — immutable annotation store (branchColors, branchWidths, nodeColors, labelOverrides, metadataStrips, cladeGroups, shapes, legend)
  - `export.js` — SVG, PNG, JPG, EPS, annotated NEXUS, session JSON
  - `panel.js` — Visual Options drawer (PearTree-inspired): node colour, label override, clade group, branch colour/width, display opts, CSV metadata import, legend, shapes
  - `app.js` — orchestrator: form → /api/run → parse → layout → render → results panel → export dialog → theme

- `index.html` changes:
  - paper.js CDN replaces phylotree scripts
  - `<canvas id="tree-canvas">` mounted dynamically by app.js
  - "Annotate" toolbar button + `#vopt-panel` drawer inside `#tree-card`
  - Export dialog: added NEXUS format button + "Save session" button
  - Old inline `<script>` replaced with `<script type="module" src="js/app.js">`

---

## Remaining work

### Phase 1 — Renderer verification (NEXT)
- [ ] Load real influenza Newick in browser, confirm tree draws
- [ ] Verify pan/zoom works (mouse drag + ctrl+wheel)
- [ ] Verify node/branch click hit-testing fires panel
- [ ] Verify radial layout toggle re-renders correctly
- [ ] Fix any layout/padding bugs (label clipping, strip overlap)

### Phase 2 — Export validation
- [ ] SVG export via `paper.project.exportSVG` — check fonts/styles preserved
- [ ] PNG/JPG via SVG→canvas path — check bgColor and scaling
- [ ] EPS — verify hex-encoded JPEG inside PS renders in Ghostscript/Illustrator
- [ ] NEXUS export — round-trip test (write → re-open in FigTree/PearTree, colours preserved)

### Phase 3 — Annotation completeness
- [ ] "Apply to clade" branch colour — recursive `setBranchColor` needs browser test
- [ ] Clade collapse triangles — confirm visual + label positioning
- [ ] Metadata strip `stripStartX` — fix potential clipping on narrow windows
- [ ] CSV TSV import with tab separator (panel.js passes sep but `parseMetadataCsv` ignores it — fix)
- [ ] Continuous colour scale for metadata (currently categorical only)

### Phase 4 — Session persistence
- [ ] Session JSON save — done
- [ ] Session JSON load/import — NOT implemented; add drag-drop or file picker to load `.json` sidecar

### Known bugs to fix
- `parseMetadataCsv` signature takes `(name, csvText, sep)` but sep is ignored — fix to honour tab
- `panel.js` `_emit(this._ann, { showSupport: ... })` passes optsPatch as second arg; verify `app.js` `onChange(newAnn, optsPatch)` wiring handles it correctly for display-only opts (radial, fontSize, showSupport) that don't mutate annotations
- `renderer.js` `drawSubtree` skips `drawnNodes` set that was stubbed but never populated — remove dead code

---

## Architecture decisions (locked)
- paper.js canvas (not SVG/D3) — needed for >5k tip influenza trees
- All rendering client-side; Flask only serves Newick + PARNAS metadata
- Plain ES modules served by Flask; no build step
- paper.js CDN: `https://cdnjs.cloudflare.com/ajax/libs/paper.js/0.12.17/paper-full.min.js`
- NEXUS export with `[&!color=hex]` metacomments (FigTree/PearTree-compatible)
- Metadata import: CSV/TSV (taxon, value[, color])
- Session persistence: JSON sidecar (tree + annotations)

### Phase 5 — Interactive Parameter Sweep (NEW)

Paper scenario: "Where Does the Elbow Actually Fall?" — CLI produces individual results; only the web UI shows live color-cluster restructuring as n changes.

**Architecture**: tree layout stays fixed between sweeps (same Newick, same coords). Only cluster colors change → each step = one POST /api/run + canvas repaint.

- [ ] **Phase 5a — Sweep controls** (`index.html` + `app.js`)
  - Sweep bar below results: [−] [n=6] [+] Diversity: 34.2%
  - `wireSweepBar()` + `sweepTo(n)` in app.js
  - `sweepTo`: reuses treeFile + form params, only changes `n`; cancels in-flight prior request
  - Disable sweep bar when `cover=true` (n is irrelevant in cover mode)
  - Hide sweep bar in evaluate mode (no n to sweep)

- [ ] **Phase 5b — Diversity sparkline** (`parnas/web/js/sweep.js` new module)
  - `SweepChart(canvasEl)` class — draws diversity-vs-n line chart
  - `addPoint(n, diversity)`, `reset()`, `draw()`
  - Appears after ≥2 data points; vertical marker on current n

- [ ] **Phase 5c — Snapshot capture**
  - "📷 Capture" button saves `{n, colors, diversity}` to `snapshots[]`
  - Snapshot shelf: row of 120×80px thumbnail canvases labeled "n=3", "n=6", etc.
  - Click thumbnail → restore that color state to main view
  - "Export snapshots" in export dialog → individual PNGs per snapshot
  - `destroy()` thumbnail renderers on clear (avoid paper.js scope leak)

**Risks**: rapid clicks need debounce/cancel-previous; snapshot shelf with many TreeRenderer instances needs explicit destroy on clear; null diversity in evaluate mode.

---

## Reference
- iTOL uses paper.js for same reason (canvas perf on large trees)
- PearTree (peartree.live / github.com/artic-network/peartree) — UI patterns borrowed: Visual Options drawer, clade collapse, CSV metadata, NEXUS with embedded settings
