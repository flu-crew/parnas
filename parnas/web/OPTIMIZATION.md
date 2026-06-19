# Parnas Web Viewer — Optimization Notes

Prioritized findings from a full audit of `parnas/web/` (JS/HTML) and `parnas/server.py`.
Implemented items are marked **✅ done**. Python/server items are documented only — not applied.

---

## Frontend — Implemented ✅

| # | File | Change | Impact |
|---|------|--------|--------|
| 1 | `renderer.js:_setupInteraction` | rAF-throttle pan mousemove: accumulate dx/dy, flush `_applyTransform()` once per frame via `requestAnimationFrame` | **med** — eliminates per-event style recalcs during drag |
| 2 | `app.js:wireRunButton` | Run button reads/writes shared `sweepCache` (LRU cap 20); sample re-runs skip server | **med-high** — instant re-run for cached params |
| 3 | `app.js:sweepCache` | LRU eviction via `_sweepCacheSet()`: Map insertion-order + delete oldest when `> 20` | **med** — bounds memory on long sweep sessions |
| 4 | `sweep.js:draw` | Cache `getBoundingClientRect()` in `this._rect`; only re-measure after `invalidateRect()` (reset/resize) | **low** — eliminates forced reflow per sweep click |
| 5 | `index.html` | `rel="preconnect"` for 3 CDN origins; `defer` on Paper.js `<script>` | **med** — unblocks HTML parse; reduces DNS/TLS latency |

---

## Frontend — Not Yet Implemented (higher risk / larger scope)

### Rendering (renderer.js)

**`_effectiveColors()` rebuilt every `_redraw()`** (`renderer.js:301`)  
Full leaf + postorder walk each render. Inside `drawSubtree`, `cladeGroups.find()` runs per-node → O(n × groups).  
Fix: hoist a `Map<nodeId, cladeGroup>` before traversal; memoize `_effectiveColors` keyed on a hash of `(ann.nodeColors, repColors)`.  
Impact: **med-high** for large trees / any clade groups.

**Linear hit-test scan on click** (`renderer.js:_hitNodes`, `_hitBranches`)  
Every click iterates all stored items checking `contains`/`hitTest`.  
Fix: use Paper's `project.hitTest(pt, {...})` (internal spatial index) or bucket by y-band.  
Impact: **med** for 1000+ leaf trees.

**`leaves()` recomputed in multiple hot paths**  
`renderer.js:_effectiveColors`, `app.js` resize/sweep callbacks, `layout.js` — each is a full preorder walk.  
Fix: cache `root._leaves` after parse; invalidate only on tree change.  
Impact: **med** (compounds on every render/resize/sweep step).

**Resize triggers full `render()` (item rebuild)**  (`app.js:413,475`)  
Resize only changes pixel scale, not topology. Full `_redraw()` clears and recreates every Paper item.  
Fix: for pure resize, call `view.viewSize` + `_applyTransform()` instead; only rebuild when width crosses label-fit thresholds.  
Impact: **med** (jank on resize of large trees).

**Renderer/canvas destroyed and recreated on every Run** (`app.js:399,433`)  
`mountRenderer` sets `treeCard.innerHTML = ""` discarding the Paper scope; a new scope is built each run.  
Fix: reuse canvas + renderer when topology is unchanged (just call `render()` again).  
Impact: **med** (GC churn, repeat-run latency).

**No `Symbol`/`SymbolDefinition` reuse for leaf dots**  
Each leaf dot is a unique `paper.Path.Circle`. With 2000 leaves = 2000 separate items.  
Fix: `paper.SymbolDefinition` for the dot; place `SymbolItem` instances.  
Impact: **med** (first-paint time, item count).

---

## Backend / Server — Documented Only (no Python changes applied)

### High impact

**No gzip on `/api/run` JSON** (`server.py:~393`)  
Newick string + colors dict for a large tree is hundreds of KB, sent uncompressed. Every sweep step pays this cost.  
Fix: `pip install flask-compress`, add `Compress(app)` in `server.py`. Trees compress ~5–10×.  
Impact: **high** — dominates sweep latency on any non-localhost connection.

**No server-side result cache** (`server.py:run_parnas`)  
Identical params re-run the full medoid solve even after page reload or from a second browser tab. Client-side `sweepCache` only covers the current session.  
Fix: `functools.lru_cache` or `cachetools.LRUCache` keyed on `(tree_hash, params_tuple)`. Tree hash = `hashlib.md5(tree_bytes).hexdigest()`.  
Impact: **high** — eliminates redundant solves across reloads and tabs.

### Medium impact

**Synchronous blocking compute on request thread** (`server.py:305, 363, 375, 317`)  
`find_n_medoids*` and the 1000-iteration random-permutation loop (evaluate mode) run inline on Flask's single default thread, blocking all other requests (including the warmup health check).  
Fix: `app.run(threaded=True)` (one line); consider vectorizing or caching `rnd_scores` per tree+params to reduce the permutation loop cost.  
Impact: **med-high** (evaluate mode latency; concurrent-user blocking).

**No cache headers for static JS/CSS** (`server.py:143 static_files`)  
`send_from_directory` uses Flask defaults (no `Cache-Control`). JS modules re-fetched/revalidated on every page load.  
Fix: `send_file(..., max_age=3600)` for `/js/*` and assets; add ETag.  
Impact: **med** (repeat-load round-trips).

**Self-host CDN assets** (`index.html:8–14`)  
Three CDN origins (jsdelivr, Google Fonts, cdnjs) require DNS + TLS on first load. For the PyInstaller-frozen offline build, CDN is unreachable.  
Fix: vendor `pico.min.css`, `paper-full.min.js`, and the font files into `parnas/web/static/`; serve them via the existing static route.  
Impact: **med** (first-load; critical for offline/frozen distribution).

### Low impact

**Exception-driven Newick/Nexus schema detection** (`server.py:214`)  
Parser tries schemas sequentially catching exceptions. Cheap sniff: first non-whitespace char `(` → Newick, `#` → Nexus.  
Impact: **low**.

**`sweepCache` key uses `hasWeights: !!weightsFile`** (`app.js:~631`) *(correctness)*  
Boolean flag — if user swaps weights file without changing other params, stale cache is served.  
Fix: include `weightsFile.name + weightsFile.size` (or a content hash) in the key.  
Impact: **low-med (correctness)**.
