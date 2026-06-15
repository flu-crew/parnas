/**
 * Export functions: PNG, SVG, JPG, EPS, NEXUS.
 */

import { toNexus } from "./newick.js";

/**
 * Build Map<leafName, annotation object> from annotation state
 * for embedding into NEXUS output.
 */
function buildNexusAnnotations(annotations) {
  const map = new Map();
  const nc  = annotations.nodeColors  || new Map();
  const bc  = annotations.branchColors || new Map();

  // We key by leaf name for nodeColors; branchColors key by nodeId.
  // NEXUS metacomments go on nodes by name.
  for (const [name, color] of nc) {
    if (!map.has(name)) map.set(name, {});
    map.get(name)["!color"] = color;
  }

  // Metadata strips: embed value per leaf
  for (const strip of (annotations.metadataStrips || [])) {
    for (const [taxon, entry] of (strip.data instanceof Map ? strip.data : new Map())) {
      if (!map.has(taxon)) map.set(taxon, {});
      map.get(taxon)[strip.name] = entry.value;
    }
  }

  return map;
}

/**
 * Export annotated NEXUS file embedding visual annotations.
 * Includes PearTree/FigTree-compatible [&!color=] metacomments.
 *
 * @param {object} root        - parsed tree root
 * @param {object} annotations - annotation state
 * @param {object} appState    - {representatives, mode, diversity, ...} from /api/run
 * @returns {string} NEXUS text
 */
export function buildAnnotatedNexus(root, annotations, appState = {}) {
  const extraAnnotations = buildNexusAnnotations(annotations);

  // Also embed PARNAS-specific data (cluster membership)
  if (appState.representatives) {
    const repSet = new Set(appState.representatives);
    for (const rep of repSet) {
      if (!extraAnnotations.has(rep)) extraAnnotations.set(rep, {});
      extraAnnotations.get(rep)["parnas_rep"] = "true";
    }
  }

  let nexus = toNexus(root, extraAnnotations);

  // Append PARNAS results block as a comment
  const meta = {
    mode:            appState.mode,
    representatives: appState.representatives,
    prior_centers:   appState.prior_centers,
    diversity:       appState.diversity,
    generated:       new Date().toISOString(),
  };
  nexus += `\n[PARNAS results: ${JSON.stringify(meta)}]\n`;

  return nexus;
}

/**
 * Export session (tree + annotations) as JSON sidecar.
 */
export function buildSessionJson(newickStr, annotations, appState = {}) {
  // Serialise Maps → arrays for JSON
  function serMap(m) {
    return m instanceof Map ? [...m] : m;
  }
  return JSON.stringify({
    version:     1,
    newick:      newickStr,
    appState:    appState,
    annotations: {
      branchColors:   serMap(annotations.branchColors),
      branchWidths:   serMap(annotations.branchWidths),
      nodeColors:     serMap(annotations.nodeColors),
      labelOverrides: serMap(annotations.labelOverrides),
      metadataStrips: (annotations.metadataStrips || []).map(s => ({
        ...s, data: serMap(s.data),
      })),
      cladeGroups:  annotations.cladeGroups,
      shapes:       annotations.shapes,
      legend:       annotations.legend,
    },
  }, null, 2);
}

/** Trigger browser download of a Blob. */
export function downloadBlob(blob, filename) {
  const url = URL.createObjectURL(blob);
  const a   = document.createElement("a");
  a.href     = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

/** Download a string as a text file. */
export function downloadText(text, filename, mime = "text/plain") {
  downloadBlob(new Blob([text], { type: mime }), filename);
}

/**
 * Convert SVG string to PNG/JPG blob via off-screen canvas.
 * @param {string} svgStr
 * @param {number} width
 * @param {number} height
 * @param {string} bgColor
 * @param {"png"|"jpg"} fmt
 * @returns {Promise<Blob>}
 */
export async function svgToRasterBlob(svgStr, width, height, bgColor, fmt = "png") {
  const blob = new Blob([svgStr], { type: "image/svg+xml;charset=utf-8" });
  const url  = URL.createObjectURL(blob);

  const canvas = document.createElement("canvas");
  canvas.width  = width;
  canvas.height = height;
  const ctx = canvas.getContext("2d");
  ctx.fillStyle = bgColor;
  ctx.fillRect(0, 0, width, height);

  await new Promise((resolve, reject) => {
    const img = new Image();
    img.onload  = () => { ctx.drawImage(img, 0, 0, width, height); URL.revokeObjectURL(url); resolve(); };
    img.onerror = e  => { URL.revokeObjectURL(url); reject(e); };
    img.src = url;
  });

  const mime = fmt === "jpg" ? "image/jpeg" : "image/png";
  return new Promise(resolve => canvas.toBlob(resolve, mime, 0.95));
}

/**
 * Build EPS from a canvas element.
 * Encodes the raster image as JPEG hex data inside PS.
 *
 * @param {HTMLCanvasElement} canvas
 * @param {number} logicalW  - logical px width (for pt sizing at 96dpi baseline)
 * @param {number} logicalH
 * @returns {string} EPS text
 */
export function canvasToEps(canvas, logicalW, logicalH) {
  const ptW  = Math.round((logicalW / 96) * 72);
  const ptH  = Math.round((logicalH / 96) * 72);
  const imgW = canvas.width;
  const imgH = canvas.height;

  const b64    = canvas.toDataURL("image/jpeg", 0.92).split(",")[1];
  const binary = atob(b64);
  let hex = "";
  for (let i = 0; i < binary.length; i++) {
    hex += binary.charCodeAt(i).toString(16).padStart(2, "0");
    if ((i + 1) % 40 === 0) hex += "\n";
  }

  return `%!PS-Adobe-3.0 EPSF-3.0
%%BoundingBox: 0 0 ${ptW} ${ptH}
%%EndComments
gsave
${ptW} ${ptH} scale
${imgW} ${imgH} 8 [${imgW} 0 0 -${imgH} 0 ${imgH}]
currentfile /ASCIIHexDecode filter /DCTDecode filter
false 3 colorimage
${hex}
>
grestore
%%EOF`;
}
