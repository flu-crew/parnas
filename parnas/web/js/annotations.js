/**
 * Immutable annotation store.
 * All mutations return a new state object.
 *
 * State shape:
 * {
 *   branchColors:  Map<nodeId, hex string>
 *   branchWidths:  Map<nodeId, number>
 *   nodeColors:    Map<nodeId, hex string>   // PARNAS cluster colours
 *   labelOverrides: Map<nodeId, string>      // rename leaf labels
 *   metadataStrips: [ { name, data: Map<leafName, {value, color}>, type } ]
 *   cladeGroups:   [ { nodeId, label, color, collapsed } ]
 *   shapes:        [ { type, x, y, w, h, color, text, id } ]
 *   legend:        [ { label, color } ]
 * }
 */

export function emptyAnnotations() {
  return {
    branchColors:   new Map(),
    branchWidths:   new Map(),
    nodeColors:     new Map(),
    labelOverrides: new Map(),
    metadataStrips: [],
    cladeGroups:    [],
    shapes:         [],
    legend:         [],
  };
}

export function setBranchColor(ann, nodeId, color) {
  const m = new Map(ann.branchColors);
  if (color === null) m.delete(nodeId);
  else m.set(nodeId, color);
  return { ...ann, branchColors: m };
}

export function setBranchWidth(ann, nodeId, width) {
  const m = new Map(ann.branchWidths);
  m.set(nodeId, width);
  return { ...ann, branchWidths: m };
}

export function setNodeColor(ann, nodeId, color) {
  const m = new Map(ann.nodeColors);
  if (color === null) m.delete(nodeId);
  else m.set(nodeId, color);
  return { ...ann, nodeColors: m };
}

/** Set PARNAS cluster colors in bulk (from /api/run response). */
export function setClusterColors(ann, colorsObj) {
  const m = new Map(ann.nodeColors);
  for (const [name, color] of Object.entries(colorsObj)) {
    m.set(name, color); // keyed by taxon name (leaves only)
  }
  return { ...ann, nodeColors: m };
}

export function setLabelOverride(ann, nodeId, label) {
  const m = new Map(ann.labelOverrides);
  if (label === null) m.delete(nodeId);
  else m.set(nodeId, label);
  return { ...ann, labelOverrides: m };
}

/** Add or replace a metadata strip. */
export function upsertMetadataStrip(ann, strip) {
  const strips = ann.metadataStrips.filter(s => s.name !== strip.name);
  return { ...ann, metadataStrips: [...strips, strip] };
}

export function removeMetadataStrip(ann, name) {
  return { ...ann, metadataStrips: ann.metadataStrips.filter(s => s.name !== name) };
}

export function upsertCladeGroup(ann, group) {
  const groups = ann.cladeGroups.filter(g => g.nodeId !== group.nodeId);
  return { ...ann, cladeGroups: [...groups, group] };
}

export function removeCladeGroup(ann, nodeId) {
  return { ...ann, cladeGroups: ann.cladeGroups.filter(g => g.nodeId !== nodeId) };
}

export function addShape(ann, shape) {
  const id = Date.now() + Math.random();
  return { ...ann, shapes: [...ann.shapes, { ...shape, id }] };
}

export function removeShape(ann, id) {
  return { ...ann, shapes: ann.shapes.filter(s => s.id !== id) };
}

export function updateShape(ann, id, patch) {
  return { ...ann, shapes: ann.shapes.map(s => s.id === id ? { ...s, ...patch } : s) };
}

export function setLegend(ann, legend) {
  return { ...ann, legend };
}

export function setLegendCaption(ann, legendCaption) {
  return { ...ann, legendCaption };
}

/**
 * Serialise annotations to a plain JSON-safe object.
 * Maps become arrays of [k, v] pairs.
 */
export function serialise(ann) {
  return {
    branchColors:   [...ann.branchColors],
    branchWidths:   [...ann.branchWidths],
    nodeColors:     [...ann.nodeColors],
    labelOverrides: [...ann.labelOverrides],
    metadataStrips: ann.metadataStrips.map(s => ({
      ...s,
      data: s.data instanceof Map ? [...s.data] : s.data,
    })),
    cladeGroups:  ann.cladeGroups,
    shapes:       ann.shapes,
    legend:       ann.legend,
  };
}

/** Deserialise from JSON-safe object back to annotation state. */
export function deserialise(obj) {
  return {
    branchColors:   new Map(obj.branchColors   || []),
    branchWidths:   new Map(obj.branchWidths   || []),
    nodeColors:     new Map(obj.nodeColors     || []),
    labelOverrides: new Map(obj.labelOverrides || []),
    metadataStrips: (obj.metadataStrips || []).map(s => ({
      ...s,
      data: s.data instanceof Array ? new Map(s.data) : s.data,
    })),
    cladeGroups: obj.cladeGroups || [],
    shapes:      obj.shapes      || [],
    legend:      obj.legend      || [],
  };
}

/**
 * Parse a CSV metadata file into a strip entry.
 * CSV format: taxon,value[,color]
 * If no color column, auto-assign from a palette.
 */
export function parseMetadataCsv(name, csvText, sep = ",") {
  const lines = csvText.trim().split(/\r?\n/);
  const header = lines[0].toLowerCase().split(sep).map(s => s.trim());
  const hasColor = header.includes("color");
  const palette = [
    "#4e9af1","#f1724e","#4ef192","#f1d34e","#c84ef1",
    "#4ef1e8","#f14ec8","#a8f14e","#f1a84e","#4e61f1",
  ];
  const valueSet = new Set();
  const rows = [];

  for (let i = 1; i < lines.length; i++) {
    const parts = lines[i].split(sep).map(s => s.trim());
    if (parts.length < 2) continue;
    const taxon = parts[0];
    const value = parts[1];
    const color = hasColor ? parts[2] : null;
    valueSet.add(value);
    rows.push({ taxon, value, color });
  }

  // Assign palette colours to unique values
  const valueColors = new Map();
  [...valueSet].forEach((v, i) => {
    valueColors.set(v, palette[i % palette.length]);
  });

  const data = new Map();
  for (const { taxon, value, color } of rows) {
    data.set(taxon, { value, color: color || valueColors.get(value) });
  }

  const legend = [...valueColors.entries()].map(([label, color]) => ({ label, color }));
  return { name, data, type: "strip", legend };
}
