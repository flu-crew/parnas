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

/** Set PARNAS cluster colors in bulk (from /api/run response). */
export function setClusterColors(ann, colorsObj) {
  const m = new Map(ann.nodeColors);
  for (const [name, color] of Object.entries(colorsObj)) {
    m.set(name, color); // keyed by taxon name (leaves only)
  }
  return { ...ann, nodeColors: m };
}

export function setLegend(ann, legend) {
  return { ...ann, legend };
}

export function setLegendCaption(ann, legendCaption) {
  return { ...ann, legendCaption };
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
