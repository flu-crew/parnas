/**
 * Tree layout engine.
 * Returns Map<id, {x, y}> where:
 *   rectangular: x = horizontal (depth), y = vertical (leaf order)
 *   radial:      x = angle (radians), y = radius
 *
 * All coordinates are in [0,1] normalised space.
 * Renderer scales to canvas size.
 */

import { leaves, postorder, depthMap, maxDepth } from "./treemodel.js";

/**
 * Rectangular phylogram layout.
 * x proportional to cumulative branch length.
 * y evenly spaced by leaf order.
 *
 * @param {object} root   - parsed tree root node
 * @returns {Map<id,{x,y}>}
 */
export function rectangularLayout(root) {
  const coords = new Map();
  const leafList = leaves(root);
  const nLeaves  = leafList.length;
  const dm       = depthMap(root);
  const maxD     = maxDepth(root) || 1;

  // Assign y to leaves in order
  leafList.forEach((leaf, i) => {
    coords.set(leaf.id, { x: dm.get(leaf.id) / maxD, y: (i + 0.5) / nLeaves });
  });

  // Assign internal nodes: y = mean of children y, x = own depth
  postorder(root, node => {
    if (node.children.length > 0) {
      const ys = node.children.map(c => coords.get(c.id)?.y ?? 0);
      coords.set(node.id, {
        x: dm.get(node.id) / maxD,
        y: (Math.min(...ys) + Math.max(...ys)) / 2,
      });
    }
  });

  return coords;
}

/**
 * Rectangular cladogram layout (equal branch lengths, tips aligned).
 */
export function cladogramLayout(root) {
  const coords = new Map();
  const leafList = leaves(root);
  const nLeaves  = leafList.length;

  // Compute node depth by edge count
  const edgeDepth = new Map();
  function edgeDfs(node, d) {
    edgeDepth.set(node.id, d);
    for (const c of node.children) edgeDfs(c, d + 1);
  }
  edgeDfs(root, 0);

  // Max depth by edge count
  const maxEdge = Math.max(...leafList.map(l => edgeDepth.get(l.id) || 0)) || 1;

  leafList.forEach((leaf, i) => {
    coords.set(leaf.id, { x: 1.0, y: (i + 0.5) / nLeaves });
  });

  postorder(root, node => {
    if (node.children.length > 0) {
      const ys = node.children.map(c => coords.get(c.id)?.y ?? 0);
      coords.set(node.id, {
        x: edgeDepth.get(node.id) / maxEdge,
        y: (Math.min(...ys) + Math.max(...ys)) / 2,
      });
    }
  });

  return coords;
}

/**
 * Radial layout.
 * Returns coords where x = angle in [0, 2π), y = radius in [0, 1].
 * Suitable for paper.js polar→cartesian transform in renderer.
 */
export function radialLayout(root) {
  const leafList = leaves(root);
  const nLeaves  = leafList.length;
  const dm       = depthMap(root);
  const maxD     = maxDepth(root) || 1;

  // Assign angle to leaves evenly around circle
  const angleCoords = new Map();
  leafList.forEach((leaf, i) => {
    angleCoords.set(leaf.id, {
      x: (i / nLeaves) * 2 * Math.PI,
      y: dm.get(leaf.id) / maxD,
    });
  });

  postorder(root, node => {
    if (node.children.length > 0) {
      const angles = node.children.map(c => angleCoords.get(c.id)?.x ?? 0);
      angleCoords.set(node.id, {
        x: (Math.min(...angles) + Math.max(...angles)) / 2,
        y: dm.get(node.id) / maxD,
      });
    }
  });

  return angleCoords;
}

/**
 * Convert polar coords (angle, radius) to cartesian (x, y) in [-1,1] space.
 * cx, cy: centre in normalised canvas coords (0.5, 0.5 for centred).
 */
export function polarToCartesian(angle, radius, cx = 0, cy = 0) {
  return {
    x: cx + radius * Math.cos(angle),
    y: cy + radius * Math.sin(angle),
  };
}
