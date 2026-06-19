/**
 * Immutable tree model helpers.
 * All functions take a root node (from newick.js) and return new data.
 */

/** DFS preorder traversal, callback(node). */
export function preorder(node, fn) {
  fn(node);
  for (const c of node.children) preorder(c, fn);
}

/** DFS postorder traversal, callback(node). */
export function postorder(node, fn) {
  for (const c of node.children) postorder(c, fn);
  fn(node);
}

/** WeakMap memo for leaves(); auto-GC when node is discarded. */
const _leavesCache = new WeakMap();

/** Return all leaf nodes (array). Result is cached per root node. */
export function leaves(root) {
  const hit = _leavesCache.get(root);
  if (hit) return hit;
  const out = [];
  preorder(root, n => { if (n.children.length === 0) out.push(n); });
  _leavesCache.set(root, out);
  return out;
}

/** Return all nodes (array, preorder). */
export function allNodes(root) {
  const out = [];
  preorder(root, n => out.push(n));
  return out;
}

/**
 * Compute max depth (sum of branch lengths) from root to each node.
 * Returns Map<id, depth>.
 */
export function depthMap(root) {
  const map = new Map();
  function dfs(node, d) {
    map.set(node.id, d);
    for (const c of node.children) dfs(c, d + (c.length ?? 1));
  }
  dfs(root, 0);
  return map;
}

/**
 * Return max total depth (root-to-leaf distance) in the tree.
 * Reuses the depthMap to avoid a second traversal.
 */
export function maxDepth(root) {
  const dm = depthMap(root);
  let max = 0;
  for (const l of leaves(root)) max = Math.max(max, dm.get(l.id) || 0);
  return max;
}
