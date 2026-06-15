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

/** Return all leaf nodes (array). */
export function leaves(root) {
  const out = [];
  preorder(root, n => { if (n.children.length === 0) out.push(n); });
  return out;
}

/** Return all nodes (array, preorder). */
export function allNodes(root) {
  const out = [];
  preorder(root, n => out.push(n));
  return out;
}

/** Build a Map<id, node> for fast lookup. */
export function indexById(root) {
  const map = new Map();
  preorder(root, n => map.set(n.id, n));
  return map;
}

/** Build a Map<name, node> for leaf/internal name lookup. */
export function indexByName(root) {
  const map = new Map();
  preorder(root, n => { if (n.name) map.set(n.name, n); });
  return map;
}

/** Build a Map<id, parentNode> for upward traversal. */
export function parentMap(root) {
  const map = new Map();
  preorder(root, n => {
    for (const c of n.children) map.set(c.id, n);
  });
  return map;
}

/**
 * Return the MRCA (most recent common ancestor) of a set of taxon names.
 * Returns null if no names found.
 */
export function mrca(root, names) {
  const nameSet = new Set(names);
  const parents = parentMap(root);
  const nameIndex = indexByName(root);

  // Get path to root for each target node
  function pathToRoot(node) {
    const path = [];
    let cur = node;
    while (cur) {
      path.push(cur.id);
      cur = parents.get(cur.id);
    }
    return path;
  }

  const targets = [...nameSet].map(n => nameIndex.get(n)).filter(Boolean);
  if (targets.length === 0) return null;
  if (targets.length === 1) return targets[0];

  // Find common ancestor by intersecting paths
  const paths = targets.map(pathToRoot);
  const sets  = paths.map(p => new Set(p));
  const common = paths[0].filter(id => sets.every(s => s.has(id)));
  if (common.length === 0) return root;

  const byId = indexById(root);
  return byId.get(common[0]) || root;
}

/**
 * Return all leaf names under a subtree rooted at `node`.
 */
export function subtreeLeafNames(node) {
  return leaves(node).map(l => l.name).filter(Boolean);
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
 */
export function maxDepth(root) {
  let max = 0;
  const dm = depthMap(root);
  for (const l of leaves(root)) max = Math.max(max, dm.get(l.id) || 0);
  return max;
}
