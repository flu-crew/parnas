# -*- coding: utf-8 -*-
from collections import deque
from typing import Dict, List, Optional, Set

import numpy as np
from dendropy import Tree, Node
from datetime import datetime

from numba import int32, float64, bool_, typed, types
from numba.experimental import jitclass

from .medoid_utils import DistFunction
from .tree_indexer import TreeIndexer
from parnas.logging import parnas_logger
from .fast_pmedian_finder import FastPMedianFinder, PMedianDP


class TreeCoverage(FastPMedianFinder):
    """
    Finding a set of leaves that cover the entire tree (under a fixed radius).

    This implementation uses numba's JIT to achieve better runtimes.
    """

    def __init__(self, tree: Tree):
        super().__init__(tree)

    def find_coverage(self, radius: float, cost_map: Dict[str, float], covered_leaves: Set[str]) -> Optional[List[str]]:
        self.cost_map = cost_map
        self.distance_functions = None
        self.prior_coverage = np.full(self.nleaves, False, dtype=np.bool_)
        for leaf in self.tree.leaf_nodes():
            if leaf.taxon.label in covered_leaves:
                self.prior_coverage[leaf.index] = True
        parnas_logger.debug(f"Prior coverage: {self.prior_coverage}")
        self._initialize_lookups()

        postorder = np.array([node.index for node in self.tree.postorder_node_iter()], dtype=np.int32)
        coverage_jit = CoverageDP(self.children, postorder, self.is_ancestor, self.leaf_lists, self.subtree_leaves,
                                  self.distance_lookup, self.index_lookup, self.nleaves, self.prior_coverage,
                                  self.costs, self.parents, self.edge_lens, self.tree.seed_node.index, radius)
        self.G, self.F = coverage_jit.run_dp()
        # parnas_logger.debug(str(self.G.tolist()))
        # parnas_logger.debug(str(self.F.tolist()))
        parnas_logger.debug(f'Finished DP {datetime.now().strftime("%H:%M:%S")}')

        root_id = self.tree.seed_node.index
        if self.G[root_id][-1] == np.inf:
            parnas_logger.debug('Impossible to cover the full tree given the constraints.')
            return None
        else:
            parnas_logger.debug(f"Optimal # of reps: {self.G[root_id][-1]}")
            # Do backtracking.
            coverage = self.coverage_backtrack(root_id)
            return coverage

    def coverage_backtrack(self, root_id: int):
        coverage_ids = []
        coverage = []
        self._backtrack_links_coverage(True, root_id, self.leaf_lists[root_id, -1], coverage_ids)
        median_ids = set(coverage_ids)
        for leaf in self.tree.leaf_nodes():
            if leaf.index in median_ids:
                coverage.append(leaf.taxon.label)
        return coverage

    def _backtrack_links_coverage(self, in_G: bool, node_id: int, radius_id: int, coverage: List[int]):
        backtrack_queue = deque()
        backtrack_queue.append((in_G, node_id, radius_id))
        while backtrack_queue:
            in_G, node_id, radius_id = backtrack_queue.popleft()
            if self.children[node_id, 0] < 0:
                # We are at a leaf.
                if (in_G and self.G[node_id, self.index_lookup[node_id, radius_id]] > 0) or\
                        (not in_G and self.F[node_id, self.index_lookup[node_id, radius_id]] > 0):
                    coverage.append(node_id)
                # if in_G or radius_id == node_id:
            else:
                # Internal node.
                if in_G:
                    cur_val = self.G[node_id, self.index_lookup[node_id, radius_id]]
                    if self.index_lookup[node_id, radius_id] > 0 and\
                            cur_val == self.G[node_id, self.index_lookup[node_id, radius_id] - 1]:
                        next_rad_id = self.leaf_lists[node_id, self.index_lookup[node_id, radius_id] - 1]
                        backtrack_queue.append((True, node_id, next_rad_id))
                    else:
                        left, right = self.children[node_id]
                        c1, c2 = (left, right) if (self.is_ancestor[left, radius_id]) else (right, left)
                        backtrack_queue.append((True, c1, radius_id))
                        backtrack_queue.append((False, c2, radius_id))
                else:
                    # In F:
                    cur_val = self.F[node_id, self.index_lookup[node_id, radius_id]]
                    if cur_val == self.G[node_id][self.index_lookup[node_id, radius_id]]:
                        backtrack_queue.append((True, node_id, radius_id))
                    else:
                        left, right = self.children[node_id]
                        backtrack_queue.append((False, left, radius_id))
                        backtrack_queue.append((False, right, radius_id))


@jitclass([('children_arr', int32[:, ::1]), ('postorder_arr', int32[::1]), ('is_ancestor_arr', bool_[:, ::1]),
           ('leaf_lists', int32[:, ::1]), ('subtree_leaves', int32[::1]), ('distance_lookup', float64[:, ::1]),
           ('index_lookup', int32[:, ::1]), ('nleaves', int32), ('nnodes', int32), ('prior_coverage', bool_[::1]),
           ('costs', float64[::1]), ('parents', int32[::1]), ('edge_lengths', float64[::1]), ('root_id', int32),
           ('radius', float64), ('G', float64[:, ::1]), ('F', float64[:, ::1])
           ])
class CoverageDP:
    """
    Numba-optimized dynamic programming solution for the tree coverage problem.
    """

    def __init__(self, children_arr, postorder_arr, is_ancestor_arr, leaf_lists, subtree_leaves, distance_lookup,
                 index_lookup, nleaves: int, prior_coverage, costs, parents, edge_lengths, root_id: int, radius: float):
        self.children_arr = children_arr
        self.postorder_arr = postorder_arr
        self.is_ancestor_arr = is_ancestor_arr
        self.leaf_lists = leaf_lists
        self.subtree_leaves = subtree_leaves
        self.distance_lookup = distance_lookup
        self.index_lookup = index_lookup
        self.prior_coverage = prior_coverage
        self.costs = costs
        self.nleaves = nleaves
        self.nnodes = len(postorder_arr)
        self.parents = parents
        self.edge_lengths = edge_lengths
        self.root_id = root_id
        self.radius = radius

        self._initialize_lookups()
        self.G = np.empty((self.nnodes, self.nleaves), dtype=np.float64)
        self.F = np.empty((self.nnodes, self.nleaves), dtype=np.float64)

    def run_dp(self):
        for node_id in self.postorder_arr:
            # Note that a cost of infinity indicates that a node cannot be selected.
            if self._is_leaf(node_id):
                self.G[node_id].fill(1 if self.costs[node_id] < np.inf else np.inf)
                if self.prior_coverage[node_id]:
                    # The leaf is already covered by prior centers.
                    self.F[node_id].fill(0)
                else:
                    for ind, rad_id in enumerate(self.leaf_lists[node_id]):
                        self_coverage = 1 if self.costs[node_id] < np.inf else np.inf
                        if self.costs[rad_id] < np.inf:
                            self.F[node_id, ind] = 0 if (self.distance_lookup[node_id, rad_id] <= self.radius) else self_coverage
                        else:
                            self.F[node_id, ind] = np.inf
            else:
                for ind, rad_id in enumerate(self.leaf_lists[node_id]):
                    if self.is_ancestor_arr[node_id, rad_id]:
                        prior_G = np.inf if ind == 0 else self.G[node_id, ind - 1]
                        if self.is_ancestor_arr[self.children_arr[node_id, 0], rad_id]:
                            c1, c2 = self.children_arr[node_id, 0], self.children_arr[node_id, 1]
                        else:
                            c1, c2 = self.children_arr[node_id, 1], self.children_arr[node_id, 0]
                        self.G[node_id, ind] = min(
                            prior_G,
                            self._get_DP_by_id(c1, rad_id, True) + self._get_DP_by_id(c2, rad_id, False)
                        )
                        self.F[node_id, ind] = self.G[node_id, ind]
                    else:
                        left, right = self.children_arr[node_id]
                        self.G[node_id, ind] = np.inf if ind == 0 else self.G[node_id, ind - 1]
                        if self.costs[rad_id] < np.inf:
                            self.F[node_id, ind] = min(
                                self.G[node_id, ind],
                                self._get_DP_by_id(left, rad_id, False) + self._get_DP_by_id(right, rad_id, False)
                            )
                        else:
                            self.F[node_id, ind] = np.inf
        return self.G, self.F

    def _is_leaf(self, node_id) -> bool:
        return self.children_arr[node_id, 0] < 0

    def _get_DP_by_id(self, node_id, radius_id, in_G):
        radius_list_index = self.index_lookup[node_id, radius_id]
        return self._get_DP_by_list_index(node_id, radius_list_index, in_G)

    def _get_DP_by_list_index(self, node_id, radius_list_index, in_G):
        source = self.G if in_G else self.F
        if radius_list_index < len(source[node_id]):
            return source[node_id][radius_list_index]
        else:
            return source[node_id][-1]

    # Code below is a copy-paste from PMedianDP. jitclass does not support subtyping.
    # TODO: refactor by taking out the methods below from jitclass structure.
    # ----------------------

    def _initialize_lookups(self):
        """
        Constructs sorted node lists for each node.
        Initializes leaf_lists, index_lookup map + other helper structures.
        TODO: re-implement to assert O(n^2) runtime. Currently O(n^2 logn).
        """
        for node1_id in range(self.nnodes):
            node_dist_pairs = typed.List()
            # Compute distances from node1 to all other nodes:
            for node2_id, dist in zip(*self._dfs_tree_traversal(node1_id)):
                if self._is_leaf(node2_id):
                    self.distance_lookup[node1_id, node2_id] = dist

            # for node2 in self.tree.find_clades(lambda x: node1.is_parent_of(x), order='postorder'):
            subtree_size = 0
            for node2_id in self._postorder_iter(node1_id):
                self.is_ancestor_arr[node1_id, node2_id] = True
                if self._is_leaf(node2_id):
                    node_dist_pairs.append((node2_id, self.distance_lookup[node1_id, node2_id]))
                    subtree_size += 1
            self.subtree_leaves[node1_id] = subtree_size
            # for node2 in self.tree.find_clades(lambda x: not node1.is_parent_of(x), order='postorder'):
            for node2_id in self._postorder_iter(self.root_id, node1_id):
                if self._is_leaf(node2_id):
                    node_dist_pairs.append((node2_id, self.distance_lookup[node1_id, node2_id]))
            node_dist_pairs = sorted(node_dist_pairs, key=lambda r: r[1])
            # node_dist_pairs.sort(key=lambda r: r[1])  # sort pairs by distance.

            self.leaf_lists[node1_id, :] = np.array([node2_id for node2_id, dist in node_dist_pairs], dtype=np.int32)

            # Mapping a node to its index in the sorted list for node1:
            for i, (node2_id, dist) in enumerate(node_dist_pairs):
                self.index_lookup[node1_id, node2_id] = i

    def _dfs_tree_traversal(self, start_node_id: int):
        # DFS stacks (node with prev_node) + distance to node:
        node_stack = np.full((self.nnodes, 2), -1, np.int32)
        dist_stack = np.zeros(self.nnodes, np.float64)  # mirrors node_stack.
        top = 0

        nodes = np.empty(self.nnodes, np.int32)  # list of nodes in dfs order.
        dists = np.empty(self.nnodes, np.float64)  # list of distances to nodes in dfs order.
        i = 0
        node_stack[top, 0], node_stack[top, 1] = start_node_id, -1
        dist_stack[top] = 0
        top = 1
        while top > 0:
            node_id, prev_node_id = node_stack[top - 1]
            cur_dist = dist_stack[top - 1]
            top -= 1
            nodes[i], dists[i] = node_id, cur_dist
            i += 1

            neighbors = typed.List()
            if self.children_arr[node_id, 0] >= 0:
                neighbors.append(self.children_arr[node_id, 0])
                neighbors.append(self.children_arr[node_id, 1])
            if self.parents[node_id] >= 0:
                neighbors.append(self.parents[node_id])
            for neighbor_id in neighbors:
                if neighbor_id != prev_node_id:
                    is_parent = self.parents[node_id] == neighbor_id
                    edge_len = self.edge_lengths[node_id] if is_parent else self.edge_lengths[neighbor_id]
                    node_stack[top, 0], node_stack[top, 1] = neighbor_id, node_id
                    dist_stack[top] = cur_dist + edge_len
                    top += 1
        return nodes, dists

    def _postorder_iter(self, start_node_id: int, filter_node_id=None):
        node_stack = np.empty((self.nnodes, 2), dtype=np.int32)
        node_stack[0, 0], node_stack[0, 1] = start_node_id, 0
        top = 1
        while top > 0:
            node_id, state = node_stack[top - 1]
            top -= 1
            if state:
                yield node_id
            elif node_id != filter_node_id:
                node_stack[top, 0], node_stack[top, 1] = node_id, 1
                top += 1
                if not self._is_leaf(node_id):
                    for child_id in [self.children_arr[node_id, 1], self.children_arr[node_id, 0]]:
                        node_stack[top, 0], node_stack[top, 1] = child_id, 0
                        top += 1
