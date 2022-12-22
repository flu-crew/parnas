# -*- coding: utf-8 -*-
from collections import deque
from typing import Dict, List, Optional

import numpy as np
from dendropy import Tree, Node
from datetime import datetime

from numba import int32, float64, bool_, typed, types
from numba.experimental import jitclass

from .medoid_utils import DistFunction
from .tree_indexer import TreeIndexer
from parnas.logging import parnas_logger


class FastPMedianFinder(object):
    """
    Class that finds the median nodes (among leaves!) given a phylogenetic tree.
    An adapted implementation of Tamir's algorithm (Tamir 1996) for generalized p-medians.

    This implementation actively uses numba's JIT to achieve better runtimes.
    Additionally, it is more memory- and time-efficient by cutting redundancies in the DP matrix.
    """

    def __init__(self, tree: Tree):
        tree_indexer = TreeIndexer(tree.taxon_namespace)
        tree_indexer.index_tree(tree)  # Now all nodes in the tree have an 'index' field. First n indices = leaves.

        self.nodes = list(tree.nodes())
        self.nnodes = len(self.nodes)
        self.nleaves = len(tree.leaf_nodes())
        self.tree = tree

        # All arrays bellow will be indexed by node ids.
        # Distance from a node to a leaf:
        self.distance_lookup = np.zeros((self.nnodes, self.nleaves), dtype=np.float64)
        # Index of a leaf in the node's sorted list:
        self.index_lookup = np.zeros((self.nnodes, self.nleaves), dtype=np.int32)
        # Is the first node an ancestor of the second node:
        self.is_ancestor = np.full((self.nnodes, self.nnodes), False, dtype=np.bool_)
        # Number of leaves in the subtree rooted at each node:
        self.subtree_leaves = np.zeros(self.nnodes, dtype=np.int32)
        # Ids of two children nodes for each node ([-1, -1] if a leaf):
        self.children = np.full((self.nnodes, 2), -1, dtype=np.int32)
        # Ids of parents nodes and lengths of parent edges:
        self.parents = np.full(self.nnodes, -1, dtype=np.int32)
        self.edge_lens = np.zeros(self.nnodes, dtype=np.float64)
        # Leaf lists for each node sorted by distance (+ additional properties - see Tamir 1996):
        self.leaf_lists = np.zeros((self.nnodes, self.nleaves), dtype=np.int32)
        # Costs array:
        self.costs = np.zeros(self.nleaves, dtype=np.float64)

        self.n_c = len(self.tree.leaf_nodes())
        self.G = np.array([], dtype=np.float64)
        self.F = np.array([], dtype=np.float64)

    def find_medoids(self, p: int, distance_functions: Dict[Node, DistFunction], cost_map: Dict[str, float]):
        parnas_logger.debug(f'Started find_medoids {datetime.now().strftime("%H:%M:%S")}')
        self.distance_functions = distance_functions
        self.cost_map = cost_map
        # Parameters of a distance function for each leaf stored in an array:
        self.dist_func_array = np.empty((self.nleaves, 5), dtype=np.float64)
        # Number of medoids to choose:
        self.n_c = min(p, len(self.tree.leaf_nodes()))

        # Initializes all the key structures (bounded by O(n^2) space) needed for DP:
        self._initialize_lookups()
        parnas_logger.debug(f'Finished initialization {datetime.now().strftime("%H:%M:%S")}')

        postorder = np.array([node.index for node in self.tree.postorder_node_iter()], dtype=np.int32)
        pmedian_jit = PMedianDP(self.children, postorder, self.is_ancestor, self.leaf_lists, self.subtree_leaves,
                                self.distance_lookup, self.index_lookup, self.n_c, self.nleaves, self.dist_func_array,
                                self.costs, self.parents, self.edge_lens, self.tree.seed_node.index)
        self.G, self.F = pmedian_jit.run_dp()
        parnas_logger.debug(f'Finished DP {datetime.now().strftime("%H:%M:%S")}')

        root_ind = self.tree.seed_node.index
        min_index = int(np.argmin(self.G[root_ind][self.n_c]))
        if self.G[root_ind][self.n_c][min_index] < self.G[root_ind][0][self.nleaves - 1]:
            obj_value = self.G[root_ind][self.n_c][min_index]
            radius_id = self.leaf_lists[root_ind, min_index]
            median_names = self.backtrack(root_ind, self.n_c, radius_id)
            # median_names = [node.taxon.label for node in median_nodes]
            parnas_logger.debug(f'Finished Backtracking {datetime.now().strftime("%H:%M:%S")}')
        else:
            # No new medoids are really needed in this case.
            obj_value = self.G[root_ind][0][self.nleaves - 1]
            median_names = []

        parnas_logger.debug(f'Optimal objective function: {obj_value}')
        return obj_value, median_names

    def get_score(self, k: int) -> Optional[float]:
        root_ind = self.tree.seed_node.index
        if not self.G or len(self.G) <= 0 or k >= len(self.G[root_ind]):
            return None
        else:
            return self.G[root_ind][k][-1]

    def _is_leaf(self, node_id):
        return self.children[node_id, 0] < 0

    def _initialize_lookups(self):
        for node1 in self.tree.postorder_node_iter():
            node_dist_pairs = []
            if node1.parent_node:
                self.edge_lens[node1.index] = node1.edge_length
                self.parents[node1.index] = node1.parent_node.index

            # Fill out dist_func and costs array.
            if node1.is_leaf():
                self.costs[node1.index] = self.cost_map[node1.taxon.label]
                if self.distance_functions:
                    func = self.distance_functions[node1]
                    if func.is_zero:
                        self.dist_func_array[node1.index, 0] = 0
                    else:
                        self.dist_func_array[node1.index, 0] = 1
                        self.dist_func_array[node1.index, 1:] = [func.min_dist, func.max_dist, func.weight,
                                                                 1 if func.is_binary else 0]

            # Fill out the children array.
            if not node1.is_leaf():
                left, right = node1.child_nodes()
                self.children[node1.index, 0] = left.index
                self.children[node1.index, 1] = right.index

            # # Compute distances from node1 to all other nodes:
            # for node2, dist in dfs_tree_traversal(node1):
            #     if node2.is_leaf():
            #         self.distance_lookup[node1.index, node2.index] = dist
            #
            # # for node2 in self.tree.find_clades(lambda x: node1.is_parent_of(x), order='postorder'):
            # subtree_size = 0
            # for node2 in node1.postorder_iter():
            #     self.is_ancestor[node1.index, node2.index] = True
            #     if node2.is_leaf():
            #         node_dist_pairs.append((node2.index, self.distance_lookup[node1.index, node2.index]))
            #         subtree_size += 1
            # self.subtree_leaves[node1.index] = subtree_size
            # # for node2 in self.tree.find_clades(lambda x: not node1.is_parent_of(x), order='postorder'):
            # for node2 in filtered_postorder_iterator(self.tree, lambda v: v is not node1):
            #     if node2.is_leaf():
            #         node_dist_pairs.append((node2.index, self.distance_lookup[node1.index, node2.index]))
            # node_dist_pairs.sort(key=lambda r: r[1])  # sort pairs by distance.
            #
            # self.leaf_lists[node1.index, :] = np.array([node2_id for node2_id, dist in node_dist_pairs], dtype=np.int32)
            #
            # # Mapping a node to its index in the sorted list for node1:
            # for i, (node2_id, dist) in enumerate(node_dist_pairs):
            #     self.index_lookup[node1.index, node2_id] = i
            # # self.distance_lookup[node1.index] = dict([(node, j) for j, (node, dist) in enumerate(node_dist_pairs)])

    def backtrack(self, root_index: int, n_c: int, radius_index: int):
        median_ids = []
        medians = []
        self._backtrack_links(True, root_index, n_c, radius_index, median_ids)
        median_ids = set(median_ids)
        for leaf in self.tree.leaf_nodes():
            if leaf.index in median_ids:
                medians.append(leaf.taxon.label)
        return medians

    def _get_DP_by_id(self, node_id, q, radius_id, in_G):
        radius_list_index = self.index_lookup[node_id, radius_id]
        return self._get_DP_by_list_index(node_id, q, radius_list_index, in_G)

    def _get_DP_by_list_index(self, node_id, q, radius_list_index, in_G):
        source = self.G if in_G else self.F
        if radius_list_index < len(source[node_id][q]):
            return source[node_id][q][radius_list_index]
        else:
            return source[node_id][q][-1]

    def _backtrack_links(self, in_G: bool, node_ind: int, q: int, radius_index: int, medians: List[int]):
        backtrack_queue = deque()
        backtrack_queue.append((in_G, node_ind, q, radius_index))
        while backtrack_queue:
            in_G, node_ind, q, radius_index = backtrack_queue.popleft()
            if not in_G and self.index_lookup[node_ind, radius_index] >= len(self.F[node_ind][q]):
                radius_index = self.leaf_lists[node_ind, len(self.F[node_ind][q]) - 1]
            if q <= 0:
                continue
            if self.children[node_ind, 0] < 0:
                # We are at a leaf.
                if in_G:
                    medians.append(node_ind)
                elif q == 1:
                    rad_list_index = self.index_lookup[node_ind, radius_index]
                    if self.F[node_ind][1][rad_list_index] == self.G[node_ind][1][rad_list_index]:
                        medians.append(node_ind)
            else:
                # Internal node.
                if in_G:
                    cur_val = self._get_DP_by_id(node_ind, q, radius_index, True)
                    #cur_val = self.G[node_ind][q][self.index_lookup[node_ind, radius_index]]
                    if self.index_lookup[node_ind, radius_index] > 0 and\
                            cur_val == self._get_DP_by_list_index(node_ind, q, self.index_lookup[node_ind, radius_index] - 1, True):
                            # cur_val == self.G[node_ind][q][self.index_lookup[node_ind, radius_index] - 1]:
                        next_rad_id = self.leaf_lists[node_ind, self.index_lookup[node_ind, radius_index] - 1]
                        backtrack_queue.append((True, node_ind, q, next_rad_id))
                    else:
                        left, right = self.children[node_ind]
                        c1, c2 = (left, right) if (self.is_ancestor[left, radius_index]) else (right, left)
                        c1_size = self.subtree_leaves[c1]
                        c2_size = self.subtree_leaves[c2]
                        for q1 in range(max(1, q - c2_size), min(c1_size + 1, q + 1)):
                            q2 = q - q1
                            ch_val = self._get_DP_by_id(c1, q1, radius_index, True) +\
                                     self._get_DP_by_id(c2, q2, radius_index, False)
                                     # self.G[c1][q1][self.index_lookup[c1, radius_index]] +\
                                     # self.F[c2][q2][self.index_lookup[c2, radius_index]]
                            ch_val += 0  # Assuming that the contribution of internal nodes is 0!
                            if ch_val == cur_val:
                                # Follow links into the children.
                                backtrack_queue.append((True, c1, q1, radius_index))
                                backtrack_queue.append((False, c2, q2, radius_index))
                                break
                else:
                    # In F:
                    cur_val = self._get_DP_by_id(node_ind, q, radius_index, False)
                    # cur_val = self.F[node_ind][q][self.index_lookup[node_ind, radius_index]]
                    # if cur_val == self.G[node_ind][q][self.index_lookup[node_ind, radius_index]]:
                    if cur_val == self._get_DP_by_id(node_ind, q, radius_index, True):
                        backtrack_queue.append((True, node_ind, q, radius_index))
                    else:
                        left, right = self.children[node_ind]
                        left_size = self.subtree_leaves[left]
                        right_size = self.subtree_leaves[right]
                        for q1 in range(max(0, q - right_size), min(left_size + 1, q + 1)):
                            q2 = q - q1
                            ch_val = self._get_DP_by_id(left, q1, radius_index, False) +\
                                     self._get_DP_by_id(right, q2, radius_index, False)
                            ch_val += 0  # Assuming that the contribution of internal nodes is 0!
                            if ch_val == cur_val:
                                # Follow links into the children.
                                backtrack_queue.append((False, left, q1, radius_index))
                                backtrack_queue.append((False, right, q2, radius_index))
                                break


@jitclass([('children_arr', int32[:, ::1]), ('postorder_arr', int32[::1]), ('is_ancestor_arr', bool_[:, ::1]),
           ('leaf_lists', int32[:, ::1]), ('subtree_leaves', int32[::1]), ('distance_lookup', float64[:, ::1]),
           ('index_lookup', int32[:, ::1]), ('n_c', int32), ('nleaves', int32), ('nnodes', int32),
           ('dist_func_arr', float64[:, ::1]), ('costs', float64[::1]), ('parents', int32[::1]),
           ('edge_lengths', float64[::1]), ('root_id', int32),
           ('G', types.ListType(types.ListType(float64[::1]))), ('F', types.ListType(types.ListType(float64[::1])))
           ])
class PMedianDP:
    """
    Numba-optimized dynamic programming solution for the p-medians (k-medoids) problem.
    Note that for numba typed.List data structure accessing elements by getitem_unchecked() achieves a 10-15 times
    faster performance than typical [] access.
    """

    def __init__(self, children_arr, postorder_arr, is_ancestor_arr, leaf_lists, subtree_leaves,
                 distance_lookup, index_lookup, n_c: int, nleaves: int, dist_func_arr, costs, parents, edge_lengths,
                 root_id: int):
        self.children_arr = children_arr
        self.postorder_arr = postorder_arr
        self.is_ancestor_arr = is_ancestor_arr
        self.leaf_lists = leaf_lists
        self.subtree_leaves = subtree_leaves
        self.distance_lookup = distance_lookup
        self.index_lookup = index_lookup
        self.dist_func_arr = dist_func_arr
        self.costs = costs
        self.n_c = n_c
        self.nleaves = nleaves
        self.nnodes = len(postorder_arr)
        self.parents = parents
        self.edge_lengths = edge_lengths
        self.root_id = root_id

        self._initialize_lookups()

        # This is a necessary hack to properly initialize typed lists by filling them with proper content.
        # We allocate arrays for q=0 for all nodes in the tree (q > 0 arrays will be allocated dynamically later):
        g_temp = typed.List()
        f_temp = typed.List()
        for i in range(self.nnodes):
            gq0 = typed.List()
            gq0.append(np.full(self.nleaves, np.inf, dtype=np.float64))
            g_temp.append(gq0)
            fq0 = typed.List()
            fq0.append(np.full(self.nleaves, np.inf, dtype=np.float64))
            f_temp.append(fq0)
        self.G = g_temp
        self.F = f_temp

    def run_dp(self):
        for node_id in self.postorder_arr:
            G_node = self.G.getitem_unchecked(node_id)
            F_node = self.F.getitem_unchecked(node_id)
            if self._is_leaf(node_id):
                # Allocate and append arrays for q = 1 (we do not need q > 1):
                G_node.append(np.full(self.nleaves, self.costs[node_id], dtype=np.float64))
                F_node.append(np.full(self.nleaves, np.inf, dtype=np.float64))
                self._initialize_G_and_F(node_id)
            else:
                max_q = min(self.subtree_leaves[node_id], self.n_c)
                # self.G[node_id].append(np.full(1, np.inf, dtype=np.float64))
                # self.F[node_id]np.full(self.nleaves, np.inf, dtype=np.float64)
                max_radius = self.nleaves - 1
                for q in range(1, max_q + 1):
                    while not self.is_ancestor_arr[node_id, self.leaf_lists[node_id, max_radius]]:
                        max_radius -= 1
                    G_node.append(np.full(max_radius + 1, np.inf, dtype=np.float64))
                    F_node.append(np.full(max_radius + 1, np.inf, dtype=np.float64))
                    if self.distance_lookup[node_id, self.leaf_lists[node_id, max_radius]] >\
                            self.distance_lookup[node_id, self.leaf_lists[node_id, max_radius - 1]]:
                        max_radius -= 1

                for q in range(max_q + 1):
                    for radius_list_index in range(len(G_node.getitem_unchecked(q))):
                        self._computeG(q, node_id, self.leaf_lists[node_id, radius_list_index])
                        self._computeF(q, node_id, self.leaf_lists[node_id, radius_list_index])
        return self.G, self.F

    def _is_leaf(self, node_id) -> bool:
        return self.children_arr[node_id, 0] < 0

    def _dist_func(self, node_id: int, dist: float):
        if self.dist_func_arr[node_id, 0] < 0.5:
            return 0
        else:
            min_dist, max_dist, weight, is_binary = self.dist_func_arr[node_id, 1:]
            if dist <= min_dist:
                return 0
            if min_dist > 0 and is_binary > 0.5:
                return weight
            # mapped_dist = dist if (dist <= max_dist) else max_dist  # This is an older interpretation -- obsolete.
            mapped_dist = min(dist, max_dist) - min_dist
            mapped_dist = mapped_dist * weight
            return mapped_dist

    def _get_DP_by_id(self, node_id, q, radius_id, in_G):
        radius_list_index = self.index_lookup[node_id, radius_id]
        return self._get_DP_by_list_index(node_id, q, radius_list_index, in_G)

    def _get_DP_by_list_index(self, node_id, q, radius_list_index, in_G):
        source = self.G if in_G else self.F
        dp_row = source.getitem_unchecked(node_id).getitem_unchecked(q)
        if radius_list_index < len(dp_row):
            return dp_row[radius_list_index]
        else:
            return dp_row[-1]

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

    def _initialize_G_and_F(self, node_id: int):
        self.G.getitem_unchecked(node_id)[0] = np.full(self.nleaves,  self._dist_func(node_id, np.inf))
        # self.G[node_id] = np.full(self.nleaves, 0)
        for radius_id in self.postorder_arr:
            if not self._is_leaf(radius_id) or self.is_ancestor_arr[node_id, radius_id]:
                continue  # skip internals and descendants.
            radius_node_list_ind = self.index_lookup[node_id, radius_id]
            # self.F[0, node_id, radius_node_list_ind] = self.distance_functions[node](
            #     self.distance_lookup[node.index, radius_node.index])
            self.F[node_id][0][radius_node_list_ind] = self._dist_func(node_id, self.distance_lookup[node_id, radius_id])
            self.F[node_id][1][radius_node_list_ind] = min(self.F[node_id][0][radius_node_list_ind],
                                                           self.G[node_id][1][radius_node_list_ind])

    def _computeG(self, q: int, node_id: int, radius_id: int):
        node_q_row = self.G.getitem_unchecked(node_id).getitem_unchecked(q)
        if q == 0:
            left_id, right_id = self.children_arr[node_id]
            node_q_row[self.index_lookup[node_id, radius_id]] =\
                self.G.getitem_unchecked(left_id).getitem_unchecked(0)[self.index_lookup[left_id, radius_id]] +\
                self.G.getitem_unchecked(right_id).getitem_unchecked(0)[self.index_lookup[right_id, radius_id]]
            return
        if node_id == radius_id:
            return  # skip.
        elif not self.is_ancestor_arr[node_id, radius_id]:
            if self.index_lookup[node_id, radius_id] > 0:
                node_q_row[self.index_lookup[node_id, radius_id]] =\
                    node_q_row[self.index_lookup[node_id, radius_id] - 1]
            else:
                node_q_row[self.index_lookup[node_id, radius_id]] = np.inf
        else:
            left_id, right_id = self.children_arr[node_id]
            if self.is_ancestor_arr[left_id, radius_id]:
                self._computeG_subtree(q, node_id, radius_id, left_id, right_id)
            else:
                self._computeG_subtree(q, node_id, radius_id, right_id, left_id)

    def _computeG_subtree(self, q: int, node_id: int, radius_id: int, c1: int, c2: int):
        c1_size = self.subtree_leaves[c1]
        c2_size = self.subtree_leaves[c2]

        min_q = max(1, q - c2_size)
        max_q = min(c1_size + 1, q + 1)
        if min_q >= max_q:
            return

        mindist = np.inf
        for q1 in range(min_q, max_q):
            q2 = q - q1
            dist = self._get_DP_by_id(c1, q1, radius_id, True) + self._get_DP_by_id(c2, q2, radius_id, False)
            # print(str(dist))
            # dist = self.G[c1.index][q1][self.index_lookup[c1.index, radius_node.index]] +\
            #        self.F[q2, c2.index, self.index_lookup[c2.index, radius_node.index]]
            if dist < mindist:
                mindist = dist
        # print('------')
        # node_distance = self.distance_functions[node](self.distance_lookup[node.index, radius_node.index]) + mindist
        node_distance = mindist  # NOTE: We assume no contribution of internal nodes!

        if self.index_lookup[node_id, radius_id] - 1 >= 0:
            self.G.getitem_unchecked(node_id).getitem_unchecked(q)[self.index_lookup[node_id, radius_id]] = min(
                node_distance,
                self._get_DP_by_list_index(node_id, q, self.index_lookup[node_id, radius_id] - 1, True)
            )
        else:
            self.G.getitem_unchecked(node_id).getitem_unchecked(q)[
                self.index_lookup[node_id, radius_id]] = node_distance

    def _computeF(self, q: int, node_id: int, radius_id: int):
        if self.is_ancestor_arr[node_id, radius_id]:
            # return  # Skip.
            self.F.getitem_unchecked(node_id).getitem_unchecked(q)[self.index_lookup[node_id, radius_id]] =\
                self.G.getitem_unchecked(node_id).getitem_unchecked(q)[self.index_lookup[node_id, radius_id]]
        left_id, right_id = self.children_arr[node_id]
        left_size = self.subtree_leaves[left_id]
        right_size = self.subtree_leaves[right_id]

        min_q = max(0, q - right_size)
        max_q = min(left_size + 1, q + 1)
        if min_q >= max_q:
            return

        mindist = np.inf
        for q1 in range(min_q, max_q):
            q2 = q - q1
            dist = self._get_DP_by_id(left_id, q1, radius_id, False) + self._get_DP_by_id(right_id, q2, radius_id, False)
            # dist = self.F[q1, left.index, self.index_lookup[left.index, radius_node.index]] + self.F[
            #     q2, right.index, self.index_lookup[right.index, radius_node.index]]
            if dist< mindist:
                mindist = dist

        node_distance = mindist  # NOTE: we assume no contribution of internal nodes!

        self.F.getitem_unchecked(node_id).getitem_unchecked(q)[self.index_lookup[node_id, radius_id]] = min(
            node_distance,
            self._get_DP_by_id(node_id, q, radius_id, True)
        )

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

# @jit('f8(i4,i4,i4,i4,i4,i4,i4,f8[:,:,:],f8[:,:,:])')
# def jit_G_DP(q, left_id, right_id, left_rad, right_rad, min_q, max_q, G, F):
#     mindist = np.inf
#     for q1 in range(min_q, max_q):
#         q2 = q - q1
#         dist = G[q1, left_id, left_rad] + F[q2, right_id, right_rad]
#         if dist < mindist:
#             mindist = dist
#     return mindist
#
#
# @jit('f8(i4,i4,i4,i4,i4,i4,i4,f8[:,:,:])')
# def jit_F_DP(q, left_id, right_id, left_rad, right_rad, min_q, max_q, F):
#     mindist = np.inf
#     for q1 in range(min_q, max_q):
#         q2 = q - q1
#         dist = F[q1, left_id, left_rad] + F[q2, right_id, right_rad]
#         if dist < mindist:
#             mindist = dist
#     return mindist
