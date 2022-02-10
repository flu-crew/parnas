# -*- coding: utf-8 -*-
from collections import deque
from typing import Dict, List

import numpy as np
from dendropy import Tree, Node
from datetime import datetime

from .medoid_utils import DistFunction
from .pmedian_utils import cost_function, filtered_preorder_iterator, filtered_postorder_iterator, dfs_tree_traversal
from .tree_indexer import TreeIndexer
from parnas.logging import parnas_logger


class PMedianFinder(object):
    """
    Class that finds the median nodes (among leaves) given a phylogenetic tree.
    An adapted implementation of Tamir's algorithm (Tamir 1996) for generalized p-medians.
    """

    def __init__(self, tree: Tree):
        tree_indexer = TreeIndexer(tree.taxon_namespace)
        tree_indexer.index_tree(tree)  # Now all nodes in the tree have an 'index' field. First n indices = leaves.
        self.cost_function = cost_function

        self.nodes = list(tree.nodes())
        self.nnodes = len(self.nodes)
        self.nleaves = len(tree.leaf_nodes())
        self.tree = tree

        self.distance_lookup = np.zeros((self.nnodes, self.nleaves), dtype=np.float64)
        self.index_lookup = np.zeros((self.nnodes, self.nleaves), dtype=np.int32)
        self.is_ancestor = np.full((self.nnodes, self.nnodes), False, dtype=np.bool_)
        self.subtree_leaves = np.zeros(self.nnodes, dtype=np.int32)
        self.children = np.full((self.nnodes, 2), -1, dtype=np.int32)
        self.leaf_lists = np.zeros((self.nnodes, self.nleaves), dtype=np.int32)

        self.n_c = len(self.tree.leaf_nodes())
        self.G = np.array([], dtype=np.float64)
        self.F = np.array([], dtype=np.float64)

    def find_medoids(self, p: int, distance_functions: Dict[Node, DistFunction]):
        parnas_logger.debug(f'Started find_medoids {datetime.now().strftime("%H:%M:%S")}')
        self.distance_functions = distance_functions
        self.n_c = min(p, len(self.tree.leaf_nodes()))
        self.G = np.full((self.n_c + 1, self.nnodes, self.nleaves), np.inf, dtype=np.float64)
        self.F = np.full((self.n_c + 1, self.nnodes, self.nleaves), np.inf, dtype=np.float64)
        parnas_logger.debug(f'Allocated arrays {datetime.now().strftime("%H:%M:%S")}')

        self._initialize_lookups()
        parnas_logger.debug(f'Finished initialization {datetime.now().strftime("%H:%M:%S")}')

        processed_nodes = 0
        i = 0
        for node in self.tree.postorder_node_iter():
            if node.is_leaf():
                self._initialize_G_and_F(node)
            else:
                for q in range(min(self.subtree_leaves[node.index], self.n_c) + 1):
                    for radius_id in self.leaf_lists[node.index]:
                        self._computeG(q, node.index, radius_id)
                        self._computeF(q, node.index, radius_id)
            processed_nodes += 1
            if processed_nodes % 100 == 0:
                parnas_logger.debug(f'\tProcessed {processed_nodes} nodes')
        parnas_logger.debug(f'Finished DP {datetime.now().strftime("%H:%M:%S")}')

        root_ind = self.tree.seed_node.index
        min_index = int(np.argmin(self.G[self.n_c][root_ind]))
        if self.G[self.n_c][root_ind][min_index] < self.G[0][root_ind][self.nleaves - 1]:
            obj_value = self.G[self.n_c][root_ind][min_index]
            radius_id = self.leaf_lists[root_ind, min_index]
            median_names = self.backtrack(root_ind, self.n_c, radius_id)
            # median_names = [node.taxon.label for node in median_nodes]
            parnas_logger.debug(f'Finished Backtracking {datetime.now().strftime("%H:%M:%S")}')
        else:
            obj_value = self.G[0][root_ind][self.nleaves - 1]
            median_names = []
        parnas_logger.debug(f'Optimal objective function: {obj_value}')
        return obj_value, median_names

    def _initialize_lookups(self):
        """
        Constructs sorted node lists for each node.
        Initializes distance_lookup and r_lookup maps + other helper structures.
        TODO: re-implement to assert O(n^2) runtime. Currently O(n^2 logn)!
        """
        for node1 in self.tree.postorder_node_iter():
            node_dist_pairs = []

            # Fill out the children array.
            if not node1.is_leaf():
                left, right = node1.child_nodes()
                self.children[node1.index, 0] = left.index
                self.children[node1.index, 1] = right.index

            # Compute distances from node1 to all other nodes:
            for node2, dist in dfs_tree_traversal(node1):
                if node2.is_leaf():
                    self.distance_lookup[node1.index, node2.index] = dist

            # for node2 in self.tree.find_clades(lambda x: node1.is_parent_of(x), order='postorder'):
            subtree_size = 0
            for node2 in node1.postorder_iter():
                self.is_ancestor[node1.index, node2.index] = True
                if node2.is_leaf():
                    node_dist_pairs.append((node2.index, self.distance_lookup[node1.index, node2.index]))
                    subtree_size += 1
            self.subtree_leaves[node1.index] = subtree_size
            # for node2 in self.tree.find_clades(lambda x: not node1.is_parent_of(x), order='postorder'):
            for node2 in filtered_postorder_iterator(self.tree, lambda v: v is not node1):
                if node2.is_leaf():
                    node_dist_pairs.append((node2.index, self.distance_lookup[node1.index, node2.index]))
            node_dist_pairs.sort(key=lambda r: r[1])  # sort pairs by distance.

            self.leaf_lists[node1.index, :] = np.array([node2_id for node2_id, dist in node_dist_pairs], dtype=np.int32)

            # Mapping a node to its index in the sorted list for node1:
            for i, (node2_id, dist) in enumerate(node_dist_pairs):
                self.index_lookup[node1.index, node2_id] = i
            # self.distance_lookup[node1.index] = dict([(node, j) for j, (node, dist) in enumerate(node_dist_pairs)])

    def _initialize_G_and_F(self, node: Node):
        self.G[0, node.index] = np.full(self.nleaves, self.distance_functions[node].get_dist(np.inf))
        self.G[1, node.index] = np.full(self.nleaves, 0)
        for radius_node in filtered_preorder_iterator(self.tree, subtree_filter=lambda v: v is not node):
            if not radius_node.is_leaf():
                continue  # skip internals.
            radius_node_list_ind = self.index_lookup[node.index, radius_node.index]
            self.F[0, node.index, radius_node_list_ind] = self.distance_functions[node].get_dist(
                self.distance_lookup[node.index, radius_node.index])
            if self.F[0, node.index, radius_node_list_ind]\
                    > self.G[1, node.index, radius_node_list_ind]:
                self.F[1, node.index, radius_node_list_ind] = self.G[1, node.index, radius_node_list_ind]
            else:
                self.F[1, node.index, radius_node_list_ind] = self.F[0, node.index, radius_node_list_ind]

    def _computeG(self, q: int, node_id: int, radius_id: int):
        if q == 0:
            left_id, right_id = self.children[node_id]
            self.G[q, node_id, self.index_lookup[node_id, radius_id]] =\
                self.G[q, left_id, self.index_lookup[left_id, radius_id]] +\
                self.G[q, right_id, self.index_lookup[right_id, radius_id]]
            return
        if node_id == radius_id:
            return  # skip.
        elif not self.is_ancestor[node_id, radius_id]:
            self.G[q, node_id, self.index_lookup[node_id, radius_id]] = self.G[
                q, node_id, self.index_lookup[node_id, radius_id] - 1]
        else:
            left_id, right_id = self.children[node_id]
            if self.is_ancestor[left_id, radius_id]:
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
        mindist = np.min(
            self.G[min_q:max_q, c1, self.index_lookup[c1, radius_id]] +
            np.flip(self.F[(q-max_q)+1:(q-min_q)+1, c2, self.index_lookup[c2, radius_id]])
        )
        # node_distance = self.distance_functions[node](self.distance_lookup[node.index, radius_node.index]) + mindist
        node_distance = mindist  # NOTE: We assume no contribution of internal nodes!

        if self.G[q, node_id, self.index_lookup[node_id, radius_id] - 1] > node_distance:
            self.G[q, node_id, self.index_lookup[node_id, radius_id]] = node_distance
        else:
            self.G[q, node_id, self.index_lookup[node_id, radius_id]] = self.G[
                q, node_id, self.index_lookup[node_id, radius_id] - 1]

    def _computeF(self, q: int, node_id: int, radius_id: int):
        if self.is_ancestor[node_id, radius_id]:
            return  # Skip.
        left_id, right_id = self.children[node_id]
        left_size = self.subtree_leaves[left_id]
        right_size = self.subtree_leaves[right_id]

        min_q = max(0, q - right_size)
        max_q = min(left_size + 1, q + 1)
        if min_q >= max_q:
            return
        mindist = np.min(
            self.F[min_q:max_q, left_id, self.index_lookup[left_id, radius_id]] +
            np.flip(self.F[(q-max_q)+1:(q-min_q)+1, right_id, self.index_lookup[right_id, radius_id]])
        )
        # self.F[q, self.index_lookup[node], self.distance_lookup[node][radius_node]] = min(
        #     self.G[q, self.index_lookup[node], self.distance_lookup[node][radius_node]],
        #     self.distance_functions[node](self.r_lookup[node][radius_node]) + mindist)
        # node_distance = self.distance_functions[node](self.distance_lookup[node.index, radius_node.index]) + mindist

        node_distance = mindist  # NOTE: we assume no contribution of internal nodes!

        if self.G[q, node_id, self.index_lookup[node_id, radius_id]] > node_distance:
            self.F[q, node_id, self.index_lookup[node_id, radius_id]] = node_distance
        else:
            self.F[q, node_id, self.index_lookup[node_id, radius_id]] = self.G[
                q, node_id, self.index_lookup[node_id, radius_id]]

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
        if radius_list_index < len(source[q, node_id]):
            return source[q, node_id, radius_list_index]
        else:
            return source[q, node_id][-1]

    def _backtrack_links(self, in_G: bool, node_ind: int, q: int, radius_index: int, medians: List[int]):
        backtrack_queue = deque()
        backtrack_queue.append((in_G, node_ind, q, radius_index))
        while backtrack_queue:
            in_G, node_ind, q, radius_index = backtrack_queue.popleft()
            if not in_G and self.index_lookup[node_ind, radius_index] >= len(self.F[q, node_ind]):
                radius_index = self.leaf_lists[node_ind, len(self.F[q, node_ind]) - 1]
            if q <= 0:
                continue
            if self.children[node_ind, 0] < 0:
                # We are at a leaf.
                if in_G:
                    medians.append(node_ind)
                elif q == 1:
                    rad_list_index = self.index_lookup[node_ind, radius_index]
                    if self.F[1][node_ind][rad_list_index] == self.G[1][node_ind][rad_list_index]:
                        medians.append(node_ind)
            else:
                # Internal node.
                if in_G:
                    cur_val = self.G[q][node_ind][self.index_lookup[node_ind, radius_index]]
                    if self.index_lookup[node_ind, radius_index] > 0 and\
                            cur_val == self.G[q][node_ind][self.index_lookup[node_ind, radius_index] - 1]:
                        next_rad_id = self.leaf_lists[node_ind, self.index_lookup[node_ind, radius_index] - 1]
                        backtrack_queue.append((True, node_ind, q, next_rad_id))
                    else:
                        left, right = self.children[node_ind]
                        c1, c2 = (left, right) if (self.is_ancestor[left, radius_index]) else (right, left)
                        c1_size = self.subtree_leaves[c1]
                        c2_size = self.subtree_leaves[c2]
                        for q1 in range(max(1, q - c2_size), min(c1_size + 1, q + 1)):
                            q2 = q - q1
                            ch_val = self.G[q1][c1][self.index_lookup[c1, radius_index]] +\
                                     self._get_DP_by_id(c2, q2, radius_index, False)
                                     # self.F[c2][q2][self.index_lookup[c2, radius_index]]
                            ch_val += 0  # Assuming that the contribution of internal nodes is 0!
                            if ch_val == cur_val:
                                # Follow links into the children.
                                backtrack_queue.append((True, c1, q1, radius_index))
                                backtrack_queue.append((False, c2, q2, radius_index))
                                break
                else:
                    # In F:
                    cur_val = self.F[q][node_ind][self.index_lookup[node_ind, radius_index]]
                    if cur_val == self.G[q][node_ind][self.index_lookup[node_ind, radius_index]]:
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
