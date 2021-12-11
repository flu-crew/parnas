# -*- coding: utf-8 -*-
from collections import deque
from typing import Dict, Callable
import numpy as np
from dendropy import Tree, Node

from .pmedian_utils import cost_function, filtered_preorder_iterator, filtered_postorder_iterator, dfs_tree_traversal
from .tree_indexer import TreeIndexer


class PMedianFinder(object):
    """
    Class that finds the median nodes (among leaves!) given a phylogenetic tree.
    An adapted implementation of Tamir's algorithm (Tamir 1996) for generalized p-medians.
    """

    def __init__(self, tree: Tree):
        tree_indexer = TreeIndexer(tree.taxon_namespace)
        tree_indexer.index_tree(tree)  # Now all nodes in the tree have an 'index'.
        self.cost_function = cost_function

        self.nodes = list(tree.nodes())
        self.nnodes = len(self.nodes)
        self.tree = tree

        # self.index_lookup = dict([(j, i) for i, j in enumerate(self.nodes)])
        self.distance_lookup = np.zeros((self.nnodes, self.nnodes), dtype=np.float64)
        self.index_lookup = np.zeros((self.nnodes, self.nnodes), dtype=np.int32)
        self.is_ancestor = np.full((self.nnodes, self.nnodes), False, dtype=np.bool_)
        self.subtree_size = np.zeros(self.nnodes, dtype=np.int32)
        self.node_lists = [None] * self.nnodes

        self.n_c = len(self.tree.leaf_nodes())
        self.G = np.array([])
        self.F = np.array([])
        self.G_back_refs = None
        self.F_back_refs = None
        # self.Gmedian_nodes = np.array([])
        # self.Fmedian_nodes = np.array([])

    def find_medoids(self, p: int, distance_functions: Dict[Node, Callable[[float], float]]):
        self.distance_functions = distance_functions
        post_order_nodes = list(self.tree.postorder_node_iter())
        self.n_c = min(p, len(self.tree.leaf_nodes()))
        self.G = np.full((self.n_c + 1, self.nnodes, self.nnodes), np.inf, dtype=np.float64)
        self.F = np.full((self.n_c + 1, self.nnodes, self.nnodes), np.inf, dtype=np.float64)

        self.G_back_refs = np.empty((self.n_c + 1, self.nnodes, self.nnodes), dtype=object)
        self.F_back_refs = np.empty((self.n_c + 1, self.nnodes, self.nnodes), dtype=object)
        for i in range(self.n_c + 1):
            for j in range(self.nnodes):
                for k in range(self.nnodes):
                    self.G_back_refs[i, j, k] = BackReference()
                    self.F_back_refs[i, j, k] = BackReference()

        # self.Gmedian_nodes = np.full((self.n_c + 1, self.nnodes, self.nnodes), set())
        # self.Fmedian_nodes = np.full((self.n_c + 1, self.nnodes, self.nnodes), set())
        self._initialize_lookups()
        # self.G[0] = np.full((self.nnodes, self.nnodes), np.inf)

        for node in post_order_nodes:
            if node.is_leaf():
                self._initialize_G_and_F(node)
            else:
                for q in range(self.n_c + 1):
                    for radius_node in self.node_lists[node.index]:
                        self._computeG(q, node, radius_node)
                        self._computeF(q, node, radius_node)

        root_ind = self.tree.seed_node.index
        if self.G[self.n_c, root_ind, self.nnodes - 1] < self.G[0, root_ind, self.nnodes - 1]:
            min_index = np.argmin(self.G[self.n_c, root_ind])
            obj_value = self.G[self.n_c, root_ind, min_index]
            median_nodes = self._backtrack(self.G_back_refs[self.n_c, root_ind, min_index])
            median_names = [node.taxon.label for node in median_nodes]
        else:
            obj_value = self.G[0, root_ind, self.nnodes - 1]
            median_names = []
        return obj_value, median_names

    def _initialize_lookups(self):
        """
        Constructs sorted node lists for each node.
        Initializes distance_lookup and r_lookup maps.
        TODO: re-implement to assert O(n^2) runtime. Currently O(n^3)!
        """
        for node1 in self.tree.postorder_node_iter():
            node_dist_pairs = []

            # Compute distances from node1 to all other nodes:
            node_dists = []
            dfs_tree_traversal(node1, None, 0, node_dists)
            for node2, dist in node_dists:
                self.distance_lookup[node1.index, node2.index] = dist

            # for node2 in self.tree.find_clades(lambda x: node1.is_parent_of(x), order='postorder'):
            subtree_size = 0
            for node2 in node1.postorder_iter():
                node_dist_pairs.append((node2, self.distance_lookup[node1.index, node2.index]))
                self.is_ancestor[node1.index, node2.index] = True
                subtree_size += 1
            self.subtree_size[node1.index] = subtree_size
            # for node2 in self.tree.find_clades(lambda x: not node1.is_parent_of(x), order='postorder'):
            for node2 in filtered_postorder_iterator(self.tree, lambda v: v is not node1):
                node_dist_pairs.append((node2, self.distance_lookup[node1.index, node2.index]))
            node_dist_pairs.sort(key=lambda r: r[1])  # sort pairs by distance.

            self.node_lists[node1.index] = [node2 for node2, dist in node_dist_pairs]

            # Mapping a node to its index in the sorted list for node1:
            for i, (node2, dist) in enumerate(node_dist_pairs):
                self.index_lookup[node1.index, node2.index] = i
            # self.distance_lookup[node1.index] = dict([(node, j) for j, (node, dist) in enumerate(node_dist_pairs)])

    def _initialize_G_and_F(self, node: Node):
        self.G[0, node.index] = np.full(self.nnodes, self.distance_functions[node](np.inf))
        self.G[1, node.index] = np.full(self.nnodes, 0)
        for i in range(self.nnodes):
            self.G_back_refs[1, node.index, i].center = node
        # self.Gmedian_nodes[1, node.index] = np.full(self.nnodes, {node})
        for radius_node in filtered_preorder_iterator(self.tree, subtree_filter=lambda v: v is not node):
            radius_node_list_ind = self.index_lookup[node.index, radius_node.index]
            self.F[0, node.index, radius_node_list_ind] = self.distance_functions[node](
                self.distance_lookup[node.index, radius_node.index])
            self.F_back_refs[0, node.index, radius_node_list_ind].center = radius_node
            # self.Fmedian_nodes[0, node.index, radius_node_list_ind] = {radius_node}
            if self.F[0, node.index, radius_node_list_ind]\
                    > self.G[1, node.index, radius_node_list_ind]:
                self.F[1, node.index, radius_node_list_ind] = self.G[1, node.index, radius_node_list_ind]
                self.F_back_refs[1, node.index, radius_node_list_ind].add_ref(1, node.index, radius_node_list_ind, True)
                # self.Fmedian_nodes[1, node.index, radius_node_list_ind] = self.Gmedian_nodes[
                #     1, node.index, radius_node_list_ind]
            else:
                self.F[1, node.index, radius_node_list_ind] = self.F[0, node.index, radius_node_list_ind]
                self.F_back_refs[1, node.index, radius_node_list_ind].add_ref(0, node.index, radius_node_list_ind, False)
                # self.Fmedian_nodes[1, node.index, radius_node_list_ind] = self.Fmedian_nodes[
                #     0, node.index, radius_node_list_ind]

    def _computeG(self, q: int, node: Node, radius_node: Node):
        if q == 0:
            left, right = node.child_nodes()
            self.G[q, node.index, self.index_lookup[node.index, radius_node.index]] =\
                self.G[q, left.index, self.index_lookup[left.index, radius_node.index]] +\
                self.G[q, right.index, self.index_lookup[right.index, radius_node.index]]
            return
        if node == radius_node:
            self.G[q, node.index, self.index_lookup[node.index, radius_node.index]] = np.inf
            # self.Gmedian_nodes[q, node.index, self.index_lookup[node.index, radius_node.index]] = self.Gmedian_nodes[
            #     q - 1, node.index, self.index_lookup[node.index, radius_node.index]].union({node})
        # elif not node.is_parent_of(radius_node):
        elif not self.is_ancestor[node.index, radius_node.index]:
            self.G[q, node.index, self.index_lookup[node.index, radius_node.index]] = self.G[
                q, node.index, self.index_lookup[node.index, radius_node.index] - 1]
            self.G_back_refs[q, node.index, self.index_lookup[node.index, radius_node.index]].add_ref(
                q, node.index, self.index_lookup[node.index, radius_node.index] - 1, True)
            # self.Gmedian_nodes[q, node.index, self.index_lookup[node.index, radius_node.index]] = self.Gmedian_nodes[
            #     q, node.index, self.index_lookup[node.index, radius_node.index] - 1]
        else:
            left, right = node.child_nodes()
            if self.is_ancestor[left.index, radius_node.index]:
                self._computeG_subtree(q, node, radius_node, left, right)
            else:
                self._computeG_subtree(q, node, radius_node, right, left)

    def _computeG_subtree(self, q: int, node: Node, radius_node: Node, c1: Node, c2: Node):
        c1_size = self.subtree_size[c1.index]
        c2_size = self.subtree_size[c2.index]
        mindist = np.inf
        min_q1 = 0
        min_q2 = 0
        for q1 in range(max(1, q - c2_size), min(c1_size + 1, q + 1)):
            q2 = q - q1
            dist = self.G[q1, c1.index, self.index_lookup[c1.index, radius_node.index]] +\
                   self.F[q2, c2.index, self.index_lookup[c2.index, radius_node.index]]
            if dist < mindist:
                mindist = dist
                min_q1 = q1
                min_q2 = q2
        node_distance = self.distance_functions[node](self.distance_lookup[node.index, radius_node.index]) + mindist
        back_reference = self.G_back_refs[q, node.index, self.index_lookup[node.index, radius_node.index]]
        if self.G[q, node.index, self.index_lookup[node.index, radius_node.index] - 1] > node_distance:
            self.G[q, node.index, self.index_lookup[node.index, radius_node.index]] = node_distance
            back_reference.add_ref(min_q1, c1.index, self.index_lookup[c1.index, radius_node.index], True)
            back_reference.add_ref(min_q2, c2.index, self.index_lookup[c2.index, radius_node.index], False)
            back_reference.center = radius_node

            # medians = self.Gmedian_nodes[min_q1, c1.index, self.index_lookup[c1.index, radius_node.index]]
            # medians = medians.union(self.Fmedian_nodes[min_q2, c2.index, self.index_lookup[c2.index, radius_node.index]])
            # medians = medians.union({radius_node})
            # self.Gmedian_nodes[q, node.index, self.index_lookup[node.index, radius_node.index]] = medians
        else:
            self.G[q, node.index, self.index_lookup[node.index, radius_node.index]] = self.G[
                q, node.index, self.index_lookup[node.index, radius_node.index] - 1]
            back_reference.add_ref(q, node.index, self.index_lookup[node.index, radius_node.index] - 1, True)

    def _computeF(self, q: int, node: Node, radius_node: Node):
        left, right = node.child_nodes()
        left_size = self.subtree_size[left.index]
        right_size = self.subtree_size[right.index]
        mindist = np.inf
        min_q1 = 0
        min_q2 = 0
        for q1 in range(max(0, q - right_size), min(left_size + 1, q + 1)):
            q2 = q - q1
            dist = self.F[q1, left.index, self.index_lookup[left.index, radius_node.index]] + self.F[
                q2, right.index, self.index_lookup[right.index, radius_node.index]]
            if mindist > dist:
                mindist = dist
                min_q1 = q1
                min_q2 = q2
        # self.F[q, self.index_lookup[node], self.distance_lookup[node][radius_node]] = min(
        #     self.G[q, self.index_lookup[node], self.distance_lookup[node][radius_node]],
        #     self.distance_functions[node](self.r_lookup[node][radius_node]) + mindist)
        back_reference = self.F_back_refs[q, node.index, self.index_lookup[node.index, radius_node.index]]
        node_distance = self.distance_functions[node](self.distance_lookup[node.index, radius_node.index]) + mindist
        if self.G[q, node.index, self.index_lookup[node.index, radius_node.index]] > node_distance:
            self.F[q, node.index, self.index_lookup[node.index, radius_node.index]] = node_distance
            back_reference.add_ref(min_q1, left.index, self.index_lookup[left.index, radius_node.index], False)
            back_reference.add_ref(min_q2, right.index, self.index_lookup[right.index, radius_node.index], False)
            back_reference.center = radius_node

            # medians = self.Fmedian_nodes[min_q1, left.index, self.index_lookup[left.index, radius_node.index]]
            # medians = medians.union(
            #     self.Fmedian_nodes[min_q2, right.index, self.index_lookup[right.index, radius_node.index]])
            # medians = medians.union({radius_node})
            # self.Fmedian_nodes[q, node.index, self.index_lookup[node.index, radius_node.index]] = medians
        else:
            self.F[q, node.index, self.index_lookup[node.index, radius_node.index]] = self.G[
                q, node.index, self.index_lookup[node.index, radius_node.index]]
            back_reference.add_ref(q, node.index, self.index_lookup[node.index, radius_node.index], True)
            # self.Fmedian_nodes[q, node.index, self.index_lookup[node.index, radius_node.index]] = self.Gmedian_nodes[
            #     q, node.index, self.index_lookup[node.index, radius_node.index]]

    def _backtrack(self, start_ref):
        medoids = set()
        ref_queue = deque()
        ref_queue.append(start_ref)
        while ref_queue:
            back_ref = ref_queue.popleft()
            if back_ref.center:
                medoids.add(back_ref.center)
            for next_ref in back_ref.get_next_references(self.G_back_refs, self.F_back_refs):
                ref_queue.append(next_ref)
        return medoids


class BackReference:

    def __init__(self, center=None):
        self.qs = []
        self.node_inds = []
        self.radius_node_inds = []
        self.for_Gs = []
        self.center = center

    def add_ref(self, q: int, node_ind: int, radius_node_ind: int, for_G: bool):
        self.qs.append(q)
        self.node_inds.append(node_ind)
        self.radius_node_inds.append(radius_node_ind)
        self.for_Gs.append(for_G)

    def get_next_references(self, G_back_refs, F_back_refs):
        next_refs = []
        for i in range(len(self.qs)):
            next_ref = G_back_refs[self.qs[i], self.node_inds[i], self.radius_node_inds[i]] if self.for_Gs[i] else \
                F_back_refs[self.qs[i], self.node_inds[i], self.radius_node_inds[i]]
            next_refs.append(next_ref)
        return next_refs
