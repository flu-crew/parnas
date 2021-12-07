# -*- coding: utf-8 -*-
from typing import Dict, Callable
import numpy as np
from Bio.Phylo import BaseTree
from Bio.Phylo.BaseTree import Clade

from .pmedian_utils import cost_function, calculate_distance, get_tree_nodes


class PMedianFinder(object):
    """
    Class that finds the median nodes (among leaves!) given a phylogenetic tree.
    An adapted implementation of Tamir's algorithm (Tamir 1996) for generalized p-medians.
    """

    def __init__(self, tree: BaseTree.Tree):
        self.cost_function = cost_function

        self.nodes = get_tree_nodes(tree)
        self.nnodes = len(self.nodes)
        self.tree = tree

        self.index_lookup = dict([(j, i) for i, j in enumerate(self.nodes)])
        self.distance_lookup = dict()
        self.r_lookup = dict()
        self.node_lists = [None] * self.nnodes

        self.n_c = self.tree.count_terminals()
        self.G = np.array([])
        self.F = np.array([])
        self.Gmedian_nodes = np.array([])
        self.Fmedian_nodes = np.array([])

    def find_medoids(self, p: int, distance_functions: Dict[Clade, Callable[[float], float]]):
        self.distance_functions = distance_functions
        post_order_nodes = list(self.tree.find_clades(order='postorder'))
        self.n_c = min(p, self.tree.count_terminals())
        self.G = np.full((self.n_c + 1, self.nnodes, self.nnodes), np.inf)
        self.F = np.full((self.n_c + 1, self.nnodes, self.nnodes), np.inf)
        self.Gmedian_nodes = np.full((self.n_c + 1, self.nnodes, self.nnodes), set())
        self.Fmedian_nodes = np.full((self.n_c + 1, self.nnodes, self.nnodes), set())
        self._initialize_lookups()
        # self.G[0] = np.full((self.nnodes, self.nnodes), np.inf)

        for node in post_order_nodes:
            if node.is_terminal():
                self._initialize_G_and_F(node)
            else:
                for q in range(self.n_c + 1):
                    for radius_node in self.node_lists[self.index_lookup[node]]:
                        self._computeG(q, node, radius_node)
                        self._computeF(q, node, radius_node)

        if self.G[self.n_c, 0, self.nnodes - 1] < self.G[0, 0, self.nnodes - 1]:
            min_index = np.argmin(self.G[self.n_c, 0])
            obj_value = self.G[self.n_c, 0, min_index]
            median_nodes = self.Gmedian_nodes[self.n_c, 0, min_index]
            median_names = [node.name for node in median_nodes]
        else:
            obj_value = self.G[0, 0, self.nnodes - 1]
            median_names = []
        return obj_value, median_names

    def _initialize_lookups(self):
        """
        Constructs sorted node lists for each node.
        Initializes distance_lookup and r_lookup maps.
        TODO: re-implement to assert O(n^2) runtime. Currently O(n^3)!
        """
        post_order_traversal = list(self.tree.find_clades(order='postorder'))
        for node1 in post_order_traversal:
            node_dist_pairs = []

            # Note: 'is_parent_of' checks ancestry and not just direct parent.
            for node2 in self.tree.find_clades(lambda x: node1.is_parent_of(x), order='postorder'):
                node_dist_pairs.append((node2, calculate_distance(node1, node2, self.tree)))
            for node2 in self.tree.find_clades(lambda x: not node1.is_parent_of(x), order='postorder'):
                node_dist_pairs.append((node2, calculate_distance(node1, node2, self.tree)))
            node_dist_pairs.sort(key=lambda r: r[1])  # sort pairs by distance.

            self.node_lists[self.index_lookup[node1]] = [node2 for node2, dist in node_dist_pairs]
            # Maps a node (Clade) to its index in the sored list for node1:
            self.distance_lookup[node1] = dict([(node, j) for j, (node, dist) in enumerate(node_dist_pairs)])
            # Maps a node (Clade) to the respective distance from node1:
            self.r_lookup[node1] = dict(node_dist_pairs)

    def _initialize_G_and_F(self, node: Clade):
        self.G[0, self.index_lookup[node]] = np.full((self.nnodes), self.distance_functions[node](np.inf))
        self.G[1, self.index_lookup[node]] = np.full((self.nnodes), 0)
        self.Gmedian_nodes[1, self.index_lookup[node]] = np.full((self.nnodes), {node})
        for radius_node in self.tree.find_clades(lambda x: not node.is_parent_of(x)):
            self.F[0, self.index_lookup[node], self.distance_lookup[node][radius_node]] = self.distance_functions[node](
                self.r_lookup[node][radius_node])
            self.Fmedian_nodes[0, self.index_lookup[node], self.distance_lookup[node][radius_node]] = {radius_node}
            if self.F[0, self.index_lookup[node], self.distance_lookup[node][radius_node]]\
                    > self.G[1, self.index_lookup[node], self.distance_lookup[node][radius_node]]:
                self.F[1, self.index_lookup[node], self.distance_lookup[node][radius_node]] =\
                    self.G[1, self.index_lookup[node], self.distance_lookup[node][radius_node]]
                self.Fmedian_nodes[1, self.index_lookup[node], self.distance_lookup[node][radius_node]] = self.Gmedian_nodes[
                    1, self.index_lookup[node], self.distance_lookup[node][radius_node]]

            else:
                self.F[1, self.index_lookup[node], self.distance_lookup[node][radius_node]] = self.F[
                    0, self.index_lookup[node], self.distance_lookup[node][radius_node]]
                self.Fmedian_nodes[1, self.index_lookup[node], self.distance_lookup[node][radius_node]] = self.Fmedian_nodes[
                    1, self.index_lookup[node], self.distance_lookup[node][radius_node]]

    def _computeG(self, q: int, node: Clade, radius_node: Clade):
        if q == 0:
            left, right = node.clades
            self.G[q, self.index_lookup[node], self.distance_lookup[node][radius_node]] =\
                self.G[q, self.index_lookup[left], self.distance_lookup[left][radius_node]] +\
                self.G[q, self.index_lookup[right], self.distance_lookup[right][radius_node]]
            return
        # for q in range(1, self.n_c + 1):
        #     for radius_node in self.distance_lookup[node]:
        if node == radius_node:
            self.G[q, self.index_lookup[node], self.distance_lookup[node][radius_node]] = np.inf
            self.Gmedian_nodes[q, self.index_lookup[node], self.distance_lookup[node][radius_node]] = self.Gmedian_nodes[
                q - 1, self.index_lookup[node], self.distance_lookup[node][radius_node]].union({node})
        elif not node.is_parent_of(radius_node):
            self.G[q, self.index_lookup[node], self.distance_lookup[node][radius_node]] = self.G[
                q, self.index_lookup[node], self.distance_lookup[node][radius_node] - 1]
            self.Gmedian_nodes[q, self.index_lookup[node], self.distance_lookup[node][radius_node]] = self.Gmedian_nodes[
                q, self.index_lookup[node], self.distance_lookup[node][radius_node] - 1]
        else:
            left, right = node.clades
            if left.is_parent_of(radius_node):
                self._computeG_subtree(q, node, radius_node, left, right)
            else:
                self._computeG_subtree(q, node, radius_node, right, left)

    def _computeG_subtree(self, q: int, node: Clade, radius_node: Clade, c1: Clade, c2: Clade):
        n1_size = len(list(c1.find_elements()))
        n2_size = len(list(c2.find_elements()))
        mindist = np.inf
        min_q1 = 0
        min_q2 = 0
        for q1 in range(max(1, q - n2_size), min(n1_size + 1, q + 1)):
            q2 = q - q1
            dist = self.G[q1, self.index_lookup[c1], self.distance_lookup[c1][radius_node]] +\
                   self.F[q2, self.index_lookup[c2], self.distance_lookup[c2][radius_node]]
            if dist < mindist:
                mindist = dist
                min_q1 = q1
                min_q2 = q2
        node_distance = self.distance_functions[node](self.r_lookup[node][radius_node]) + mindist
        if self.G[q, self.index_lookup[node], self.distance_lookup[node][radius_node] - 1] > node_distance:
            self.G[q, self.index_lookup[node], self.distance_lookup[node][radius_node]] = node_distance
            medians = self.Gmedian_nodes[min_q1, self.index_lookup[c1], self.distance_lookup[c1][radius_node]]
            medians = medians.union(self.Fmedian_nodes[min_q2, self.index_lookup[c2], self.distance_lookup[c2][radius_node]])
            medians = medians.union({radius_node})
            self.Gmedian_nodes[q, self.index_lookup[node], self.distance_lookup[node][radius_node]] = medians
        else:
            self.G[q, self.index_lookup[node], self.distance_lookup[node][radius_node]] = self.G[
                q, self.index_lookup[node], self.distance_lookup[node][radius_node] - 1]
            self.Gmedian_nodes[q, self.index_lookup[node], self.distance_lookup[node][radius_node]] = self.Gmedian_nodes[
                q, self.index_lookup[node], self.distance_lookup[node][radius_node] - 1]

    def _computeF(self, q: int, node: Clade, radius_node: Clade):
        left, right = node.clades
        left_size = len(list(left.find_elements()))  # TODO: Precompute subtree sizes.
        right_size = len(list(right.find_elements()))
        mindist = np.inf
        min_q1 = 0
        min_q2 = 0
        for q1 in range(max(0, q - right_size), min(left_size + 1, q + 1)):
            q2 = q - q1
            dist = self.F[q1, self.index_lookup[left], self.distance_lookup[left][radius_node]] + self.F[
                q2, self.index_lookup[right], self.distance_lookup[right][radius_node]]
            if mindist > dist:
                mindist = dist
                min_q1 = q1
                min_q2 = q2
        # self.F[q, self.index_lookup[node], self.distance_lookup[node][radius_node]] = min(
        #     self.G[q, self.index_lookup[node], self.distance_lookup[node][radius_node]],
        #     self.distance_functions[node](self.r_lookup[node][radius_node]) + mindist)
        node_distance = self.distance_functions[node](self.r_lookup[node][radius_node]) + mindist
        if self.G[q, self.index_lookup[node], self.distance_lookup[node][radius_node]] > node_distance:
            self.F[q, self.index_lookup[node], self.distance_lookup[node][radius_node]] = node_distance
            medians = self.Fmedian_nodes[min_q1, self.index_lookup[left], self.distance_lookup[left][radius_node]]
            medians = medians.union(
                self.Fmedian_nodes[min_q2, self.index_lookup[right], self.distance_lookup[right][radius_node]])
            medians = medians.union({radius_node})
            self.Fmedian_nodes[q, self.index_lookup[node], self.distance_lookup[node][radius_node]] = medians
        else:
            self.F[q, self.index_lookup[node], self.distance_lookup[node][radius_node]] = self.G[
                q, self.index_lookup[node], self.distance_lookup[node][radius_node]]
            self.Fmedian_nodes[q, self.index_lookup[node], self.distance_lookup[node][radius_node]] = self.Gmedian_nodes[
                q, self.index_lookup[node], self.distance_lookup[node][radius_node]]
