# -*- coding: utf-8 -*-
import numpy as np
from Bio.Phylo import BaseTree
from Bio.Phylo.BaseTree import Clade

from pmedian_utils import cost_function, distance_function, calculate_distance, get_tree_nodes


class PMedianFinder(object):
    """
    Class that finds the median nodes (among leaves!) given a phylogenetic tree.
    An adapted implementation of Tamir's algorithm (Tamir 1996) for generalized p-medians.
    """

    def __init__(self, tree: BaseTree.Tree):
        self.cost_function = cost_function
        self.distance_function = distance_function

        self.nodes = get_tree_nodes(tree)
        self.nnodes = len(self.nodes)
        self.tree = tree

        self.index_lookup = dict([(j, i) for i, j in enumerate(self.nodes)])
        self.distance_lookup = dict()
        self.r_lookup = dict()

        self.n_c = self.tree.count_terminals()
        self.G = np.array([])
        self.F = np.array([])
        self.Gmedian_nodes = np.array([])
        self.Fmedian_nodes = np.array([])

    def find_medoids(self, p: int):
        post_order_nodes = list(self.tree.find_clades(order='postorder'))
        self.n_c = min(p, self.tree.count_terminals())
        self.G = np.zeros((self.n_c + 1, self.nnodes, self.nnodes))
        self.F = np.zeros((self.n_c + 1, self.nnodes, self.nnodes))
        self.Gmedian_nodes = np.full((self.n_c + 1, self.nnodes, self.nnodes), set())
        self.Fmedian_nodes = np.full((self.n_c + 1, self.nnodes, self.nnodes), set())
        self.initialize_lookups()
        self.G[0] = np.full((self.nnodes, self.nnodes), np.inf)
        for i in post_order_nodes:
            if i.is_terminal():
                self.initialize_G_and_F(i)
            else:
                for q in range(self.n_c + 1):
                    for j in self.distance_lookup[i]:
                        self.computeF(q, i, j)
                        if q > 0:
                            self.computeG(q, i, j)

        min_index = np.argmin(self.G[self.n_c, 0])
        return self.G[self.n_c, 0, min_index], self.Gmedian_nodes[self.n_c, 0, min_index]

    def initialize_lookups(self):
        post_order_traversal = list(self.tree.find_clades(order='postorder'))
        for i in post_order_traversal:
            self.distance_lookup[i] = []
            self.r_lookup[i] = []

            for j in self.tree.find_clades(lambda x: i.is_parent_of(x),
                                           order='postorder'):  # is_parent_of checks ancestry and not just direct parent
                self.distance_lookup[i].append((j, calculate_distance(i, j, self.tree)))
                self.r_lookup[i].append((j, calculate_distance(i, j, self.tree)))

            for j in self.tree.find_clades(lambda x: not i.is_parent_of(x), order='postorder'):
                self.distance_lookup[i].append((j, calculate_distance(i, j, self.tree)))
                self.r_lookup[i].append((j, calculate_distance(i, j, self.tree)))

            self.distance_lookup[i].sort(key=lambda r: r[1])
            self.r_lookup[i].sort(key=lambda r: r[1])
            self.distance_lookup[i] = dict([(i[0], j) for j, i in enumerate(self.distance_lookup[i])])
            self.r_lookup[i] = dict(self.r_lookup[i])

    def initialize_G_and_F(self, node: Clade):
        self.G[1, self.index_lookup[node]] = np.full((self.nnodes), 0)
        self.Gmedian_nodes[1, self.index_lookup[node]] = np.full((self.nnodes), {node})
        for j in self.tree.find_clades(lambda x: not node.is_parent_of(x)):
            self.F[0, self.index_lookup[node], self.distance_lookup[node][j]] = distance_function(
                calculate_distance(node, j, self.tree))
            self.Fmedian_nodes[0, self.index_lookup[node], self.distance_lookup[node][j]] = {j}
            if self.F[0, self.index_lookup[node], self.distance_lookup[node][j]]\
                    > self.G[1, self.index_lookup[node], self.distance_lookup[node][j]]:
                self.F[1, self.index_lookup[node], self.distance_lookup[node][j]] =\
                    self.G[1, self.index_lookup[node], self.distance_lookup[node][j]]
                self.Fmedian_nodes[1, self.index_lookup[node], self.distance_lookup[node][j]] = self.Gmedian_nodes[
                    1, self.index_lookup[node], self.distance_lookup[node][j]]

            else:
                self.F[1, self.index_lookup[node], self.distance_lookup[node][j]] = self.F[
                    0, self.index_lookup[node], self.distance_lookup[node][j]]
                self.Fmedian_nodes[1, self.index_lookup[node], self.distance_lookup[node][j]] = self.Fmedian_nodes[
                    1, self.index_lookup[node], self.distance_lookup[node][j]]

    def computeG(self, q: int, node: Clade, radius_node: Clade):
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
            left_size = len(list(left.find_elements()))
            right_size = len(list(right.find_elements()))
            if left.is_parent_of(radius_node):
                self.computeG_subtree(q, node, radius_node, left, right)
            else:
                self.computeG_subtree(q, node, radius_node, right, left)

    def computeG_subtree(self, q: int, node: Clade, radius_node: Clade, c1: Clade, c2: Clade):
        n1_size = len(list(c1.find_elements()))
        n2_size = len(list(c2.find_elements()))
        mindist = np.inf
        min_q1 = 0
        min_q2 = 0
        for q1 in range(max(1, q - n2_size), min(n1_size + 1, q + 1)):
            q2 = q - q1
            dist = self.G[q1, self.index_lookup[c1], self.distance_lookup[c1][radius_node]] + \
                   self.F[q2, self.index_lookup[c2], self.distance_lookup[c2][radius_node]]
            if (mindist > dist):
                mindist = dist
                min_q1 = q1
                min_q2 = q2
        node_distance = distance_function(self.r_lookup[node][radius_node]) + mindist
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

    def computeF(self, q: int, node, radius_node):
        left, right = node.clades
        left_size = len(list(left.find_elements()))
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
        self.F[q, self.index_lookup[node], self.distance_lookup[node][radius_node]] = min(
            self.G[q, self.index_lookup[node], self.distance_lookup[node][radius_node]],
            distance_function(self.r_lookup[node][radius_node]) + mindist)
        node_distance = distance_function(self.r_lookup[node][radius_node]) + mindist
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
