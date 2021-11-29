# -*- coding: utf-8 -*-
from typing import List
from dendropy import Tree, Node
from Bio.Phylo import BaseTree


def annotate_with_closest_centers(tree: Tree, centers: List[str], radius=None):
    """
    For each node in the tree finds its closest center. The nodes are then annotated with the 'center' keyword.
    A simple O(kn) algorithm where k is the number of centers and n is the number of leaves in the tree.
    :param tree: input tree to be annotated
    :param centers: tip labels of the chosen representatives
    :param radius: optional: if provided, only the vertices within the radius are covered by a center.
    """

    def tree_traversal(node: Node, prev_node: Node, cur_dist: float, distances: dict):
        distances[node] = cur_dist
        neighbors = []
        if node.child_nodes():
            neighbors += node.child_nodes()
        if node.parent_node:
            neighbors.append(node.parent_node)
        for neighbor in neighbors:
            if neighbor is not prev_node:
                is_parent = node.parent_node is neighbor
                edge_len = node.edge_length if is_parent else neighbor.edge_length
                tree_traversal(neighbor, node, cur_dist + edge_len, distances)

    distances_by_center = []
    for center_label in centers:
        start_leaf = tree.find_node_with_taxon_label(center_label)
        distances = {}
        tree_traversal(start_leaf, None, 0, distances)
        distances_by_center.append(distances)

    for node in tree.nodes():
        assert isinstance(node, Node)
        dist_to_centers = [distances[node] for distances in distances_by_center]
        min_dist = min(dist_to_centers)
        if radius and min_dist >= radius:
            continue
        closest_center = dist_to_centers.index(min_dist)
        node.annotations.add_new('center', closest_center)


def build_distance_functions(tree: BaseTree, radius=None, prior_centers=None, taxa_weights=None):
    # TODO: implement handling of prior_centers.
    distance_functions = {}
    for node in tree.find_clades(order='preorder'):
        if node.is_terminal():
            weight = taxa_weights[node.name] if taxa_weights else 1  # default weight is 1.
            if radius:
                function = lambda dist: (0 if dist < radius else dist) * weight
            else:
                function = lambda dist: dist * weight
        else:
            function = lambda dist: 0  # non-terminal nodes do not contribute to the distance function.

        distance_functions[node] = function

    return distance_functions
