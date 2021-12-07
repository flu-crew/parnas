# -*- coding: utf-8 -*-
from typing import List, Dict, Tuple
from dendropy import Tree, Node
from Bio.Phylo import BaseTree
import math


def annotate_with_closest_centers(tree: Tree, centers: List[str], prior_centers=None, radius=None):
    """
    For each node in the tree finds its closest center and annotates the nodes with the 'center' keyword.
    :param tree: input tree to be annotated
    :param centers: tip labels of the chosen representatives
    :param prior_centers: an optional list of prior centers - closest nodes will be annotated with 'center=-1'.
    :param radius: optional: if provided, only the vertices within the radius are covered by a center.
    """
    all_centers = centers
    if prior_centers:
        all_centers = centers + prior_centers
    closest_centers = find_closest_centers(tree, all_centers)

    for node in tree.nodes():
        closest_center, dist = closest_centers[node]
        if radius and dist > radius:
            continue
        closest_center = all_centers.index(closest_center)
        if closest_center >= len(centers):
            closest_center = -1
        node.annotations.add_new('center', closest_center)


def find_closest_centers(tree: Tree, centers: List[str]) -> Dict[Node, Tuple[str, float]]:
    """
    For each node in the tree finds its closest center.
    A simple O(kn) algorithm where k is the number of centers and n is the number of leaves in the tree.
    :param tree: input tree to be annotated
    :param centers: tip labels of the chosen representatives
    :return: a dict that maps each node to the closest center/leaf and the respective distance on the tree.
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

    center_map = {}
    for node in tree.nodes():
        dist_to_centers = [distances[node] for distances in distances_by_center]
        min_dist = min(dist_to_centers)
        closest_center = dist_to_centers.index(min_dist)
        center_map[node] = (centers[closest_center], min_dist)
    return center_map


def build_distance_functions(tree: BaseTree.Tree, radius=None, prior_centers=None, taxa_weights=None):
    if prior_centers:
        dendropy_tree = Tree.get(data=tree.format(fmt='newick'), schema='newick',
                                 preserve_underscores=True)  # convert to dendropy.
        closest_centers = find_closest_centers(dendropy_tree, prior_centers)
        # For each leaf label stores the distance to the closest prior center:
        closest_prior_dist = dict([(leaf.taxon.label, closest_centers[leaf][1]) for leaf in dendropy_tree.leaf_nodes()])

    distance_functions = {}
    for node in tree.find_clades(order='preorder'):
        if node.is_terminal():
            weight = taxa_weights[node.name] if taxa_weights else 1  # default weight is 1.
            max_dist = closest_prior_dist[node.name] if prior_centers else math.inf
            min_dist = radius if radius else 0
            def dist_func(min_d, max_d, w):
                return lambda dist: (0 if dist <= min_d else (dist if dist <= max_d else max_d)) * w
            if max_dist <= min_dist:
                function = lambda dist: 0
            else:
                function = dist_func(min_dist, max_dist, weight)
            # function = lambda dist: (0 if dist <= min_dist else (dist if dist <= max_dist else max_dist)) * weight
            # if radius:
            #     function = lambda dist: (0 if dist < radius else dist) * weight
            # else:
            #     function = lambda dist: dist * weight
        else:
            function = lambda dist: 0  # non-terminal nodes do not contribute to the objective function.

        distance_functions[node] = function

    return distance_functions
