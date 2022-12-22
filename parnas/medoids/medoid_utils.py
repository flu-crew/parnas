# -*- coding: utf-8 -*-
from typing import List, Dict, Tuple
from dendropy import Tree, Node, Taxon
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


class DistFunction(object):
    def __init__(self, is_zero, min_dst=None, max_dist=None, weight=None, is_binary=False):
        if weight is not None and weight < 1e-8:
            is_zero = True
        self.is_zero = is_zero
        self.is_binary = is_binary
        if not is_zero:
            assert None not in (min_dst, max_dist, weight)
            assert min_dst <= max_dist
        self.min_dist = min_dst
        self.max_dist = max_dist
        self.weight = weight

    def get_dist(self, dist: float):
        if self.is_zero:
            return 0

        if dist <= self.min_dist:
            return 0
        if self.min_dist > 0 and self.is_binary:
            return self.weight  # binary means that if a leaf is not covered (dist > min_dist), then contribution is the weight.

        mapped_dist = min(dist, self.max_dist) - self.min_dist  # min_dist (radius) is subtracted from the distance.
        # mapped_dist = dist if (dist <= self.max_dist) else self.max_dist
        mapped_dist *= self.weight
        return mapped_dist


def build_distance_functions(tree: Tree, radius=None, is_binary=False, prior_centers=None, fully_excluded=None,
                             taxa_weights: Dict[str, float]=None) -> Dict[Node, DistFunction]:
    if prior_centers:
        # dendropy_tree = Tree.get(data=tree.format(fmt='newick'), schema='newick',
        #                          preserve_underscores=True)  # convert to dendropy.
        closest_centers = find_closest_centers(tree, prior_centers)
        # For each leaf label stores the distance to the closest prior center:
        closest_prior_dist = dict([(leaf.taxon.label, closest_centers[leaf][1]) for leaf in tree.leaf_nodes()])

    distance_functions = {}
    for node in tree.preorder_node_iter():
        if node.is_leaf():
            excluded = (node.taxon.label in fully_excluded) if fully_excluded else False  # If excluded, dist function is 0.
            weight = taxa_weights.get(node.taxon.label, 1) if taxa_weights else 1  # default weight is 1.
            max_dist = closest_prior_dist[node.taxon.label] if prior_centers else math.inf
            min_dist = radius if radius else 0
            if excluded or max_dist <= min_dist:
                function = DistFunction(is_zero=True)
            else:
                function = DistFunction(False, min_dist, max_dist, weight, is_binary)
        else:
            function = DistFunction(is_zero=True)  # non-terminal nodes do not contribute to the objective function.

        distance_functions[node] = function

    return distance_functions


def get_costs(tree: Tree, excluded=None, fully_excluded=None) -> Dict[str, float]:
    cost_map = {}
    for taxon in tree.taxon_namespace:
        if (excluded and taxon.label in excluded) or (fully_excluded and taxon.label in fully_excluded):
            cost_map[taxon.label] = math.inf
        else:
            cost_map[taxon.label] = 0
    return cost_map


def binarize_tree(tree: Tree, edge_length=0):
    """
    Adds/removes nodes from the tree to make it fully binary (added edges will have length 'edge_length')
    :param tree: Dendropy tree to be made bifurcating.
    """

    # First suppress unifurcations.
    tree.suppress_unifurcations()

    # Now binarize multifurcations.
    for node in tree.postorder_node_iter():
        assert isinstance(node, Node)
        if node.child_nodes() and len(node.child_nodes()) > 2:
            num_children = len(node.child_nodes())
            children = node.child_nodes()
            interim_node = node
            # Creates a caterpillar structure with children on the left of the trunk:
            for child_ind in range(len(children) - 2):
                new_node = Node(edge_length=edge_length)
                interim_node.set_child_nodes([children[child_ind], new_node])
                interim_node = new_node
            interim_node.set_child_nodes(children[num_children - 2:])
