# -*- coding: utf-8 -*-
from typing import List, Tuple, Dict, Optional
from dendropy import Tree

from .fast_pmedian_finder import FastPMedianFinder
from .pmedian_finder import PMedianFinder
from .medoid_utils import find_closest_centers
from .tree_coverage import TreeCoverage
from ..logging import parnas_logger


def find_n_medoids(tree: Tree, n: int, distance_functions: Dict, cost_map: Dict[str, float], max_dist=None)\
        -> Tuple[List[str], float]:
    """
    Finds n medoids on a tree by a modification of Tamir's algorithm for p-median.
    If max_dist is specified, the method finds n representatives that cover
    as much diversity (within max_dist) as possible.

    :param tree: phylogenetic tree in the dendropy format.
    :param n: number of representatives to be chosen.
    :param distance_functions: a map that links a node (dendropy.Node) to the distance function of that node.
    :param max_dist: an optional parameter that specifies the maximum coverage distance by a single representative.
    :return: (1) a list of tip labels that have been chosen as representatives;
             (2) the minimal objective function value.
    """
    medoidFinder = FastPMedianFinder(tree)
    objective, medoids = medoidFinder.find_medoids(n, distance_functions, cost_map)

    return medoids, objective


def find_n_medoids_with_diversity(tree: Tree, n: int, distance_functions: Dict, cost_map: Dict[str, float], max_dist=None)\
        -> Tuple[List[str], float, List[float]]:
    """
    Mirrors "find_n_medoids", but adds a diversity score to the output.
    """
    medoidFinder = FastPMedianFinder(tree)
    objective, medoids = medoidFinder.find_medoids(n, distance_functions, cost_map)
    diversity_scores = []
    if n > 1:
        obj1 = medoidFinder.get_score(1)  # Objective function value for n=1 (1 medoid).
        if obj1 > 0:
            for k in range(2, n + 1):
                objk = medoidFinder.get_score(k)
                diversity_score = (obj1 - objk) / obj1 * 100  # Calculate the % of covered diversity.
                diversity_scores.append(diversity_score)
    return medoids, objective, diversity_scores


def find_coverage(tree: Tree, radius: float, cost_map: Dict[str, float], prior_centers: List[str] = None,
                  fully_excluded: List[str] = None, obj_excluded: List[str] = None) -> Optional[List[str]]:
    """
    Find an optimal set of taxa that fully cover all leaves within the specified radius.
    """
    tree_copy = tree.clone()
    labels_to_prune = []
    if fully_excluded:
        labels_to_prune.extend(fully_excluded)
    # Find leaves that are already covered by prior centers:
    covered_leaves = set()
    if prior_centers:
        closest_priors = find_closest_centers(tree, prior_centers)
        for node, (prior_center, dist) in closest_priors.items():
            if dist <= radius and node.taxon:
                # labels_to_prune.append(node.taxon.label)  # Can't prune covered leaves -- they can be used as reps.
                covered_leaves.add(node.taxon.label)
    if obj_excluded:  # Also ignore taxa that are excluded from the objective
        for label in obj_excluded:
            covered_leaves.add(label)
    # Prune all fully excluded leaves:
    if labels_to_prune:
        if len(labels_to_prune) == len(tree_copy.leaf_nodes()):
            return []
        else:
            tree_copy.prune_taxa_with_labels(labels_to_prune)
            # parnas_logger.debug("After pruning:" + tree_copy.as_string("newick"))
    if covered_leaves and len(covered_leaves) == len(tree_copy.leaf_nodes()):
        return []

    tree_coverer = TreeCoverage(tree_copy)
    coverage = tree_coverer.find_coverage(radius, cost_map, covered_leaves)  # Compute the coverage.
    return coverage


# if __name__ == "__main__":
#    tree = get_Tree_Phylo(input_string="((A:2,B:3):4,(C:5,(D:7,E:1):7):11);")
#    #tree = get_Tree_Phylo(input_string="((A:23,B:27):47,(C:35,(D:76,E:18):28):31);")
#
#    print(find_n_medoids(tree, 3, None))
