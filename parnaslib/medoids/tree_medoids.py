# -*- coding: utf-8 -*-
from typing import List, Tuple, Dict
from dendropy import Tree

from .fast_pmedian_finder import FastPMedianFinder
from .pmedian_finder import PMedianFinder


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


# if __name__ == "__main__":
#    tree = get_Tree_Phylo(input_string="((A:2,B:3):4,(C:5,(D:7,E:1):7):11);")
#    #tree = get_Tree_Phylo(input_string="((A:23,B:27):47,(C:35,(D:76,E:18):28):31);")
#
#    print(find_n_medoids(tree, 3, None))
