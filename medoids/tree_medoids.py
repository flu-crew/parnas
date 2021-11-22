# -*- coding: utf-8 -*-
from typing import List, Tuple

from Bio.Phylo import BaseTree


def find_n_medoids(tree: BaseTree, n: int, max_dist=None) -> Tuple[List[str], float]:
    """
    Finds n medoids on a tree by a modification of Tamir's algorithm for p-median.
    If max_dist is specified, the method finds n representatives that cover
    as much diversity (within max_dist) as possible.

    :param tree: phylogenetic tree in the biopython fromat.
    :param n: number of representatives to be chosen.
    :param max_dist: an optional parameter that specifies the maximum coverage distance by a single representative.
    :return: (1) a list of tip labels that have been chosen as representatives;
             (2) the minimal objective function value.
    """

    # TODO: Sanket and Sid implementation starts here
