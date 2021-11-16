# -*- coding: utf-8 -*-
from typing import List

from Bio.Phylo import BaseTree


def find_n_medoids(tree: BaseTree, n: int, tip_clusters: dict, max_dist=None) -> List[str]:
    """
    Finds n medoids on a tree by a modification of Tamir's algorithm for p-median.
    If max_dist is specified, the method finds n representatives that cover
    as much diversity (within max_dist) as possible.

    :param tree: phylogenetic tree in the biopython fromat.
    :param n: number of representatives to be chosen.
    :param tip_clusters: an empty dictionary: on return it will contain for each tip label an index of a representative
                         that 'covers' that tip.
    :param max_dist: an optional parameter that specifies the maximum coverage distance by a single representative.
    :return: a list of tip labels that have been chosen as representatives
    """

    # TODO: Sanket and Sid implementation starts here
