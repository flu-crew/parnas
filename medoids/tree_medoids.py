# -*- coding: utf-8 -*-
from io import StringIO
from typing import List
import numpy as np
from Bio import Phylo

from Bio.Phylo import BaseTree

from input import get_Tree_Phylo
from medoidFinder import MedoidFinder


def find_n_medoids(tree: BaseTree, n: int, tip_clusters: dict =None, max_dist=None) -> List[str]:
    """
    Finds n medoids on a tree by a modification of Tamir's algorithm for p-median.
    If max_dist is specified, the method finds n representatives that cover
    as much diversity (within max_dist) as possible.

    :param tree: phylogenetic tree in the biopython fromat.
    :param p: number of representatives to be chosen.
    :param tip_clusters: an empty dictionary: on return it will contain for each tip label an index of a representative
                         that 'covers' that tip.
    :param max_dist: an optional parameter that specifies the maximum coverage distance by a single representative.
    :return: a list of tip labels that have been chosen as representatives
    """
    medoidFinder = MedoidFinder(tree)
    medoids = medoidFinder.find_medoids(n)
    return medoids[1]


if __name__ == "__main__":
   tree = get_Tree_Phylo(input_string="((A:2,B:3):4,(C:5,(D:7,E:1):7):11);")
   print(find_n_medoids(tree,2))
