# -*- coding: utf-8 -*-
from io import StringIO
from typing import List
import numpy as np
from Bio import Phylo
from typing import List, Tuple, Dict

from Bio.Phylo import BaseTree

from input import get_Tree_Phylo
from medoidFinder import MedoidFinder

def find_n_medoids(tree: BaseTree, n: int, distance_functions: Dict, max_dist=None) -> Tuple[List[str], float]:
    """
    Finds n medoids on a tree by a modification of Tamir's algorithm for p-median.
    If max_dist is specified, the method finds n representatives that cover
    as much diversity (within max_dist) as possible.

    :param tree: phylogenetic tree in the biopython fromat.
    :param n: number of representatives to be chosen.
    :param max_dist: an optional parameter that specifies the maximum coverage distance by a single representative.
    :param distance_functions: a map that links a node (Phylo.Clade) to the distance function of that node.
    :return: (1) a list of tip labels that have been chosen as representatives;
             (2) the minimal objective function value.
    """
    medoidFinder = MedoidFinder(tree)
    medoids = medoidFinder.find_medoids(n)
    return medoids[1]


if __name__ == "__main__":
   tree = get_Tree_Phylo(input_string="((A:2,B:3):4,(C:5,(D:7,E:1):7):11);")
   #tree = get_Tree_Phylo(input_string="((A:23,B:27):47,(C:35,(D:76,E:18):28):31);")

   print(find_n_medoids(tree,3,None))
