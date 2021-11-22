# -*- coding: utf-8 -*-
from io import StringIO
from typing import List
import numpy as np
from Bio import Phylo

from Bio.Phylo import BaseTree

from helper_methods import calculate_distance
from helper_methods import sortAcc1
from helper_methods import distance_function
from helper_methods import cost_function
from helper_methods import getTreeNodes

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
    n_c= min(n,tree.count_terminals()) # makes sure that we dont have more centers than leaves
    centers =[]# nitially the list of centers is empty
    nodes = getTreeNodes(tree)#get nodes of tree
    nnodes = len(nodes)# gets number of nodes in tree
    index_lookup = dict([(j, i) for i, j in enumerate(nodes)])
    ## We use distance_lookup to store all distances for a node.
    distance_lookup = dict()
    r_lookup = dict()
    post_order_nodes = list(tree.find_clades(order='postorder'))#goes through the nodes of tree in post order
    for i in post_order_nodes:
        distance_lookup[i] = []#initially the distance list for each node is empty
        r_lookup[i] = []#initially the distance list for each node is empty
        #the loop below will sort the list of distances for each node , where sortAcc1 is the method if the distance of 2 nodes are the same from node i
        for j in post_order_nodes:
            distance_lookup[i].append((j, calculate_distance(i, j, tree)))
            r_lookup[i].append((j, calculate_distance(i, j, tree)))
        distance_lookup[i].sort(key=lambda r: r[1])
        distance_lookup[i] = dict([(i[0], j) for j, i in enumerate(sortAcc1(distance_lookup[i], tree))])
        r_lookup[i] = dict(sortAcc1(r_lookup[i], tree))
        # G (node i, number of centers in subtree rooted at i, r_i^j) where G is the optimal value at subtree at i
        G = np.zeros((n_c + 1, nnodes, nnodes))
        # F (node i, number of centers , r_i^j) where F is the optimal value at for tree
        F = np.zeros((n_c + 1, nnodes, nnodes))
        G_hit = np.full((n_c + 1, nnodes, nnodes), False)
        F_hit = np.full((n_c + 1, nnodes, nnodes), False)
        # fill initial value
        G[0] = np.full((nnodes, nnodes), np.inf)
        G_hit[0] = np.full((nnodes, nnodes), True)
        # fill value for each nodes
        for i in post_order_nodes:
            if i.is_terminal():
                G[1, index_lookup[i]] = np.full((nnodes), 1)
                G_hit[1, index_lookup[i]] = True
                for j in tree.find_clades(lambda x: not i.is_parent_of(x)):
                    F[0, index_lookup[i], distance_lookup[i][j]] = distance_function(calculate_distance(i, j, tree))
                    F_hit[0, index_lookup[i], distance_lookup[i][j]] = True
                    F[1, index_lookup[i], distance_lookup[i][j]] = min(F[0, index_lookup[i], distance_lookup[i][j]],
                                                                       G[1, index_lookup[i], distance_lookup[i][j]])
                    F_hit[1, index_lookup[i], distance_lookup[i][j]] = True
            else:
                for q in range(1, n_c + 1):
                  #  print("q change")
                    for j in distance_lookup[i]:
                   #     print("node change")
                        if i == j:
                            G[q, index_lookup[i], distance_lookup[i][j]] = np.inf
                            G_hit[q, index_lookup[i], distance_lookup[i][j]] = True
                        elif not i.is_parent_of(j):
                            G[q, index_lookup[i], distance_lookup[i][j]] = G[
                                q, index_lookup[i], distance_lookup[i][j] - 1]
                            G_hit[q, index_lookup[i], distance_lookup[i][j]] = True
                        else:
                            left, right = i.clades
                            left_size = len(list(left.find_elements()))
                            right_size = len(list(right.find_elements()))
                            if left.is_parent_of(j):
                                mindist = np.inf
                                for q1 in range(1, left_size + 1):
                                    for q2 in range(0, right_size + 1):
                                        if q1 + q2 == q:
                                            dist = G[q1, index_lookup[left], distance_lookup[left][j]] + F[
                                                q2, index_lookup[right], distance_lookup[right][j]]
                                            if (mindist > dist):
                                                mindist = dist
                                G[q, index_lookup[i], distance_lookup[i][j]] = min(
                                    G[q, index_lookup[i], distance_lookup[i][j] - 1],
                                    distance_function(r_lookup[i][j]) + mindist)
                                G_hit[q, index_lookup[i], distance_lookup[i][j]] = True
                            else:
                                mindist = np.inf
                                for q1 in range(1, right_size + 1):
                                    for q2 in range(0, left_size + 1):
                                        if q1 + q2 == q:
                                            dist = G[q1, index_lookup[right], distance_lookup[right][j]] + F[
                                                q2, index_lookup[left], distance_lookup[left][j]]
                                            if (mindist > dist):
                                                mindist = dist
                                G[q, index_lookup[i], distance_lookup[i][j]] = min(
                                    G[q, index_lookup[i], distance_lookup[i][j] - 1],
                                    distance_function(r_lookup[i][j]) + mindist)
                                G_hit[q, index_lookup[i], distance_lookup[i][j]] = True
                            mindist = np.inf
                            for q1 in range(0, left_size + 1):
                                for q2 in range(0, right_size + 1):
                                    if q1 + q2 == q:
                                        dist = F[q1, index_lookup[left], distance_lookup[left][j]] + F[
                                            q2, index_lookup[right], distance_lookup[right][j]]
                                        if (mindist > dist):
                                            mindist = dist
                            F[q, index_lookup[i], distance_lookup[i][j]] = min(
                                G[q, index_lookup[i], distance_lookup[i][j]],
                                distance_function(r_lookup[i][j]) + mindist)
                            F_hit[q, index_lookup[i], distance_lookup[i][j]] = True

    print("The number of centers is", n)
    print(G)
    print(F)
    print(G_hit)
    print(F_hit)
    print(G[n_c, 0, nnodes - 1])
    return centers


