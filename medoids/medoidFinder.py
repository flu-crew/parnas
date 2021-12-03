# -*- coding: utf-8 -*-
import numpy as np

from helper_methods import cost_function
from helper_methods import distance_function
from helper_methods import calculate_distance
from helper_methods import getTreeNodes


class MedoidFinder(object):
    """
    Class that finds the medoids given a certain tree.
    An adapted implementation of Tamir's algorithm (Tamir 1996) for generalized p-medians.
    """

    def __init__(self,tree):
        self.cost_function = cost_function
        self.distance_function = distance_function

        self.nodes = getTreeNodes(tree)
        self.nnodes = len(self.nodes)
        self.tree = tree

        self.index_lookup = dict([(j,i) for i,j in enumerate(self.nodes)])
        self.distance_lookup = dict()
        self.r_lookup = dict()
        
        self.n_c =  self.tree.count_terminals()
        self.G = np.array([])
        self.F = np.array([])
        self.Gmedian_nodes = np.array([])
        self.Fmedian_nodes = np.array([])

    def find_medoids(self,p):
        post_order_nodes= list(self.tree.find_clades(order='postorder'))
        self.n_c = min(p,self.tree.count_terminals())
        self.G = np.zeros((self.n_c + 1,self.nnodes, self.nnodes))
        self.F = np.zeros((self.n_c + 1,self.nnodes, self.nnodes))
        self.Gmedian_nodes = np.full((self.n_c+1,self.nnodes,self.nnodes),set())
        self.Fmedian_nodes = np.full((self.n_c+1,self.nnodes,self.nnodes),set())
        self.initialize_lookups()
        self.G[0] = np.full((self.nnodes,self.nnodes),np.inf)
        for i in post_order_nodes:
            if i.is_terminal():
                self.initialize_G_and_F(i)
            else:
                for q in range(self.n_c+1):
                    for j in self.distance_lookup[i]:
                        self.computeF(q,i,j)
                        if q > 0:
                            self.computeG(q,i,j)
            
        min_index = np.argmin(self.G[self.n_c,0])
        return self.G[self.n_c,0,min_index],self.Gmedian_nodes[self.n_c,0,min_index]

    def initialize_lookups(self):
        post_order_traversal = list(self.tree.find_clades(order='postorder'))
        for i in post_order_traversal:
            self.distance_lookup[i] = []
            self.r_lookup[i] = []

            for j in self.tree.find_clades(lambda x:i.is_parent_of(x),order='postorder'): ## is_parent_of checks ancestry and not just direct parent
              self.distance_lookup[i].append((j,calculate_distance(i, j, self.tree)))
              self.r_lookup[i].append((j, calculate_distance(i, j, self.tree)))

            for j in self.tree.find_clades(lambda x:not i.is_parent_of(x),order='postorder'):
              self.distance_lookup[i].append((j,calculate_distance(i, j, self.tree)))
              self.r_lookup[i].append((j, calculate_distance(i, j, self.tree)))

            self.distance_lookup[i].sort(key=lambda r: r[1])
            self.r_lookup[i].sort(key=lambda r: r[1])
            self.distance_lookup[i] = dict([(i[0], j) for j, i in enumerate(self.distance_lookup[i])])
            self.r_lookup[i] = dict(self.r_lookup[i])

    def initialize_G_and_F(self,i):
        self.G[1,self.index_lookup[i]] = np.full((self.nnodes),0)
        self.Gmedian_nodes[1,self.index_lookup[i]] = np.full((self.nnodes),set([i]))
        for j in self.tree.find_clades(lambda x: not i.is_parent_of(x)):
            self.F[0, self.index_lookup[i], self.distance_lookup[i][j]] = distance_function(calculate_distance(i, j, self.tree))
            self.Fmedian_nodes[0, self.index_lookup[i], self.distance_lookup[i][j]] = set([j])
            if self.F[0,self.index_lookup[i],self.distance_lookup[i][j]]>self.G[1,self.index_lookup[i],self.distance_lookup[i][j]]:
                self.F[1,self.index_lookup[i],self.distance_lookup[i][j]] = self.G[1,self.index_lookup[i],self.distance_lookup[i][j]]
                self.Fmedian_nodes[1,self.index_lookup[i],self.distance_lookup[i][j]] = self.Gmedian_nodes[1,self.index_lookup[i],self.distance_lookup[i][j]]

            else:
                self.F[1,self.index_lookup[i],self.distance_lookup[i][j]] = self.F[0,self.index_lookup[i],self.distance_lookup[i][j]] 
                self.Fmedian_nodes[1,self.index_lookup[i],self.distance_lookup[i][j]] = self.Fmedian_nodes[1,self.index_lookup[i],self.distance_lookup[i][j]]

    def computeG(self,q,i,j):
        for q in range(1,self.n_c + 1):
            for j in self.distance_lookup[i]:
                if i == j:
                    self.G[q,self.index_lookup[i],self.distance_lookup[i][j]] = np.inf
                    self.Gmedian_nodes[q,self.index_lookup[i],self.distance_lookup[i][j]] = self.Gmedian_nodes[q-1,self.index_lookup[i],self.distance_lookup[i][j]].union(set([i]))
                elif not i.is_parent_of(j):
                    self.G[q,self.index_lookup[i], self.distance_lookup[i][j]] = self.G[q,self.index_lookup[i],self.distance_lookup[i][j]-1]
                    self.Gmedian_nodes[q,self.index_lookup[i],self.distance_lookup[i][j]] = self.Gmedian_nodes[q,self.index_lookup[i],self.distance_lookup[i][j]-1]
                else:
                    left, right = i.clades
                    left_size = len(list(left.find_elements()))
                    right_size = len(list(right.find_elements()))
                    if left.is_parent_of(j):
                        self.computeG_subtree(q,i,j,left,right)
                    else:
                        self.computeG_subtree(q,i,j,right,left)


    def computeG_subtree(self,q,i,j,n1,n2):
        n1_size = len(list(n1.find_elements()))
        n2_size = len(list(n2.find_elements()))
        mindist = np.inf
        min_q1 = 0
        min_q2 = 0
        for q1 in range(max(1,q-n2_size),min(n1_size + 1,q+1)):
            q2 = q - q1
            dist = self.G[q1, self.index_lookup[n1], self.distance_lookup[n1][j]] + \
                    self.F[q2, self.index_lookup[n2], self.distance_lookup[n2][j]]
            if (mindist > dist):
                mindist = dist
                min_q1 = q1
                min_q2 = q2
        node_distance = distance_function(self.r_lookup[i][j]) + mindist
        if(self.G[q, self.index_lookup[i],self.distance_lookup[i][j]-1] > node_distance):
            self.G[q,self.index_lookup[i],self.distance_lookup[i][j]] = node_distance
            medians = self.Gmedian_nodes[min_q1,self.index_lookup[n1],self.distance_lookup[n1][j]]
            medians = medians.union(self.Fmedian_nodes[min_q2,self.index_lookup[n2],self.distance_lookup[n2][j]])
            medians = medians.union(set([j]));
            self.Gmedian_nodes[q,self.index_lookup[i],self.distance_lookup[i][j]] = medians
        else:
            self.G[q,self.index_lookup[i],self.distance_lookup[i][j]] = self.G[q,self.index_lookup[i],self.distance_lookup[i][j]-1]
            self.Gmedian_nodes[q,self.index_lookup[i],self.distance_lookup[i][j]] = self.Gmedian_nodes[q,self.index_lookup[i],self.distance_lookup[i][j]-1]

    def computeF(self,q,i,j):
        left, right = i.clades
        left_size = len(list(left.find_elements()))
        right_size = len(list(right.find_elements()))
        mindist = np.inf
        min_q1 = 0
        min_q2 = 0
        for q1 in range(max(1,q-right_size),min(left_size + 1,q+1)):
            q2 = q - q1
            dist = self.F[q1, self.index_lookup[left], self.distance_lookup[left][j]] + self.F[q2, self.index_lookup[right], self.distance_lookup[right][j]]
            if (mindist > dist):
                mindist = dist
                min_q1 = q1
                min_q2 = q2
        self.F[q, self.index_lookup[i], self.distance_lookup[i][j]] = min(
            self.G[q, self.index_lookup[i], self.distance_lookup[i][j]],
            distance_function(self.r_lookup[i][j]) + mindist)
        node_distance = distance_function(self.r_lookup[i][j]) + mindist
        if(self.G[q,self.index_lookup[i], self.distance_lookup[i][j]] > node_distance):
            self.F[q,self.index_lookup[i],self.distance_lookup[i][j]] =  node_distance
            medians = self.Fmedian_nodes[min_q1, self.index_lookup[left], self.distance_lookup[left][j]]
            medians = medians.union(self.Fmedian_nodes[min_q2, self.index_lookup[right], self.distance_lookup[right][j]])
            medians = medians.union(set([j]))
            self.Fmedian_nodes[q,self.index_lookup[i],self.distance_lookup[i][j]] = medians
        else:
            self.F[q,self.index_lookup[i],self.distance_lookup[i][j]] = self.G[q,self.index_lookup[i],self.distance_lookup[i][j]]
            self.Fmedian_nodes[q,self.index_lookup[i],self.distance_lookup[i][j]] = self.Gmedian_nodes[q,self.index_lookup[i],self.distance_lookup[i][j]]
           
            
