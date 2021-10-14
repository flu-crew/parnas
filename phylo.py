from io import StringIO
import phylotreelib
from Bio import Phylo
import numpy as np
import math

tree = Phylo.read(StringIO("(A:0.1,B:0.2,E(C:0.3,D:0.4):0.5)F"), "newick");

inner_nodes = tree.get_nonterminals()# gets internal nodes
leaf_nodes = tree.get_terminals()# gets leaf nodes


def getAllNodes(): #simple function to get all nodes in a list
    all_nodes = []
    for x in inner_nodes:
        all_nodes.append(x)
    for x in leaf_nodes:
        all_nodes.append(x)

    return all_nodes

all_nodes = getAllNodes()

def setConstantCostDict(value):# create a dictionary for costs and give all nodes a cost of 'value'
    res = {}
    for key in all_nodes:
        res[key.name] = value
    return res

currentcosts = setConstantCostDict(1)


def getNodeCost(node):
   return currentcosts[node]



def getNodeDistList(node):
    retlist =[]
    for x in all_nodes:
        retlist.append((x.name+node,tree.distance(node,x)))
    retlist.sort(key=lambda x: x[1])
    return retlist

#print(getNodeDistList("B"))


def G(node,n,closest_distance):
    ##base cases
    if (n == 0):
     return np.inf
    elif (node.is_terminal() and n==1):
     return getNodeCost(node)
    ##need to do other cases


def F(node,n,closest_distance):
    ##base cases
    if (n == 0):
     return closest_distance
    elif ( n==1):
     return min(F(node,0,closest_distance), G(node,1,closest_distance))
    ##need to do other cases


#jg;klfgj;oenk







