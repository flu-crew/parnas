from io import StringIO
#import phylotreelib
from Bio import Phylo
import numpy as np
import math

"""
1) method to get all nodes---getAllNodes(input:tree)--- output list of all nodes
2) set a constant cost to all nodes --setConstantCostDict(value to set cost, tree)-- no output just sets the cost function
3) method to get a list of tupples with sorted distances from a node --getNodeDict(node,treem) , this takes a node ad returns a sorted list if tupples with (vertex,distance fromnode to vertex)
4)sortAcc1--- method called in get NodeDict to sort according to paper

"""


checktree = Phylo.read(StringIO("ROOT(A(A1:0.5,A2:0.456):0.1,B:0.2,E(C:0.3,D(F:0.7,G:0.643):0.4):0.5)"), "newick");

print (checktree)


def getAllNodes(treem): #simple function to get all nodes in a list
    inner_nodes = treem.get_nonterminals()  # gets internal nodes
    leaf_nodes = treem.get_terminals()  # gets leaf nodes
    all_nodes = []
    for x in inner_nodes:
        all_nodes.append(x)
    for x in leaf_nodes:
        all_nodes.append(x)
    return all_nodes

x = getAllNodes(checktree)
print (x)
print(x[4].is_parent_of(x[0]))



def setConstantCostDict(value,treem):# create a dictionary for costs and give all nodes a cost of 'value'
    res = {}
    for key in getAllNodes(treem) :
        res[key.name] = value
    return res

def getNodeDict(node,treem):
    retlist =[]
    for x in getAllNodes(treem):
      retlist.append((x.name,treem.distance(node,x)))###adds tupple of form( vertex ,distance from node to vertex)
    retlist.sort(key = lambda r:r[1])### returns a numerically sorted list of tupples
   # retlist = sortAcc(retlist,treem) ### sends list to sort according to the paper
    retlist = sortAcc1(retlist, treem) ### sends list to sort according to the paper
    return retlist

def sortAcc1(nodelist,treem):#### sort according to condition in paper
    for i in range(len(nodelist)):
        if nodelist[i - 1][1] == nodelist[i][1]:  ##checks for equidistant vertices vk and vm where vk appears before vm
            if treem.common_ancestor(nodelist[i - 1][0],nodelist[0][0]) == nodelist[0][0]:###check if vk is in subtree induced by vj
                if treem.common_ancestor(nodelist[i][0],nodelist[0][0]) == nodelist[0][0]:###check if vm is in subtree induced by vj
                    if treem.common_ancestor(nodelist[i-1][0],nodelist[i][0]) == nodelist[i][0]:##IF vm is the parent of vk
                        # then switch postion of vk and vm
                      temp = nodelist[i]
                      nodelist[i] = nodelist[i - 1]
                      nodelist[i - 1] = temp

        if treem.common_ancestor(nodelist[i][0], nodelist[0][0]) == nodelist[0][0]:  ###check if vm is in subtree induced by vj
            if treem.common_ancestor(nodelist[i-1][0], nodelist[0][0]) != nodelist[0][0]:  ### if vk is in not in induced by vj
                # then switch postion of vk and vm
                    temp = nodelist[i]
                    nodelist[i] = nodelist[i - 1]
                    nodelist[i - 1] = temp

    return nodelist


Dl= getNodeDict("C",checktree)
print(Dl)

def F(node,n,closest_distance):
    ##base cases
    if (n == 0):
     return closest_distance
    elif ( n==1):
     return min(F(node,0,closest_distance), G(node,1,closest_distance))
    ##need to do other cases

"""
if __name__ == "__main__":
   # tree = get_tree(input_string="((A:2,B:3):4,(C:5,(D:7,E:1):7):11);")
    tree = Phylo.read(StringIO("(A:2,B:3):4,(C:5,(D:7,E:1):7):11)"), "newick");
    p = min(10, len(tree.get_terminals()))
    all_nodes = []
    for x in tree.get_nonterminals:
     all_nodes.append(x)
    for x in tree.get_terminals:
     all_nodes.append(x)
    nnodes = len(all_nodes)
    index_lookup = dict([(j, i) for i, j in enumerate(all_nodes)])
     ## We use distance_lookup to store all distances for a node.
   distance_lookup = dict()
    for i in all_nodes:
     distance_lookup[i] = []
    for j in all_nodes:
        distance_lookup[i].append((j,tree.distance(i, j)))#### adds the distance of that node to all other nodes

    distance_lookup[i].sort(key=lambda x: x[i])####sorts the distance in ascending order

    G = np.zeros((p, nnodes, nnodes))
    F = np.zeros((p, nnodes, nnodes))
    G[0] = np.full((nnodes, nnodes), np.inf)
    G[1] = np.full((nnodes, nnodes), 1)

"""
