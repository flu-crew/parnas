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
#print(checktree)


def getAllNodes(treem): #simple function to get all nodes in a list
    inner_nodes = treem.get_nonterminals()  # gets internal nodes
    leaf_nodes = treem.get_terminals()  # gets leaf nodes
    all_nodes = []
    for x in inner_nodes:
        all_nodes.append(x)
    for x in leaf_nodes:
        all_nodes.append(x)
    return all_nodes


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

def sortAcc1(list,treem):#### sort according to condition in paper
    for i in range(len(list)):
        if list[i - 1][1] == list[i][1]:  ##checks for equidistant vertices vk and vm where vk appears before vm
            if treem.common_ancestor(list[i - 1][0],list[0][0]) == list[0][0]:###check if vk is in subtree induced by vj
                if treem.common_ancestor(list[i][0],list[0][0]) == list[0][0]:###check if vm is in subtree induced by vj
                    if treem.common_ancestor(list[i-1][0],list[i][0]) == list[i][0]:##IF vm is the parent of vk
                        # then switch postion of vk and vm
                      temp = list[i]
                      list[i] = list[i - 1]
                      list[i - 1] = temp

        if treem.common_ancestor(list[i][0], list[0][0]) == list[0][0]:  ###check if vm is in subtree induced by vj
            if treem.common_ancestor(list[i-1][0], list[0][0]) != list[0][0]:  ### if vk is in not in induced by vj
                # then switch postion of vk and vm
                    temp = list[i]
                    list[i] = list[i - 1]
                    list[i - 1] = temp

    return list


Dl= getNodeDict("C",checktree)
print(Dl)
"""
def sortAcc(list,treem):### sorts list according to the paper
    for i in range(len(list)):
        if list[i - 1][1] == list[i][1]:##checks for equidistant vertices vk and vm where vk appears before vm
          treem1 =treem.from_clade(list[0][0])#gets the subtree rooted at the vertex rj
         # print("cehcking if", list[i - 1][0] ,"is in tree")
         # print(treem1)
         # print(treem1.is_parent_of("E"))
         # print(treem1.is_parent_of("A"))
          if treem1.is_parent_of(list[i - 1][0]):#### check if vk appears in subtree at rj

              if treem1.is_parent_of(list[i][0]):### check if vm appears in subtree at rj and if it does then
                  pathl = treem.trace(list[0][0], list[i - 1][0])####  get vertices in path from rj to vk
                  checkl = []
                  del pathl[0]
                  del pathl[len(pathl) - 1]
                  for x in pathl:
                      checkl.append(x.name)
                      ##check for parent constraint
                      ##### if vm is in path to vk , it means vm precedes vk thus vm should be placed earlier in the list
                      if list[i][0] in checkl:
                          temp = list[i]
                          list[i] = list[i - 1]
                          list[i - 1] = temp

          if treem1.is_parent_of(list[i][0]):####if vm appears in subtree rooted at rj
              if treem1.is_parent_of(list[i-1][0]) == False:####if vk does not appear in subtree of rj
                  ### then vm should appear before vk
                  temp = list[i]
                  list[i] = list[i - 1]
                  list[i - 1] = temp

    return list
"""

"""
bl= getNodeDict("E",checktree)
for i in range(len(bl)):
    print(bl[i])
"""









"""
      ### now check if equidistant vertex belongs to
          treem1_all= getAllNodes(treem1)
          pathl = treem.trace(list[0][0], list[i - 1][0])  ####get vertices in path from original node to the node that has equal dist
        ### will turn this in mthod that takes nide as well
        checkl = []
        del pathl[0]
        del pathl[len(pathl) - 1]
         ### check if node vi is in the path to the node vk , where both vi and vk are eqidistant from node r
        for x in pathl:
            checkl.append(x.name)
            ##check for parent constraint
            ##### if vi is in path to vk , it means vi precedes vk thus v1 should be placed earlier in the list
            if list[i][0] in checkl:
                temp = list[i][0]
                list[i][0] = list[i - 1][0]
                list[i - 1][0] = temp
"""





















"""
Dl= getNodeDict("E",checktree)
print(Dl)
### checking the constraints
for i in range(len(Dl)):
    if Dl[i-1][1] == Dl[i][1]:
        pathl = treem.trace("E", Dl[i - 1][0])  ####get vertices in path from original node to the node that has equal dist
        ### will turn this in mthod that takes nide as well
        checkl = []
        del pathl[0]
        del pathl[len(pathl) - 1]
       ### check if node vi is in the path to the node vk , where both vi and vk are eqidistant from node r
        for x in pathl:
            checkl.append(x.name)
##### if vi is in path to vk , it means vi precedes vk thus v1 should be placed earlier in the list
        if Dl[i][0] in checkl:
            temp = Dl[i][0]
            Dl[i][0] = Dl[i - 1][0]
            Dl[i - 1][0] = temp

           #print ("checking constraints as",Dl[i-1][0],"equals", Dl[i][0] )

"""





















"""
(ii) If vk and vm are two distinct nodes in Vj, and

Vm is a descendant of vk, then the element of Lj rep-
resenting vk will precede the one representing v,,. In

particular, vj corresponds to r).
"""



"""
def G(node,n,closest_distance):
    ##base cases
    if (n == 0):
     return np.inf
    elif (node.is_terminal() and n==1):
     return getNodeCost(node)
    ##need to do other cases
"""



def F(node,n,closest_distance):
    ##base cases
    if (n == 0):
     return closest_distance
    elif ( n==1):
     return min(F(node,0,closest_distance), G(node,1,closest_distance))
    ##need to do other cases


#jg;klfgj;oenk

##################################################################################################
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









































