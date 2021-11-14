from input import get_Tree_Phylo
import numpy as np
import math

def getTreeNodes(tree):
    nodes = tree.get_nonterminals()
    nodes.extend(tree.get_terminals())
    return nodes

def cost_function():
    """ Returns the setup cost c_j  for a node
    :returns: Number 

    """
    return 0

def distance_function(distance):
    """ Returns the value of distance function f_j
        distance -- the distance between two nodes

    :returns: Number

    """
    return distance

def sortNodeList(nodelist,treem):#### sort according to condition in paper
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

## TODO Calculate the cost between two nodes
def calculate_distance(u,v,tree):
    """ Calculates the distance between two nodes
        u -- the node from which distance is calculated
        v -- the node from which distance is to be calculated

        :returns: Number 
    """
    return tree.distance(u,v)

if __name__ == "__main__":
    tree=  get_Tree_Phylo(input_string="((A:2,B:3):4,(C:5,(D:7,E:1):7):11);")
    p = min(10,tree.count_terminals())
    nodes = getTreeNodes(tree) 
    nnodes = len(nodes)
    index_lookup = dict([(j,i) for i,j in enumerate(nodes)]) 
    ## We use distance_lookup to store all distances for a node.
    distance_lookup = dict()
    r_lookup = dict()
    for i in nodes:
        distance_lookup[i] = []
        r_lookup[i] = []
        for j in nodes:
            ## TODO sort this list according to the algorithm given in the book
            distance_lookup[i].append((j,calculate_distance(i,j,tree)))
            r_lookup[i].append((j,calculate_distance(i,j,tree)))
        distance_lookup[i] = dict([(i,j) for j,i in enumerate(sortAcc1(distance_lookup[i],tree))])
        r_lookup[i] = dict(sortAcc1(r_lookup[i],tree))
    G = np.zeros((p,nnodes,nnodes))
    F = np.zeros((p,nnodes,nnodes))
    G[0] = np.full((nnodes,nnodes),np.inf)
    post_order_nodes = list(tree.find_elements(order='postorder'))[:-1]
    for i in post_order_nodes:
        if i.is_terminal():
            G[1,index_lookup[i]] = np.full((nnodes),1)
        for j in tree.find_elements(lambda x:not i.is_parent_of(x)):
            F[0,index_lookup[i],distance_lookup[i][j]] = distance_function(calculate_distance(i,j,tree))
            F[1,index_lookup[i],distance_lookup[i][j]] = min(F[0,index_lookup[i],distance_lookup[i][j]],G[1,index_lookup[i],distance_lookup[i][j]])
        else:
            for q in range(1,p):
                for j in post_order_nodes:
                    if i==j:
                        G[q,index_lookup[i],distance_lookup[i][j]] = np.inf
                    elif not i.is_parent_of(j):
                        G[q,index_lookup[i],distance_lookup[i][j]] = G[q,index_lookup[i],distance_lookup[i][j]-1]
                    else:
                        left,right = i.clades
                        left_size = len(list(left.find_elements()))
                        right_size = len(list(right.find_elements()))
                        if left.is_parent_of(j):
                            mindist = np.inf
                            for q1 in range(min(1,q-right_size),left_size):
                                q2 = q - q1
                                dist = G[q1,index_lookup[left],distance_lookup[left][j]] + F[q2,index_lookup[right],distance_lookup[right][j]]
                                if(mindist>dist):
                                    mindist = dist
                            G[q,index_lookup[i],distance_lookup[i][j]] = min(G[q,index_lookup[i],distance_lookup[i][j]-1],distance_function(r_lookup[i][j])+mindist)
                        else:
                            mindist = np.inf
                            for q1 in range(min(1,q-left_size),right_size):
                                q2 = q - q1
                                dist = G[q1,index_lookup[right],distance_lookup[right][j]] + F[q2,index_lookup[left],distance_lookup[left][j]]
                                if(mindist>dist):
                                    mindist = dist
                            G[q,index_lookup[i],distance_lookup[i][j]] = min(G[q,index_lookup[i],distance_lookup[i][j]-1],distance_function(r_lookup[i][j])+mindist)
                        mindist = np.inf
                        for q1 in range(max(0,q-right_size),left_size):
                            q2 = q -q1
                            dist = F[q1,index_lookup[left],distance_lookup[left][j]] + F[q2,index_lookup[right],distance_lookup[right][j]]
                            if(mindist>dist):
                                mindist = dist
                        F[q,index_lookup[i],distance_lookup[i][j]] = min(G[q,index_lookup[i],distance_lookup[i][j]],distance_function(r_lookup[i][j]) + mindist)
    print(G)


