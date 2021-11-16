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
    tree=  get_Tree_Phylo(input_string="((F:2,E:1):6,D:7)")
    p = min(2,tree.count_terminals())
    nodes = getTreeNodes(tree) 
    nnodes = len(nodes)
    index_lookup = dict([(j,i) for i,j in enumerate(nodes)]) 
    ## We use distance_lookup to store all distances for a node.
    distance_lookup = dict()
    r_lookup = dict()
    post_order_nodes = list(tree.find_clades(order='postorder'))
    for i in post_order_nodes:
        distance_lookup[i] = []
        r_lookup[i] = []
        for j in post_order_nodes:
            ## TODO sort this list according to the algorithm given in the book
            distance_lookup[i].append((j,calculate_distance(i,j,tree)))
            r_lookup[i].append((j,calculate_distance(i,j,tree)))
        distance_lookup[i].sort(key=lambda r:r[1])
        distance_lookup[i] = dict([(i[0],j) for j,i in enumerate(sortAcc1(distance_lookup[i],tree))])
        r_lookup[i] = dict(sortAcc1(r_lookup[i],tree))
    G = np.zeros((p+1,nnodes,nnodes))
    F = np.zeros((p+1,nnodes,nnodes))
    G_hit = np.full((p+1,nnodes,nnodes),False)
    F_hit = np.full((p+1,nnodes,nnodes),False)
    G[0] = np.full((nnodes,nnodes),np.inf)
    G_hit[0] = np.full((nnodes,nnodes),True)
    for i in post_order_nodes:
        if i.is_terminal():
            G[1,index_lookup[i]] = np.full((nnodes),1)
            G_hit[1,index_lookup[i]] = True
            for j in tree.find_clades(lambda x:not i.is_parent_of(x)):
                F[0,index_lookup[i],distance_lookup[i][j]] = distance_function(calculate_distance(i,j,tree))
                F_hit[0,index_lookup[i],distance_lookup[i][j]] = True
                F[1,index_lookup[i],distance_lookup[i][j]] = min(F[0,index_lookup[i],distance_lookup[i][j]],G[1,index_lookup[i],distance_lookup[i][j]])
                F_hit[1,index_lookup[i],distance_lookup[i][j]] = True
        else:
            for q in range(1,p+1):
                print("q change")
                for j in distance_lookup[i]:
                    print("node change")
                    if i==j:
                        G[q,index_lookup[i],distance_lookup[i][j]] = np.inf
                        G_hit[q,index_lookup[i],distance_lookup[i][j]]=True 
                    elif not i.is_parent_of(j):
                        G[q,index_lookup[i],distance_lookup[i][j]] = G[q,index_lookup[i],distance_lookup[i][j]-1]
                        G_hit[q,index_lookup[i],distance_lookup[i][j]]=True
                    else:
                        left,right = i.clades
                        left_size = len(list(left.find_elements()))
                        right_size = len(list(right.find_elements()))
                        if left.is_parent_of(j):
                            mindist = np.inf
                            for q1 in range(1,left_size+1):
                                for q2 in range(0,right_size+1):
                                    if q1+q2==q:
                                        dist = G[q1,index_lookup[left],distance_lookup[left][j]] + F[q2,index_lookup[right],distance_lookup[right][j]]
                                        if(mindist>dist):
                                            mindist = dist
                            G[q,index_lookup[i],distance_lookup[i][j]] = min(G[q,index_lookup[i],distance_lookup[i][j]-1],distance_function(r_lookup[i][j])+mindist)
                            G_hit[q,index_lookup[i],distance_lookup[i][j]]=True
                        else:
                            mindist = np.inf
                            for q1 in range(1,right_size+1):
                                for q2 in range(0,left_size+1):
                                    if q1+q2==q:
                                        dist = G[q1,index_lookup[right],distance_lookup[right][j]] + F[q2,index_lookup[left],distance_lookup[left][j]]
                                        if(mindist>dist):
                                            mindist = dist
                            G[q,index_lookup[i],distance_lookup[i][j]] = min(G[q,index_lookup[i],distance_lookup[i][j]-1],distance_function(r_lookup[i][j])+mindist)
                            G_hit[q,index_lookup[i],distance_lookup[i][j]]=True
                        mindist = np.inf
                        for q1 in range(0,left_size+1):
                            for q2 in range(0,right_size+1):
                                if q1+q2==q:
                                    dist = F[q1,index_lookup[left],distance_lookup[left][j]] + F[q2,index_lookup[right],distance_lookup[right][j]]
                                    if(mindist>dist):
                                        mindist = dist
                        F[q,index_lookup[i],distance_lookup[i][j]] = min(G[q,index_lookup[i],distance_lookup[i][j]],distance_function(r_lookup[i][j]) + mindist)
                        F_hit[q,index_lookup[i],distance_lookup[i][j]]=True
    print(G)
    print(F)
    print(G_hit)
    print(F_hit)
    print(G[p,0,nnodes-1])

