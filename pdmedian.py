from input import get_tree
import numpy as np
import math

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

## TODO Calculate the cost between two nodes
def calculate_distance(u,v):
    """ Calculates the distance between two nodes
        u -- the node from which distance is calculated
        v -- the node from which distance is to be calculated

        :returns: Number
    """
    return 1

if __name__ == "__main__":
   tree=  get_tree(input_string="((A:2,B:3):4,(C:5,(D:7,E:1):7):11);")
   p = min(10,len(tree.leaf_nodes()))
   nodes = tree.nodes()
   nnodes = len(tree.nodes())
   index_lookup = dict([(j,i) for i,j in enumerate(nodes)]) 
   ## We use distance_lookup to store all distances for a node.
   distance_lookup = dict()
   for i in nodes:
       distance_lookup[i] = []
       for j in nodes:
           ## TODO sort this list according to the algorithm given in the book
           distance_lookup[i].append((calculate_distance(i,j),j))
   print(tree.as_ascii_plot())
   G = np.zeros((p,nnodes,nnodes))
   F = np.zeros((p,nnodes,nnodes))
   G[0] = np.full((nnodes,nnodes),np.inf)
   G[1] = np.full((nnodes,nnodes),cost_function())
   for i in nodes:
       for j in tree.find_nodes(lambda x:x not in i.child_nodes()):
            F[0,index_lookup[i],index_lookup[j]] = distance_function(calculate_distance(i,j))
   ## Computation for G
   for q in range(2,p): 
        for v in tree.postorder_iter():
            for (i,r) in enumerate(distance_lookup[v][1:]):
                if r[1] not in v.child_nodes(): ## aka v_i^j \in V - V_j
                    G[q,index_lookup[v],i] = G[q,index_lookup[v],i-1]
                else: 
                    ## TODO Write a method that differentiates between right and left children of a tree
                    child = v.child_node_iter(lambda x:r[1] in x.child_nodes())[0]
                    min_child_distance = math.inf
                    r_child = 0
                    i_child = -1
                    # To find r^k_j(1)
                    for (i_c,r_c) in enumerate(distance_lookup[child]):
                        if r_c[1] == r[1]:
                            r_child = r_c
                            i_child = i_c
                    ## To compute min(G,F)
                    for q1 in range(q):
                        for q2 in range(q):
                            if (q1+q2)==q:
                                ## !! The third parameter for F should be the other child not the same one
                                child_distance = G[q1,index_lookup[child],i_child] + F[q1,index_lookup[child],i_child]
                                if child_distance<min_child_distance:
                                    min_child_distance = child_distance
                    G[q,index_lookup[v],i] = min(G[q,index_lookup[v],i-1],distance_function(r[0])+min_child_distance)
                    


   print(G)
   print(F)
