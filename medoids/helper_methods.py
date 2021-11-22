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

def resolveTies(nodelist,treem):#### sort according to condition in paper
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
