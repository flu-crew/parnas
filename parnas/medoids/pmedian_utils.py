# -*- coding: utf-8 -*-
from typing import Callable, Iterable, Optional, List, Tuple
from dendropy import Tree, Node


def dfs_tree_traversal(node: Node) -> List[Tuple[Node, float]]:
    node_dist_pairs = []
    node_stack = []  # node with prev_node and distance.
    node_stack.append((node, None, 0))
    while node_stack:
        node, prev_node, cur_dist = node_stack.pop()
        node_dist_pairs.append((node, cur_dist))
        neighbors = []
        if node.child_nodes():
            neighbors += node.child_nodes()
        if node.parent_node:
            neighbors.append(node.parent_node)
        for neighbor in neighbors:
            if neighbor is not prev_node:
                is_parent = node.parent_node is neighbor
                edge_len = node.edge_length if is_parent else neighbor.edge_length
                node_stack.append((neighbor, node, cur_dist + edge_len))
    return node_dist_pairs


def filtered_preorder_iterator(tree: Tree, subtree_filter: Optional[Callable[[Node], bool]]) -> Iterable[Node]:
    """
    Preorder iterator that will avoid subtrees rooted at nodes, s.t. subtree_filter(node) is False.
    """
    stack = [tree.seed_node]
    while stack:
        node = stack.pop()
        if subtree_filter is None or subtree_filter(node):
            yield node
            stack.extend(n for n in reversed(node._child_nodes))


def filtered_postorder_iterator(tree: Tree, subtree_filter: Optional[Callable[[Node], bool]]) -> Iterable[Node]:
    stack = [(tree.seed_node, False)]
    while stack:
        node, state = stack.pop()
        if state:
            if subtree_filter is None or subtree_filter(node):
                yield node
        elif subtree_filter is None or subtree_filter(node):
            stack.append((node, True))
            stack.extend([(n, False) for n in reversed(node._child_nodes)])


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


# def resolve_ties(nodelist, treem):  # sort according to condition in paper
#     for i in range(len(nodelist)):
#         if nodelist[i - 1][1] == nodelist[i][1]:  # checks for equidistant vertices vk and vm where vk appears before vm
#             if treem.common_ancestor(nodelist[i - 1][0], nodelist[0][0]) == nodelist[0][0]:  # check if vk is in subtree induced by vj
#                 if treem.common_ancestor(nodelist[i][0], nodelist[0][0]) == nodelist[0][0]:  # check if vm is in subtree induced by vj
#                     if treem.common_ancestor(nodelist[i - 1][0], nodelist[i][0]) == nodelist[i][0]:  ##IF vm is the parent of vk then switch postion of vk and vm
#                         temp = nodelist[i]
#                         nodelist[i] = nodelist[i - 1]
#                         nodelist[i - 1] = temp
#
#             if treem.common_ancestor(nodelist[i][0], nodelist[0][0]) == nodelist[0][0]:  # check if vm is in subtree induced by vj
#                 if treem.common_ancestor(nodelist[i - 1][0], nodelist[0][0]) != nodelist[0][0]:  # if vk is in not in induced by vj
#                     # then switch postion of vk and vm
#                     temp = nodelist[i]
#                     nodelist[i] = nodelist[i - 1]
#                     nodelist[i - 1] = temp
#
#     return nodelist


def calculate_distance(u, v, tree):
    """ Calculates the distance between two nodes
        u -- the node from which distance is calculated
        v -- the node from which distance is to be calculated
        :returns: Number 
    """
    return tree.distance(u, v)
