# import numpy as np
# from dendropy import Tree
# from Bio import Phylo
# from io import StringIO
#
# def get_tree(filepath=None,input_string=None):
#     ## TODO add parameters for different input schemas
#     """ Returns a tree based on the input
#
#     Keyword argument:
#     filepath -- the path to the tree file (default None)
#     string -- a string containing the tree
#
#
#     :returns: Dendropy.Tree
#
#     """
#     if filepath:
#         return Tree.get(path=filepath,schema="newick")
#     elif input_string:
#         return Tree.get(data=input_string,schema="newick")
#     else:
#         return Tree()
#
#
# def get_Tree_Phylo(filepath=None,input_string=None):
#     if filepath:
#         return Phylo.parse(filepath,"newick")
#     elif input_string:
#         return Phylo.read(StringIO(input_string),"newick")
#     else:
#         return Tree()
#
# if __name__ == "__main__":
#     t1 = get_tree(input_string="((A,B),(C,(D,E)));") ## TODO add a test input file
#     print(t1.as_ascii_plot())
