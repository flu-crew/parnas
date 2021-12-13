# -*- coding: utf-8 -*-

import unittest
from dendropy import Tree

from medoids.medoid_utils import binarize_tree
from medoids.pmedian_utils import filtered_preorder_iterator, filtered_postorder_iterator

class TestMedoidUtils(unittest.TestCase):

    def test_binarize_tree(self):
        tree = Tree.get(data='((a:3,b:3,c:5,d:5):2,((e:2):1):2);', schema='newick')
        expected_result = '((a:3.0,(b:3.0,(c:5.0,d:5.0):0.0):0.0):2.0,e:5.0)'
        binarize_tree(tree)
        self.assertEqual(expected_result, str(tree))

    def test_filtered_preorder(self):
        tree = Tree.get(data='(((a,b)n1,(c,d)n2)n3,((e,f)n4,g)n5)n6;', schema='newick')
        n1 = tree.nodes(filter_fn=lambda v: v.label == 'n1')[0]
        preorder = ''
        for node in filtered_preorder_iterator(tree, subtree_filter=lambda v: v is not n1):
            if node.is_leaf():
                preorder += node.taxon.label
            else:
                preorder += node.label
        expected_preorder = 'n6n3n2cdn5n4efg'
        self.assertEqual(expected_preorder, preorder)

    def test_filtered_postorder(self):
        tree = Tree.get(data='(((a,b)n1,(c,d)n2)n3,((e,f)n4,g)n5)n6;', schema='newick')
        n1 = tree.nodes(filter_fn=lambda v: v.label == 'n1')[0]
        postorder = ''
        for node in filtered_postorder_iterator(tree, subtree_filter=lambda v: v is not n1):
            if node.is_leaf():
                postorder += node.taxon.label
            else:
                postorder += node.label
        expected_preorder = 'cdn2n3efn4gn5n6'
        self.assertEqual(expected_preorder, postorder)


if __name__ == '__main__':
    unittest.main()
