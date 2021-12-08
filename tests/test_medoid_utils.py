# -*- coding: utf-8 -*-

import unittest
from dendropy import Tree

from medoids.medoid_utils import binarize_tree


class TestMedoidUtils(unittest.TestCase):

    def test_binarize_tree(self):
        tree = Tree.get(data='((a:3,b:3,c:5,d:5):2,((e:2):1):2);', schema='newick')
        expected_result = '((a:3.0,(b:3.0,(c:5.0,d:5.0):0.0):0.0):2.0,e:5.0)'
        binarize_tree(tree)
        self.assertEqual(expected_result, str(tree))


if __name__ == '__main__':
    unittest.main()