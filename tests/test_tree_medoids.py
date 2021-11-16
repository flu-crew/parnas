# -*- coding: utf-8 -*-

import unittest
from Bio import Phylo
from io import StringIO

from medoids import find_n_medoids


class TestTreeMedoids(unittest.TestCase):

    def test_one_medoid(self):
        tree = Phylo.read(StringIO("((a:2,b:3):2,c:5);"), 'newick')
        medoids = find_n_medoids(tree, 1, {})
        self.assertEqual(len(medoids), 1)
        self.assertEqual(medoids[0], 'a')
        # TODO: define some actual tests...


if __name__ == '__main__':
    unittest.main()
