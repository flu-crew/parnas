# -*- coding: utf-8 -*-

import unittest
from dendropy import Tree

from parnas.medoids.tree_medoids import find_coverage
from parnas.medoids import get_costs, build_distance_functions


class TestTreeCoverage(unittest.TestCase):

    def test_coverage1(self):
        tree = Tree.get(data="(((a:2,b:1):2,e:1):1,(c:1,d:2):1);", schema='newick')
        cost_map = get_costs(tree)
        coverage5 = find_coverage(tree, 5, cost_map)
        coverage4 = find_coverage(tree, 4, cost_map)
        coverage3 = find_coverage(tree, 3, cost_map)
        self.assertEqual(len(coverage5), 1)
        self.assertEqual(len(coverage4), 2)
        self.assertEqual(len(coverage3), 3)

    def test_impossible_cover1(self):
        tree = Tree.get(data="(((a:2,b:1):2,e:1):1,(c:1,d:2):1);", schema='newick')
        cost_map = get_costs(tree, excluded=['a', 'b'])
        coverage4 = find_coverage(tree, 4, cost_map)
        print(coverage4)
        self.assertTrue(coverage4 is None)

    def test_impossible_cover2(self):
        tree = Tree.get(data="(((a:2,b:1):2,e:1):1,(c:1,d:2):1);", schema='newick')
        cost_map = get_costs(tree, excluded=['a'])
        coverage = find_coverage(tree, 2, cost_map)
        print(coverage)
        self.assertTrue(coverage is None)

    def test_prior_full(self):
        tree = Tree.get(data="((a:1,(b:2.5,c:2.5):1):3,(d:0.5,(e:2.5,f:3.5):1):2);", schema='newick')
        radius = 6  # everything is covered by a and f
        prior = ['a', 'f']
        cost_map = get_costs(tree)
        coverage = find_coverage(tree, radius, cost_map=cost_map, prior_centers=prior)
        print(coverage)
        self.assertEqual(len(coverage), 0)

    def test_prior_plus_one(self):
        tree = Tree.get(data="((a:1,(b:2.5,c:2.5):1):3,(d:0.5,(e:2.5,f:3.5):1):2);", schema='newick')
        radius = 6  # everything can covered by a and f
        prior = ['a']
        cost_map = get_costs(tree)
        coverage = find_coverage(tree, radius, cost_map=cost_map, prior_centers=prior)
        print(coverage)
        self.assertEqual(len(coverage), 1)

    def test_coverage_w_obj_rep_exclusion(self):
        tree = Tree.get(data="(((a:1,(R1:1,(R2:1,R3:1):1):0.5):0.5,b:1):2,((c:0.5,R4:2):0.5,d:1):1);", schema='newick')
        radius = 3.7
        cost_map = get_costs(tree, excluded=['a', 'b', 'c', 'd'])
        coverage = find_coverage(tree, radius, cost_map=cost_map, obj_excluded=['R1', 'R2', 'R3', 'R4'])
        self.assertEqual(len(coverage), 2)
        self.assertEqual(set(coverage), {'R1', 'R4'})