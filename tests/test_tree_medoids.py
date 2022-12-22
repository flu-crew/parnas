# -*- coding: utf-8 -*-

import unittest
from dendropy import Tree

from parnas.medoids import find_n_medoids, build_distance_functions, get_costs, find_n_medoids_with_diversity


class TestTreeMedoids(unittest.TestCase):

    def test_one_medoid(self):
        tree = Tree.get(data="(((a:2,b:1):2,e:1):1,(c:1,d:2):1);", schema='newick')
        distance_funcs = build_distance_functions(tree)
        cost_map = get_costs(tree)
        medoids, obj = find_n_medoids(tree, 1, distance_funcs, cost_map)
        print(medoids, obj)
        self.assertEqual(len(medoids), 1)
        self.assertEqual(medoids[0], 'e')
        self.assertEqual(obj, 18)

    def test_two_medoids(self):
        tree = Tree.get(data="((a:1,(b:2.5,c:2.5):1):3,(d:0.5,(e:2.5,f:3):1):2);", schema='newick')
        distance_funcs = build_distance_functions(tree)
        cost_map = get_costs(tree)
        medoids, obj, diversity = find_n_medoids_with_diversity(tree, 2, distance_functions=distance_funcs, cost_map=cost_map)
        self.assertEqual(len(medoids), 2)
        self.assertEqual(obj, 17.5)
        self.assertEqual(set(medoids), {'a', 'd'})
        self.assertAlmostEqual(diversity[-1], 46.96969, places=4)

    def test_two_medoids_w_radius1(self):
        tree = Tree.get(data="((a:1,(b:2.5,c:2.5):1):3,(d:0.5,(e:2.5,f:3):1):2);", schema='newick')
        radius = 4  # d should cover e.
        distance_funcs = build_distance_functions(tree, radius=radius)
        cost_map = get_costs(tree)
        medoids, obj = find_n_medoids(tree, 2, distance_functions=distance_funcs, cost_map=cost_map)
        self.assertEqual(len(medoids), 2)
        self.assertEqual(obj, 1.5)
        self.assertEqual(set(medoids), {'a', 'd'})

    def test_two_medoids_w_radius2(self):
        tree = Tree.get(data="((a:1,(b:2.5,c:2.5):1):3,(d:0.5,(e:2.5,f:3.5):1):2);", schema='newick')
        radius = 4.5  # d should cover e; a should cover b and c.
        distance_funcs = build_distance_functions(tree, radius=radius)
        cost_map = get_costs(tree)
        medoids, obj = find_n_medoids(tree, 2, distance_functions=distance_funcs, cost_map=cost_map)
        self.assertEqual(len(medoids), 2)
        self.assertEqual(obj, 0.5)
        self.assertEqual(set(medoids), {'a', 'd'})

    def test_two_medoids_w_radius2_binary(self):
        tree = Tree.get(data="((a:1,(b:2.5,c:2.5):1):3,(d:0.5,(e:2.5,f:3.5):1):2);", schema='newick')
        radius = 4.5  # d should cover e; a should cover b and c.
        distance_funcs = build_distance_functions(tree, radius=radius, is_binary=True)
        cost_map = get_costs(tree)
        medoids, obj = find_n_medoids(tree, 2, distance_functions=distance_funcs, cost_map=cost_map)
        self.assertEqual(len(medoids), 2)
        self.assertEqual(obj, 1)
        self.assertEqual(set(medoids), {'a', 'd'})

    def test_medoids_w_obj_rep_exclusion(self):
        tree = Tree.get(data="(((a:1,(R1:1,(R2:1,R3:1):1):0.5):0.5,b:1):2,((c:0.5,R4:2):0.5,d:1):1);", schema='newick')
        distance_funcs = build_distance_functions(tree, fully_excluded=['R1', 'R2', 'R3', 'R4'])
        cost_map = get_costs(tree, excluded=['a', 'b', 'c', 'd'])
        medoids, obj = find_n_medoids(tree, 1, distance_functions=distance_funcs, cost_map=cost_map)
        self.assertEqual(len(medoids), 1)
        self.assertEqual(obj, 17.5)
        self.assertEqual(set(medoids), {'R1'})
        medoids, obj = find_n_medoids(tree, 2, distance_functions=distance_funcs, cost_map=cost_map)
        self.assertEqual(len(medoids), 2)
        self.assertEqual(obj, 11.5)
        self.assertEqual(set(medoids), {'R1', 'R4'})

    def test_zero_medoids_w_radius_and_prior(self):
        tree = Tree.get(data="((a:1,(b:2.5,c:2.5):1):3,(d:0.5,(e:2.5,f:3.5):1):2);", schema='newick')
        radius = 6  # everything is covered by a and f
        prior = ['a', 'f']
        distance_funcs = build_distance_functions(tree, radius=radius, prior_centers=prior)
        cost_map = get_costs(tree)
        medoids, obj = find_n_medoids(tree, 1, distance_functions=distance_funcs, cost_map=cost_map)
        print(medoids)
        self.assertEqual(obj, 0)

    def test_one_medoid_w_prior(self):
        tree = Tree.get(data="(((a:2,b:1):2,e:1):1,(c:1,d:2):1);", schema='newick')
        prior_centers = ['b']
        distance_funcs = build_distance_functions(tree, prior_centers=prior_centers)
        cost_map = get_costs(tree)
        medoids, obj = find_n_medoids(tree, 1, distance_funcs, cost_map)
        print(medoids, obj)
        self.assertEqual(len(medoids), 1)
        self.assertTrue(medoids[0] in {'c', 'd'})
        # self.assertEqual(obj, 18)

    def test_rep_exclusion(self):
        # Excluding 'e' from possible representatives.
        tree = Tree.get(data="(((a:2,b:1):2,e:1):1,(c:1,d:2):1);", schema='newick')
        distance_funcs = build_distance_functions(tree)
        cost_map = get_costs(tree, excluded=['e'])  # Exclude 'e' from reps.
        medoids, obj = find_n_medoids(tree, 1, distance_funcs, cost_map)
        print(medoids, obj)
        self.assertEqual(len(medoids), 1)
        self.assertTrue(medoids[0] in {'b', 'c'})
        self.assertEqual(obj, 20)

    def test_full_exclusion(self):
        # Fully excluding 'e' from the tree.
        tree = Tree.get(data="(((a:2,b:1):2,e:1):1,(c:1,d:2):1);", schema='newick')
        distance_funcs = build_distance_functions(tree, fully_excluded=['e'])
        cost_map = get_costs(tree, excluded=['e'])
        medoids, obj = find_n_medoids(tree, 1, distance_funcs, cost_map)
        print(medoids, obj)
        self.assertEqual(len(medoids), 1)
        self.assertTrue(medoids[0] in {'b', 'c'})
        self.assertEqual(obj, 16)

    def test_zero_weight(self):
        tree = Tree.get(data="(((a:2,b:1):2,e:1):1,(c:0.5,d:2):1);", schema='newick')
        weights = {'a': 0}
        distance_funcs = build_distance_functions(tree, taxa_weights=weights)
        cost_map = get_costs(tree)
        medoids, obj = find_n_medoids(tree, 1, distance_funcs, cost_map)
        self.assertAlmostEqual(obj, 11.5, 10)
        self.assertEqual(medoids[0], 'c')


if __name__ == '__main__':
    unittest.main()
