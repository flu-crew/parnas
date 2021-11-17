#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
from argparse import RawTextHelpFormatter
from io import StringIO

from dendropy import Tree
from Bio import Phylo
from medoids import find_n_medoids


if __name__ == '__main__':
    # Program interface:
    parser = argparse.ArgumentParser(description='Phylogenetic mAximum RepreseNtAtion Sampling (PARNAS)',
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument('-t', '--tree', type=str, action='store', dest='tree',
                        help='path to the input tree in newick or nexus format', required=True)
    parser.add_argument('-n', type=int, action='store', dest='samples',
                        help='number of samples (representatives) to be chosen.\n' +
                             'This argument is required unless the --cover option is specified', required=False)
    parser.add_argument('--threshold', type=float, action='store', dest='percent',
                        help='sequences similarity threshold: the algorithm will choose best representatives that cover as much\n' +
                             'diversity as possible within the given similarity threshold.' +
                             '--nt or --aa must be specified with this option', required=False)
    parser.add_argument('--nt', type=str, action='store', dest='nt_alignment',
                        help='path to nucleotide sequences associated with the tree tips', required=False)
    parser.add_argument('--aa', type=str, action='store', dest='aa_alignment',
                        help='path to amino acid sequences associated with the tree tips', required=False)
    parser.add_argument('--cover', action='store_true',
                        help="choose the best representatives (smallest number) that cover all the tips within the specified threshold.\n" +
                        "If specified, the --threshold argument must be specified as well",
                        required=False)
    parser.add_argument('--prior', metavar='TAXON', type=str, nargs='+',
                        help='space-separated list of taxa that have been previously chosen as centers.\n' +
                             'The algorithm will choose new representatives that cover the "new" diversity in the tree')
    args = parser.parse_args()

    # Validate the tree.
    tree = None
    try:
        tree = Tree.get(path=args.tree, schema='newick')
    except:
        try:
            tree = Tree.get(path=args.tree, schema='nexus')
        except:
            parser.error('Cannot read the specified tree file "%s". ' % args.tree +
                         'Make sure the tree is in the newick or nexus format.')
            sys.exit(-1)

    # Validate n.
    n = args.n
    if n < 1 or n >= len(tree.taxon_namespace):
        parser.error('n should be at least 1 and smaller than the number of taxa in the tree.')

    # TODO: implement handling of the other parameters.

    bio_tree = Phylo.read(StringIO(str(tree) + ';'), 'newick')  # convert the denropy tree to biopython tree.
    representatives, tip_clusters = find_n_medoids(tree, n)
