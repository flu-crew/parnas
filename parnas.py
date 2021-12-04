#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import subprocess
import sys
from argparse import RawTextHelpFormatter
from io import StringIO
from typing import List
from dendropy import Tree, Node, Edge
from Bio import Phylo, AlignIO
from Bio.Align import MultipleSeqAlignment
from math import floor
import re
import colorsys
from random import shuffle

from medoids import find_n_medoids, annotate_with_closest_centers, build_distance_functions


# Computes the coverage radius (# of substitutions) that satisfies the similarity threshold.
def threshold_to_substitutions(sim_threshold: float, alignment: MultipleSeqAlignment) -> int:
    subs = floor((1 - sim_threshold / 100) * len(alignment))
    print("%3f similarity threshold implies that a single representative will cover all tips "
          "in the %d-substitution radius." % (sim_threshold, subs))
    return subs


def reweigh_tree_ancestral(tree_path: str, alignment_path: str, aa=False) -> Tree:
    """
    Re-weighs the tree edges according to the number of substitutions per edge.
    The ancestral substitutions are inferred using TreeTime (Sargulenko et al. 2018).
    :param tree_path: path to the tree to be re-weighed.
    :param alignment_path: path to MSA associated with the tree tips.
    :param aa: whether the sequences consist of amino acid residues (default: nucleotide).
    :return: The re-weighed tree
    """
    # Run ancestral inference with treetime.
    treetime_outdir = 'treetime_ancestral_%s' % tree_path
    command = ['treetime', 'ancestral', '--aln', alignment_path, '--tree', tree_path, '--outdir', treetime_outdir,
               '--gtr', 'infer']
    if aa:
        command += ['--aa']
    treetime_out = '%s/ancestral_tree.nexus' % treetime_outdir
    subprocess.call(command)

    # Update the treetime output file.
    treetime_for_dendropy = '%s/ancestral_updated.nexus' % treetime_outdir
    with open(treetime_out, 'r') as treetime_nexus:
        with open(treetime_for_dendropy, 'w') as dendropy_nexus:
            for line in treetime_nexus:
                for mutations in re.findall(r'mutations=".*?"', line):
                    upd_mutations = mutations.replace(',', '___')  # DendroPy does not like commas in the annotations.
                    line = line.replace(mutations, upd_mutations)
                dendropy_nexus.write(line)

    # Read the treetime output and weight the edges according to the number of subs.
    try:
        ancestral_tree = Tree.get(path=treetime_out, schema='nexus')
    except:
        print('Failed to infer an ancestral tree with TreeTime.'
              'Please see the TreeTime output log and consider inferring the ancestral states manually.')
        sys.exit(-1)
    reweighed_tree = ancestral_tree
    for node in reweighed_tree.node():
        edge_length = 0
        mutations_str = node.annotations.get_value('mutations')
        if mutations_str and mutations_str.strip():
            edge_length = mutations_str.count('___') + 1
        node.edge_length = edge_length
    return reweighed_tree


def color_by_clusters(tree: Tree, centers: List[str], prior_centers=None, radius=None):
    # Annotate the tree nodes with the closest centers.
    annotate_with_closest_centers(tree, centers, prior_centers=prior_centers, radius=radius)

    # Choose n different colors. If prior_centers list is not empty, add an additional reserved color.
    n = len(centers) if not prior_centers else len(centers) + 1
    hsv_colors = [(i / n, 0.7, 0.7) for i in range(n)]
    rgb_colors = [colorsys.hsv_to_rgb(*hsv) for hsv in hsv_colors]
    rgb_colors = [[round(255 * color[i]) for i in range(len(color))] for color in rgb_colors]
    hex_colors = ['#%02x%02x%02x' % (color[0], color[1], color[2]) for color in rgb_colors]
    if prior_centers:
        prior_color = hex_colors[0]  # Reserve the first color for prior centers.
        hex_colors = hex_colors[1:]

    # Color the edges.
    for edge in tree.preorder_edge_iter():
        # If the two ends of the edge are covered by the same center -- color it.
        if edge.head_node and edge.tail_node and\
                edge.head_node.annotations.get_value('center') == edge.tail_node.annotations.get_value('center'):
            if edge.head_node.annotations.get_value('center') < 0:  # negative if it is covered by a prior center.
                edge.head_node.annotations.add_new('!color', prior_color)  # use the reserved color.
            else:
                center_ind = edge.head_node.annotations.get_value('center')
                color = hex_colors[center_ind]
                edge.head_node.annotations.add_new('!color', color)

    # Color the centers.
    for center_ind, center in enumerate(centers):
        taxon = tree.taxon_namespace.get_taxon(center)
        color = hex_colors[center_ind]
        taxon.annotations.add_new('!color', color)

    # Color prior centers if applicable.
    for prior_center in prior_centers:
        taxon = tree.taxon_namespace.get_taxon(prior_center)
        taxon.annotations.add_new('!color', prior_color)


if __name__ == '__main__':
    # Program interface:
    parser = argparse.ArgumentParser(description='Phylogenetic mAximum RepreseNtAtion Sampling (PARNAS)',
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument('-t', '--tree', type=str, action='store', dest='tree',
                        help='path to the input tree in newick or nexus format', required=True)
    parser.add_argument('-n', type=int, action='store', dest='samples',
                        help='number of samples (representatives) to be chosen.\n' +
                             'This argument is required unless the --cover option is specified', required=False)
    parser.add_argument('--color', type=str, action='store', dest='out_path',
                        help='PARNAS will save a colored tree, where the chosen representatives are highlighted '
                        'and the tree is color-partitioned respective to the representatives.\n'
                        'If prior centers are specified, they (and the subtrees they represent) will be colored red.')
    parser.add_argument('--prior-regex', type=str, action='store', dest='prior_regex',
                        help='indicate the previous centers (if any) with a regex. '
                             'The regex should match the full taxon name.\n'
                             'PARNAS will then select centers that represent diversity'
                             'not covered by the previous centers.', required=False)
    # parser.add_argument('--threshold', type=float, action='store', dest='percent',
    #                     help='sequences similarity threshold: the algorithm will choose best representatives that cover as much\n' +
    #                          'diversity as possible within the given similarity threshold.' +
    #                          '--nt or --aa must be specified with this option', required=False)
    # parser.add_argument('--nt', type=str, action='store', dest='nt_alignment',
    #                     help='path to nucleotide sequences associated with the tree tips', required=False)
    # parser.add_argument('--aa', type=str, action='store', dest='aa_alignment',
    #                     help='path to amino acid sequences associated with the tree tips', required=False)
    # parser.add_argument('--cover', action='store_true',
    #                     help="choose the best representatives (smallest number) that cover all the tips within the specified threshold.\n" +
    #                     "If specified, the --threshold argument must be specified as well",
    #                     required=False)
    # parser.add_argument('--prior', metavar='TAXON', type=str, nargs='+',
    #                     help='space-separated list of taxa that have been previously chosen as centers.\n' +
    #                          'The algorithm will choose new representatives that cover the "new" diversity in the tree')
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

    # Validate n.
    n = args.samples
    if n < 1 or n >= len(tree.taxon_namespace):
        parser.error('n should be at least 1 and smaller than the number of taxa in the tree.')

    # Handle --prior-regex.
    prior_centers = None
    if args.prior_regex:
        prior_regex = args.prior_regex
        prior_centers = []
        print('Prior centers that match the regex:')
        for taxon in tree.taxon_namespace:
            if re.match(prior_regex, taxon.label):
                prior_centers.append(taxon.label)
                print('\t%s' % taxon.label)
        if not prior_centers:
            print('\tNone.')

    # Validate threshold-related parameters and re-weigh the tree
    if args.percent:
        if args.percent <= 0 or args.percent >= 100:
            parser.error('Invalid "--threshold %.3f" option. The threshold must be between 0 and 100 (exclusive)'
                         % args.percent)
        if args.nt_alignment and args.aa_alignment:
            parser.error('Please specify EITHER the nucleotide or amino acid alignment - not both.')
        if not (args.nt_alignment or args.aa_alignment):
            parser.error('To use the --threshold parameter please specify a nucleotide' +
                         'or amino acid alignment associated with the tree tips.')
        alignment_path = args.nt_alignment if args.nt_alignment else args.aa_alignment
        is_aa = args.aa_alignment is not None
        try:
            alignment = AlignIO.read(alignment_path, 'fasta')
        except:
            parser.error('Cannot read the specified FASTA alignment in "%s".' % alignment_path)
        radius = threshold_to_substitutions(args.percent, alignment)
        query_tree = reweigh_tree_ancestral(args.tree, alignment_path, is_aa)
    else:
        query_tree = tree
        radius = None

    bio_tree = Phylo.read(StringIO(str(query_tree) + ';'), 'newick')  # convert the denropy tree to a biopython tree.
    dist_functions = build_distance_functions(bio_tree, prior_centers=prior_centers, radius=radius)
    representatives, value = find_n_medoids(bio_tree, n, dist_functions, max_dist=radius)

    print('Chosen representatives:')
    for rep in representatives:
        print('\t%s' % rep)

    # Choose random centers for testing.
    # taxa = [taxon.label for taxon in query_tree.taxon_namespace]
    # shuffle(taxa)
    # rnd_representatives = taxa[:n]
    # print(rnd_representatives)

    if args.out_path:
        color_by_clusters(query_tree, representatives, prior_centers=prior_centers, radius=radius)
        try:
            query_tree.write(path=args.out_path, schema='nexus')
            print('Colored tree was saved to "%s".' % args.out_path)
        except:
            parser.error('Cant write to the specified path "%s".' % args.out_path)
