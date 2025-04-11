# -*- coding: utf-8 -*-

import colorsys
from datetime import datetime
from typing import List
import sys
import random as rnd
# import math

# from Bio import SeqRecord
from dendropy import Tree
from scipy.stats import percentileofscore

from parnas import parnas_logger
from parnas.options import parser, parse_and_validate
# from parnas.sequences import SequenceSimilarityMatrix
from parnas.medoids import find_n_medoids, annotate_with_closest_centers, build_distance_functions, binarize_tree,\
    get_costs, find_n_medoids_with_diversity, find_coverage
from parnas.medoids.medoid_utils import get_centers_score, compute_percent_coverage


# os.environ["NUMBA_DUMP_ANNOTATION"] = "1"
sys.setrecursionlimit(100000)


def color_by_clusters(tree: Tree, centers: List[str], prior_centers=None, fully_excluded=None, radius=None):
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
        if edge.head_node:
            edge.head_node.annotations.drop(name='!color')  # remove any previous color scheme.
        # If the two ends of the edge are covered by the same center -- color it.
        if edge.head_node and edge.tail_node and edge.head_node.annotations.get_value('center') is not None and\
                edge.head_node.annotations.get_value('center') == edge.tail_node.annotations.get_value('center'):
            if edge.head_node.annotations.get_value('center') < 0:  # negative if it is covered by a prior center.
                edge.head_node.annotations.add_new('!color', prior_color)  # use the reserved color.
            else:
                center_ind = edge.head_node.annotations.get_value('center')
                color = hex_colors[center_ind]
                edge.head_node.annotations.add_new('!color', color)

    # Drop previous color annotations and color 'regular' taxa black (to avoid color mixing in FigTree).
    black = '#000000'
    grey = '#a9a9a9'
    special_taxa = centers
    if prior_centers:
        special_taxa = special_taxa + prior_centers
    special_taxa = set(special_taxa)
    for taxon in tree.taxon_namespace:
        taxon.annotations.drop(name='!color')  # Remove any previous colors if any.
        if not taxon.label in special_taxa:
            if fully_excluded and taxon.label in fully_excluded:
                taxon.annotations.add_new('!color', grey)
            else:
                taxon.annotations.add_new('!color', black)

    # Color the centers.
    for center_ind, center in enumerate(centers):
        taxon = tree.taxon_namespace.get_taxon(center)
        color = hex_colors[center_ind]
        taxon.annotations.add_new('!color', color)

    # Color prior centers if applicable.
    if prior_centers:
        for prior_center in prior_centers:
            taxon = tree.taxon_namespace.get_taxon(prior_center)
            taxon.annotations.add_new('!color', prior_color)


def save_clusters(clusters_path: str, tree: Tree, centers: List[str], prior_centers=None, fully_excluded=None,
                  radius=None, annotated=False):
    if not annotated:
        annotate_with_closest_centers(tree, centers, prior_centers=prior_centers, radius=radius)

    clusters = {}
    for leaf in tree.leaf_nodes():  # Go through the annotations and collect clusters together.
        if fully_excluded and leaf.taxon.label in fully_excluded:
            continue
        if leaf.annotations.get_value('center') is not None and leaf.annotations.get_value('center') >= 0:
            center: int = leaf.annotations.get_value('center')
            if center in clusters.keys():
                cluster: List[str] = clusters.get(center)
                cluster.append(leaf.taxon.label)
            else:
                clusters[center] = [leaf.taxon.label]

    with open(clusters_path, 'w') as clusters_file:  # Save clusters (partition) to a user-specified path.
        for center in sorted(clusters.keys()):
            for taxon in clusters[center]:
                clusters_file.write(f'{taxon}\t{center}\n')  # <taxon-name><tab><partition-number>
    parnas_logger.info(f'Saved clusters to {clusters_path}')


# def assess_clade_similarity(annotated_tree: Tree, centers: List[str], records: List[SeqRecord.SeqRecord],
#                             sim_matrix: SequenceSimilarityMatrix):
#     for i, center in enumerate(centers):
#         print('Clade of %s:' % center)
#         max_within_dist = 0
#         min_outside_dist = math.inf
#         for leaf1 in annotated_tree.leaf_nodes():
#             if leaf1.taxon.label != center:
#                 continue
#             for leaf2 in annotated_tree.leaf_nodes():
#                 leaf1_ind = [i for i, rec in enumerate(records) if rec.id == leaf1.taxon.label][0]
#                 leaf2_ind = [i for i, rec in enumerate(records) if rec.id == leaf2.taxon.label][0]
#                 dist = 1 - sim_matrix.matrix[leaf1_ind][leaf2_ind]
#                 if leaf1.annotations.get_value('center') == i and leaf2.annotations.get_value('center') == i:
#                     if dist > max_within_dist:
#                         max_within_dist = dist
#                 if leaf1.annotations.get_value('center') == i and leaf2.taxon.label in centers and leaf1 != leaf2:
#                     if dist < min_outside_dist:
#                         min_outside_dist = dist
#         print('\tmax within clade: %.3f%%' % (max_within_dist * 100))
#         print('\tclosest other representative: %.3f%%' % (min_outside_dist * 100))


def run_parnas_cli():
    args, query_tree, n, radius, is_binary, prior_centers, excluded_taxa, obj_excluded, fully_excluded, taxa_weights = parse_and_validate()

    # Binarize the query tree:
    binarize_tree(query_tree, edge_length=0)

    cost_map = get_costs(query_tree, excluded_taxa, fully_excluded)
    if not args.evaluate:
        # The main procedure for finding best representatives.
        dist_functions = build_distance_functions(query_tree, prior_centers=prior_centers, is_binary=is_binary,
                                                  fully_excluded=fully_excluded + obj_excluded, radius=radius,
                                                  taxa_weights=taxa_weights)
        parnas_logger.info("Inferring best representatives...")
        if not args.cover:
            representatives, value, diversity_scores, _ = find_n_medoids_with_diversity(query_tree, n, dist_functions,
                                                                                        cost_map, max_dist=radius)
        else:
            # coverage = None
            coverage = find_coverage(query_tree, radius, cost_map, prior_centers, fully_excluded, obj_excluded)
            if coverage is None:
                parnas_logger.warning("The tree cannot be fully covered given the exclusion constraints.")
                parnas_logger.warning("Falling back onto a slower method that would cover the tree as much as possible.")
                opt_n = -1
                prev_value = -1
                for n in range(1, len(query_tree.leaf_nodes()) + 1):
                    representatives, value = find_n_medoids(query_tree, n, dist_functions, cost_map, max_dist=radius)
                    if value == prev_value:
                        # The best possible value was achieved on the previous n (full coverage is impossible).
                        opt_n = n - 1
                        break
                    elif value == 0:
                        # Achieved full coverage.
                        opt_n = n
                        break
                    prev_value = value
            else:
                representatives = coverage

        if len(representatives) == 0:
            parnas_logger.info('The diversity on the tree is already fully covered by the prior centers - no new '
                               'representatives needed.')
        else:
            parnas_logger.info('Chosen representatives:')
            for rep in representatives:
                print('%s' % rep)
        if not args.cover and len(representatives) > 1 and diversity_scores:
            parnas_logger.info('Chosen representatives account for %.2f%% of ' % diversity_scores[-1] +
                               f'{"(new) " if prior_centers else "overall "}diversity.')

        annotated = False
        if args.out_path:
            color_by_clusters(query_tree, representatives, prior_centers=prior_centers, fully_excluded=fully_excluded,
                              radius=radius)
            annotated = True
            parnas_logger.debug(f'Finished coloring {datetime.now().strftime("%H:%M:%S")}')
            try:
                query_tree.write(path=args.out_path, schema='nexus')
                parnas_logger.info('Colored tree was saved to "%s".' % args.out_path)
            except Exception:
                parser.error('Cant write to the specified path "%s".' % args.out_path)

        if args.clusters_path:
            save_clusters(args.clusters_path, query_tree, representatives, prior_centers=prior_centers,
                          fully_excluded=fully_excluded, radius=radius, annotated=annotated)

        if args.csv_path:
            if args.cover:
                parnas_logger.warning('"--diversity" cannot be used in combination with "--cover".')
            elif n <= 1:
                parnas_logger.warning('No diversity scores to report as n < 2.')
            elif diversity_scores:
                try:
                    with open(args.csv_path, 'w') as diversity_log:
                        diversity_log.write('Representatives, Diversity_covered\n')
                        for i, score in enumerate(diversity_scores):
                            diversity_log.write('%d, %.2f\n' % (i + 2, score))
                    parnas_logger.info('Diversity scores were saved to "%s".' % args.csv_path)
                except IOError:
                    parser.error('Cant write to the specified path "%s".' % args.csv_path)

        if args.sample_tree_path:
            retain_taxa = representatives.copy()
            if args.include_prior:
                retain_taxa.extend(prior_centers)
            sample_tree: Tree = query_tree.extract_tree_with_taxa_labels(retain_taxa)
            sample_tree.purge_taxon_namespace()
            for taxon in sample_tree.taxon_namespace:
                taxon.annotations.drop(name='!color')
            sample_tree.write(path=args.sample_tree_path, schema='nexus')
            parnas_logger.info('The subtree was saved to "%s".' % args.sample_tree_path)

            # if args.nt_alignment:
            #     # Parsing assuming its swIAV HA sequences.
            #     alignment_path = args.nt_alignment
            #     records = list(SeqIO.parse(alignment_path, 'fasta'))
            #     sim_matrix = SequenceSimilarityMatrix(None, None, records, records, aligned=True, ignore_tails=True, nt_alphabet=True)
            #     assess_clade_similarity(query_tree, representatives, records, sim_matrix)
    else:
        if args.cover:
            # Evaluate prior reps based on the % of taxa they cover
            coverage = compute_percent_coverage(query_tree, prior_centers, radius, fully_excluded + obj_excluded)
            parnas_logger.info(
                f'Prior representative(s) cover {round(coverage * 100, 2)}% of taxa within the specified radius.')

            if args.out_path:
                color_by_clusters(query_tree, prior_centers, prior_centers=[],
                                  fully_excluded=fully_excluded,
                                  radius=radius)
                try:
                    query_tree.write(path=args.out_path, schema='nexus')
                    parnas_logger.info('Colored tree was saved to "%s".' % args.out_path)
                except Exception:
                    parser.error('Cant write to the specified path "%s".' % args.out_path)
        else:
            # The procedure for evaluating the prior centers for the non-coverage case.

            dist_functions = build_distance_functions(query_tree, is_binary=is_binary,
                                                      fully_excluded=fully_excluded + obj_excluded, radius=radius,
                                                      taxa_weights=taxa_weights)
            prior_score = get_centers_score(query_tree, prior_centers, dist_functions)

            # Find the best n representatives
            n = len(prior_centers)
            parnas_logger.info(f'Finding best {n} representative(s) to compare with prior...')
            representatives, value, diversity_scores, obj1 = find_n_medoids_with_diversity(query_tree, n, dist_functions,
                                                                                           cost_map, max_dist=radius)

            if obj1 <= 0:
                parser.error('It seems there is no diversity to be covered. '
                             'Please check your tree and the --radius/--threshold setting (if any).')

            if n > 1:
                best_diversity = diversity_scores[-1]
                # Compare the prior reps with the best n reps.
                prior_diversity = (obj1 - prior_score) / obj1 * 100

                # decrease_from_optimal = (1 - value / prior_score) * 100
                # parnas_logger.info('Prior representatives are %.2f%% less representative than the optimal set.'
                #                    % decrease_from_optimal)

                # rnd_diversity = [(obj1 - rnd_score) / obj1 * 100 for rnd_score in rnd_scores]

                parnas_logger.plain('')
                parnas_logger.info('The amount of diversity covered by the prior and best sets:')
                if prior_diversity > 0:
                    parnas_logger.info('PRIOR:\t%.2f%%' % prior_diversity)
                else:
                    parnas_logger.info('PRIOR:\t%.2f%%   (Note that this value is negative because the prior set is less representative than a single best taxon)' % prior_diversity)
                parnas_logger.info('BEST:\t\t%.2f%%' % best_diversity)
                # parnas_logger.info('RANDOM (mean):\t%.2f%%', sum(rnd_diversity) / len(rnd_diversity))
                # parnas_logger.plain('')
            else:
                # Evaluate the single representative.
                print(f'Best score: {obj1}. Prior score: {prior_score}')
                parnas_logger.info(f'The prior representative is {round(100 - (prior_score - obj1) / prior_score * 100, 3)}% '
                                   f'as representative as the best one.')

            # Compare prior to random sets.
            replicates = 1000
            parnas_logger.info(f'Comparing the prior representatives with random reps ({replicates} replicates)...')
            taxa_labels = [leaf.taxon.label for leaf in query_tree.leaf_nodes() if
                           leaf.taxon.label not in fully_excluded and leaf.taxon.label not in excluded_taxa]
            rnd_scores = []
            for i in range(replicates):
                rnd.shuffle(taxa_labels)
                rnd_reps = taxa_labels[:n]
                rnd_score = get_centers_score(query_tree, rnd_reps, dist_functions)
                rnd_scores.append(rnd_score)
            rnd_scores = sorted(rnd_scores)
            prior_percentile = percentileofscore(rnd_scores, prior_score, kind='strict')
            parnas_logger.info(f'Prior strains are more representative than {100 - prior_percentile}% of random sets.')
