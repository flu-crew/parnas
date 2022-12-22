# -*- coding: utf-8 -*-
import argparse
import os
import re
import subprocess
from argparse import RawTextHelpFormatter, RawDescriptionHelpFormatter
from math import floor

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from dendropy import Tree

from parnas import parnas_logger

# Program interface:
parser = argparse.ArgumentParser(description='Phylogenetic mAximum RepreseNtAtion Sampling (PARNAS)',
                                 formatter_class=RawDescriptionHelpFormatter)
parser._optionals.title = "Arguments"
parser.add_argument('-t', '--tree', type=str, action='store', dest='tree',
                    help='Path to the input tree in newick or nexus format.', required=True)
parser.add_argument('-n', type=int, action='store', dest='samples',
                    help='Number of representatives to be chosen.\n' +
                         'This argument is required unless the --cover option is specified')
parser.add_argument('--prior', type=str, action='store', dest='prior_regex',
                    help='Indicate the previous representatives (if any) with a regex. '
                         'The regex should match a full taxon name.\n'
                         'PARNAS will then select centers that represent diversity '
                         'not covered by the previous representatives.', required=False)
# parser.add_argument('--threshold', type=float, action='store', dest='percent',
#                     help='sequences similarity threshold: the algorithm will choose best representatives that cover as much\n' +
#                          'diversity as possible within the given similarity threshold. ' +
#                          '--nt or --aa must be specified with this option', required=False)
parser.add_argument('--weights', type=str, action='store', dest='weights_csv',
                    help='A CSV file specifying a weight for some or all taxa. '
                         'The column names must be "taxon" and "weight".\n'
                         'If a taxon is not listed in the file, its weight is assumed to be 1. '
                         'Maximum allowed weight is 1000 and weights below 1e-8 are considered 0.')
parser.add_argument('--radius', type=float, action='store', dest='radius',
                    help='Each representative will "cover" all leaves within the specified radius on the tree. '
                         'PARNAS will then choose representatives so that the amount of uncovered diversity is minimized.',
                    required=False)
parser.add_argument('--cover', action='store_true',
                    help="Choose the best representatives (smallest number) that cover all the tips within the specified radius/threshold.\n" +
                    "If specified, a --radius or --threshold argument must be specified as well.",
                    required=False)
parser.add_argument('--binary', action='store_true',
                    help="To be used with --radius. Instead of covering as much diversity as possible, "
                         "PARNAS will cover as many tips as possible within the radius. "
                         "Each leaf will have a binary contribution to the objective: 0 if covered, else its weight.",
                    required=False)

output_options = parser.add_argument_group('Output options')
output_options.add_argument('--color', type=str, action='store', dest='out_path',
                    help='PARNAS will save a colored tree, where the chosen representatives are highlighted '
                    'and the tree is color-partitioned respective to the representatives.\n'
                    'If prior centers are specified, they (and the subtrees they represent) will be colored red.')
output_options.add_argument('--diversity', type=str, action='store', dest='csv_path',
                    help='Save diversity scores for all k (number of representatives) from 2 to n.\n'
                         'Can be used to choose the right number of representatives for a dataset.')
output_options.add_argument('--subtree', type=str, action='store', dest='sample_tree_path',
                    help='Prune the tree to the sampled taxa and save to the specified file in NEXUS format.')

taxa_handler = parser.add_argument_group('Excluding taxa')
taxa_handler.add_argument('--exclude-rep', type=str, action='store', dest='exclude_regex',
                          help='Prohibits parnas to choose representatives from the taxa matching this regex. '
                               'However, the excluded taxa will still contribute to the objective function.')
taxa_handler.add_argument('--exclude-obj', type=str, action='store', dest='eobj_regex',
                          help='Matching taxa can be selected, but will not contribute to the objective function. '
                               'Can be used if one wants to select taxa from a reference set.')
taxa_handler.add_argument('--exclude-fully', type=str, action='store', dest='full_regex',
                          help='Completely ignore the taxa matching this regex.')
taxa_handler.add_argument('--constrain-fully', type=str, action='store', dest='constrain_regex',
                          help='Completely ignore the taxa NOT matching this regex.')

alignment_parser = parser.add_argument_group('Controlling sequence divergence')
alignment_parser.add_argument('--threshold', type=float, action='store', dest='percent',
                              help='The sequence similarity threshold that works as --radius. ' +
                                   'For example, "95" will imply that each representative ' +
                                   'covers all leaves within 5%% divergence on the tree.\n'
                                   'To account for sequence divergence, parnas will run TreeTime to infer ancestral substitutions '
                                   'along the tree edges and re-weigh the edges based on the number of sustitutions on them.\n'
                                   'A sequence alignment (--nt or --aa) must be specified with this option',
                              required=False)
alignment_parser.add_argument('--nt', type=str, action='store', dest='nt_alignment',
                    help='Path to nucleotide sequence alignment associated with the tree tips.', required=False)
alignment_parser.add_argument('--aa', type=str, action='store', dest='aa_alignment',
                    help='Path to amino acid sequence alignment associated with the tree tips.', required=False)
# parser.add_argument('--prior', metavar='TAXON', type=str, nargs='+',
#                     help='space-separated list of taxa that have been previously chosen as centers.\n' +
#                          'The algorithm will choose new representatives that cover the "new" diversity in the tree')


# Computes the coverage radius (# of substitutions) that satisfies the similarity threshold.
def threshold_to_substitutions(sim_threshold: float, alignment: MultipleSeqAlignment) -> int:
    subs = floor((1 - sim_threshold / 100) * len(alignment[0]))
    parnas_logger.info("%.3f%% similarity threshold implies that a single representative will cover all tips "
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
    treetime_outdir = 'treetime_ancestral_%s' % tree_path.split(os.sep)[-1]
    if not os.path.exists(treetime_outdir):
        os.mkdir(treetime_outdir)
    treetime_log_path = '%s/treetime.log' % treetime_outdir
    parnas_logger.info('Inferring ancestral substitutions with TreeTime. The log will be written to "%s".' % treetime_log_path)
    command = ['treetime', 'ancestral', '--aln', alignment_path, '--tree', tree_path, '--outdir', treetime_outdir,
               '--gtr', 'infer']
    if aa:
        command += ['--aa']
    treetime_out = '%s/annotated_tree.nexus' % treetime_outdir
    with open(treetime_log_path, 'w') as treetime_log:
        subprocess.call(command, stdout=treetime_log, stderr=subprocess.STDOUT)

    # Read the treetime output and weight the edges according to the number of subs.
    try:
        # Update the treetime output file.
        treetime_for_dendropy = '%s/ancestral_updated.nexus' % treetime_outdir
        with open(treetime_out, 'r') as treetime_nexus:
            with open(treetime_for_dendropy, 'w') as dendropy_nexus:
                for line in treetime_nexus:
                    for mutations in re.findall(r'mutations=".*?"', line):
                        upd_mutations = mutations.replace(',',
                                                          '||')  # DendroPy does not like commas in the annotations.
                        line = line.replace(mutations, upd_mutations)
                    dendropy_nexus.write(line)
        ancestral_tree = Tree.get(path=treetime_for_dendropy, schema='nexus', preserve_underscores=True)
    except Exception:
        parser.error('Failed to infer an ancestral substitutions with TreeTime. '
                     'Please see the TreeTime output log.')

    parnas_logger.info('Re-weighing the tree based on ancestral substitutions.')
    reweighed_tree_path = '%s/ancestral_reweighed.tre' % treetime_outdir
    reweighed_tree = ancestral_tree
    for node in reweighed_tree.nodes():
        edge_length = 0
        mutations_str = node.annotations.get_value('mutations')
        if mutations_str and mutations_str.strip():
            edge_length = mutations_str.count('||') + 1
        node.edge_length = edge_length
    reweighed_tree.write(path=reweighed_tree_path, schema='newick')
    return reweighed_tree


def validate_weights(path: str):
    taxa_weights = {}
    try:
        with open(path, 'r') as weights_csv:
            headers = [header.strip() for header in weights_csv.readline().split(',')]
            if len(headers) != 2 or headers[0] != 'taxon' or headers[1] != 'weight':
                parser.error('Invalid weights file: the first line of the scv file must be "taxon, weight".')
            while True:
                l = weights_csv.readline().strip()
                if not l:
                    break
                taxon, weight_str = [value.strip() for value in l.split(',')]
                weight = float(weight_str)
                if weight < 0 or weight > 1000:
                    parser.error('Invalid weight %.8f for taxon %s. All weights must be non-negative and not exceed 1000.'
                                 % (weight, taxon))
                taxa_weights[taxon] = weight
    except Exception:
        parser.error('Cannot open/read the specified weights file "%s".' % path)
    return taxa_weights


def find_matching_taxa(tree: Tree, regex: str, title: str, none_message: str, print_taxa=True):
    matching_taxa = []
    for taxon in tree.taxon_namespace:
        if re.match(regex, taxon.label):
            matching_taxa.append(taxon.label)

    if print_taxa:
        if matching_taxa:
            parnas_logger.info(title)
            for t in matching_taxa:
                parnas_logger.plain('\t%s' % t)
        else:
            parnas_logger.info(none_message)
        parnas_logger.plain('')
    return matching_taxa


def parse_and_validate():
    args = parser.parse_args()

    # Validate the tree.
    tree = None
    try:
        tree = Tree.get(path=args.tree, schema='newick', preserve_underscores=True)
    except Exception:
        try:
            tree = Tree.get(path=args.tree, schema='nexus', preserve_underscores=True)
        except Exception:
            parser.error('Cannot read the specified tree file "%s". ' % args.tree +
                         'Make sure the tree is in the newick or nexus format.')

    # Validate n.
    n = -1
    if not args.samples and not args.cover:
        parser.error('Please either specify the number of representatives with "-n" or use the --cover option.')
    if args.samples:
        n = args.samples
        if n < 1 or n >= len(tree.taxon_namespace):
            parser.error('n should be at least 1 and smaller than the number of taxa in the tree.')

    # Handle --prior-regex.
    prior_centers = None
    if args.prior_regex:
        prior_centers = find_matching_taxa(tree, args.prior_regex, 'Prior centers that match the regex:',
                                           'No taxa matched PRIOR_REGEX', True)

    # Validate exclusions.
    excluded_taxa = []
    fully_excluded = []
    obj_excluded = []
    if args.exclude_regex:
        excluded_taxa = find_matching_taxa(tree, args.exclude_regex,
                                           'Not considering the following as representatives (matched EXCLUDE_REGEX):',
                                           'No taxa matched EXCLUDE_REGEX', True)
    if args.full_regex:
        fully_excluded = find_matching_taxa(tree, args.full_regex, 'Ignoring the following taxa (matched FULL_REGEX):',
                                            'No taxa matched FULL_REGEX', True)
    if args.eobj_regex:
        obj_excluded = find_matching_taxa(tree, args.eobj_regex,
                                          'Not contributing to the objective (matched EOBJ_REGEX):',
                                          'No taxa matched EOBJ_REGEX')
    if args.constrain_regex:
        constrained_taxa = find_matching_taxa(tree, args.constrain_regex,
                                              'Constraining to the following taxa (matched CONSTRAIN_REGEX):',
                                              'No taxa matched CONSTRAIN_REGEX', True)
        fully_excluded += list({taxon.label for taxon in tree.taxon_namespace}.symmetric_difference(constrained_taxa))

    exclude_intersection = set(excluded_taxa).intersection(set(fully_excluded))
    if exclude_intersection:
        for taxon in exclude_intersection:
            parnas_logger.warning(f'{taxon} matches both EXCLUDE_REGEX and FULL_REGEX. PARNAS will fully exclude it.')

    # Validate radius.
    radius = None
    if args.radius is not None:
        if args.radius <= 0:
            parser.error('radius should be positive.')
        else:
            radius = args.radius

    # Validate alignment.
    if args.nt_alignment or args.aa_alignment:
        if args.nt_alignment and args.aa_alignment:
            parser.error('Please specify EITHER the nucleotide or amino acid alignment - not both.')
        alignment_path = args.nt_alignment if args.nt_alignment else args.aa_alignment
        is_aa = args.aa_alignment is not None
        try:
            alignment = list(AlignIO.read(alignment_path, 'fasta'))
        except Exception:
            parser.error('Cannot read the specified FASTA alignment in "%s".' % alignment_path)
        alignment_present = True
    else:
        alignment_present = False

    # Validate threshold-related parameters and re-weigh the tree
    if args.percent:
        if args.percent <= 0 or args.percent >= 100:
            parser.error('Invalid "--threshold %.3f" option. The threshold must be between 0 and 100 (exclusive)'
                         % args.percent)
        if not alignment_present:
            parser.error('To use the --threshold parameter, please specify a nucleotide' +
                         'or amino acid alignment associated with the tree tips.')
        else:
            radius = threshold_to_substitutions(args.percent, alignment)
            query_tree = reweigh_tree_ancestral(args.tree, alignment_path, is_aa)
    else:
        query_tree = tree

    is_binary = args.binary
    if is_binary and not radius:
        is_binary = False

    # Validate cover
    if args.cover:
        if not args.percent and not args.radius:
            parser.error('To use --cover parameter, please specify --threshold or --radius option.')

    # Validate weights.
    taxa_weights = None
    if args.weights_csv:
        taxa_weights = validate_weights(args.weights_csv)

    return args, query_tree, n, radius, is_binary, prior_centers, excluded_taxa, obj_excluded, fully_excluded, taxa_weights
