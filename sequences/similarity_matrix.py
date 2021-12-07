# -*- coding: utf-8 -*-
from typing import Tuple, List
from Levenshtein import ratio
from Bio.SeqRecord import SeqRecord

NT_ALPHABET = 'ACGT-_'
AA_ALPHABET = 'ACDEFGHIKLMNPQRSTVWY-_'


class SequenceSimilarityMatrix:

    def __init__(self, row_sequences: List[SeqRecord], col_sequences: List[SeqRecord], aligned=False, ignore_tails=True,
                 dna_alphabet=False, aa_alphabet=False):
        """
        :param protein: Name of the protein for which the similarity matrix is constructed.
        :param path: Path to the matrix (to load from or to save to)
        :param force_compute: Recompute the matrix even if a matrix was already saved
        :param aligned: The given sequences are aligned
        :param ignore_tails: Only applicable if aligned is True. Ignores missing nt/aa at the tails.
        """
        assert len(row_sequences) > 0 and len(col_sequences) > 0
        self.row_sequences = list(row_sequences)
        self.col_sequences = list(col_sequences)

        self.matrix = []
        for swine_seq in row_sequences:
            row = []
            for human_seq in col_sequences:
                # Compute the similarity between sequences and append to row.
                if aligned:
                    # Aligned distance
                    sim = 1 - self.aligned_dist(str(human_seq.seq), str(swine_seq.seq), normalized=True,
                                                ignore_tails=ignore_tails, dna_alphabet=dna_alphabet,
                                                aa_alphabet=aa_alphabet)
                else:
                    # Levenshtein ratio
                    sim = ratio(str(human_seq.seq), str(swine_seq.seq))
                row.append(sim)
            self.matrix.append(row)

    def find_n_closest_cols(self, row_ind: int, n=1, col_filter=None, latest_first=True)\
            -> [Tuple[float, SeqRecord]]:
        """
        Returns a list on n closest column sequences in the [(similarity, sequence), ...] format.
        :type col_filter: Callable[[int, IAVSeqence], bool]
        :param col_filter: Takes column index and the respective sequence and returns True
                           if the column needs to be considered
        :param latest_first: If multiple sequences have the same similarity to 'row_ind',
                             the method will sort them from earliest to latest by default.
                             Specify True to reverse that order.
        """
        # For each column that passes col_filter get (similarity, sequence) tuple:
        col_sims = [(sim, self.col_sequences[col_ind]) for col_ind, sim in enumerate(self.matrix[row_ind])
                    if (not col_filter or col_filter(col_ind, self.col_sequences[col_ind]))]
        # Sort tuples by similarity in the descending order:
        if latest_first:
            col_sims = sorted(col_sims, key=lambda x: (x[0], x[1].date), reverse=True)
        else:
            col_sims = sorted(col_sims, key=lambda x: (-x[0], x[1].date), reverse=False)
        return col_sims[:n]

    def get_row_id(self, row_seq: SeqRecord) -> int:
        return self.row_sequences.index(row_seq)

    def get_col_id(self, col_seq: SeqRecord) -> int:
        return self.col_sequences.index(col_seq)

    def get_row_id_by_name(self, name: str) -> int:
        for i, seq in enumerate(self.row_sequences):
            if seq.id == name:
                return i
        return -1

    def get_col_id_by_name(self, name: str) -> int:
        for i, seq in enumerate(self.col_sequences):
            if seq.id == name:
                return i
        return -1

    @staticmethod
    def aligned_dist(seq1: str, seq2: str, normalized=True, ignore_tails=True, dna_alphabet=False, aa_alphabet=False):
        """
        Computes the hamming distance between two aligned sequences
        :param normalized: if True, the distance will be in the [0,1] range.
                           Otherwise, the number of substitutions will be returned.
        :param ignore_tails: if True, the flanking '---' parts will be ignored on both sequences.
        :param dna_alphabet: only consider two bases different if they both are in the standard DNA alphabet.
        :param aa_alphabet: only consider two bases different if they both are in the standard AA alphabet.
        """
        # Assert that 1. sequences have the same lengths and
        # 2. both sequences contain at least one informative (not '-') symbol.
        assert len(seq1) == len(seq2)
        assert seq1.count('-') != len(seq1)
        assert seq2.count('-') != len(seq2)
        length = len(seq1)
        seq1 = seq1.upper()
        seq2 = seq2.upper()

        # Compute the number of differing sites.
        difference = 0
        for site, symbol1 in enumerate(seq1):
            symbol2 = seq2[site]
            if dna_alphabet and (symbol1 not in NT_ALPHABET or symbol2 not in NT_ALPHABET):
                continue
            if aa_alphabet and (symbol1 not in AA_ALPHABET or symbol2 not in AA_ALPHABET):
                continue
            if symbol1 != symbol2:
                difference += 1

        # Correct the count based on missing tails.
        if ignore_tails:
            seq1_start, seq2_start = 0, 0
            seq1_end, seq2_end = length - 1, length - 1
            while seq1[seq1_start] == '-':
                seq1_start += 1
            while seq2[seq2_start] == '-':
                seq2_start += 1
            while seq1[seq1_end] == '-':
                seq1_end -= 1
            while seq2[seq2_end] == '-':
                seq2_end -= 1
            for i in range(max(seq1_start, seq2_start)):
                if seq1[i] != seq2[i]:
                    difference -= 1
            for i in range(length - 1, min(seq1_end, seq2_end), -1):
                if seq1[i] != seq2[i]:
                    difference -= 1
            length = min(seq1_end, seq2_end) - max(seq1_start, seq2_start) + 1
        if normalized:
            return difference / length
        else:
            return difference
