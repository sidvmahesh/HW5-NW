# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    nw = NeedlemanWunsch("substitution_matrices/BLOSUM62.mat", -10, -1)
    score, seq1_align, seq2_align = nw.align(seq1, seq2)
    assert score == 4
    assert seq1_align == "MYQR"
    assert seq2_align == "M-QR"
    # Check the three matrices
    assert nw._gapB_matrix[2][3] == -11 # Since we know the second character of seq2 is a gap, after subsequently aligning the two M's
    for i in range(1, len(seq1)):
        for j in range(1, len(seq2)):
            gapA = -11 if nw._gapA_matrix[i-1][j] == 0 else - 1
            gapB = -11 if nw._gapB_matrix[i][j-1] == 0 else - 1
            match_mismatch = nw.sub_dict[(seq1[i-1], seq2[j-1])]
            assert nw._align_matrix[i][j] == max(nw._align_matrix[i-1][j] + gapA, nw._align_matrix[i][j-1] + gapB, nw._align_matrix[i-1][j-1] + match_mismatch)

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    nw = NeedlemanWunsch("substitution_matrices/BLOSUM62.mat", -10, -1)
    score, seq3_align, seq4_align = nw.align(seq3, seq4)
    assert score == 17
    assert seq3_align == "MAVHQLIRRP"
    assert seq4_align == "M---QLIRHP"
    seqA_list = [] # Each char will be | if there is no gap, and - for a gap
    seqB_list = [] # Each char will be | if there is no gap, and - for a gap
    i, j = (len(seq3), len(seq4))
    while (i != 0) and (j != 0):
        if nw._gapA_matrix[i][j] != 0:
            seqA_list.append("|")
            seqB_list.append("-")
            i = i - 1
        elif nw._gapB_matrix[i][j] != 0:
            seqA_list.append("-")
            seqB_list.append("|")
            j = j - 1
        else:
            seqA_list.append("|")
            seqB_list.append("|")
            i = i - 1
            j = j - 1
    seqA_list.reverse()
    seqB_list.reverse()
    assert len(seqA_list) == len(seqB_list) # Assert the length of the alignment strings are equal for both sequences
    for i in range(len(seqA_list)):
        assert seqA_list[i] == "-" if seq3_align[i] == "-" else "|"
        assert seqB_list[i] == "-" if seq4_align[i] == "-" else "|"

