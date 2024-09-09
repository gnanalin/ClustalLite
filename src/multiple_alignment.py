"""Program which aims at computing the multiple alignments

Usage:
------
------

"""

__authors__ = ("Stéphanie Gnanalingam")
__contact__ = ("stephanie.gnanalingam@etu.u-paris.fr")
__date__ = "2024-09-09"

import pandas as pd
import numpy as np 
import copy 
import needleman_wunsch_method

BLOSUM_MATRIX = pd.read_csv("../data/BLOSUM_MATRIX.csv", index_col=[0])
GAP = 8

def calculate_score_multiple_alignment(matrix, current_line, current_col, cluster1_seq_list, cluster2_seq_list):
    """Calculate score for the current box

    Parameters
    ----------
    matrix : a score matrix 
    current_line, current_col : the coordinates of the current box
    cluster1_seq_list, cluster2_seq_list : the list of sequences from cluster 1 and 2

    Returns
    ------
    tuple : (score, index) where index will be the box which gave 
    the maximal score (0: diagonal, 1: left, 2: up). If there are equals with
    the diagonal, the diagonal box will be selected and otherwise, it left box is selected
    """
    mean_diag_score = 0
    counter = 0
    for sequence1 in cluster1_seq_list:
        for sequence2 in cluster2_seq_list:
            counter += 1
            aa1, aa2 = sequence1[current_line-1], sequence2[current_col-1]
            if aa1 != '-' and aa2 != '-':
                mean_diag_score += BLOSUM_MATRIX.loc[aa1, aa2]
            else:
                mean_diag_score -= GAP
    mean_diag_score /= counter
    score_diag = matrix[current_line-1, current_col-1]+mean_diag_score
    score_up = matrix[current_line-1, current_col]-GAP
    score_left = matrix[current_line, current_col-1]-GAP
    scores_array = np.array([score_diag, score_left, score_up])
    score_max = np.max(scores_array)
    index_max = np.argmax(scores_array)
    return score_max, index_max

def compute_multiple_alignment(matrix_align, cluster1_list, cluster2_list):
    """Compute mutliple alignment of two clusters of sequences

    Parameters
    ----------
    matrix_align : the alignement matrix containing the indexes
    from where the max was get
    cluster1_list : the list of sequences from the first cluster
    cluster2_list : the list of sequences from the second cluster

    Returns
    ------
    list : all the alignments of both clusters
    """
    matrix_shape = matrix_align.shape
    current_box_i, current_box_j = matrix_shape[0] - 1, matrix_shape[1] - 1
    cluster1_alignments, cluster2_alignments = ["" for i in range(len(cluster1_list))], ["" for i in range(len(cluster2_list))]
    while (current_box_i >= 1) or (current_box_j >= 1) :
        match matrix_align[current_box_i, current_box_j]:
            # diag
            case 0:
                for i in range(len(cluster1_alignments)):
                    cluster1_alignments[i] = cluster1_list[i][current_box_i-1] + cluster1_alignments[i]
                for i in range(len(cluster2_alignments)):
                    cluster2_alignments[i] = cluster2_list[i][current_box_j-1] + cluster2_alignments[i]
                current_box_i -= 1
                current_box_j -= 1
            #left
            case 1:
                for i in range(len(cluster1_alignments)):
                    cluster1_alignments[i] = "-" + cluster1_alignments[i]
                for i in range(len(cluster2_alignments)):
                    cluster2_alignments[i] = cluster2_list[i][current_box_j-1]+cluster2_alignments[i]
                current_box_j -= 1
            #up
            case 2:
                for i in range(len(cluster1_alignments)):
                    cluster1_alignments[i] = cluster1_list[i][current_box_i-1]+cluster1_alignments[i]
                for i in range(len(cluster2_alignments)):
                    cluster2_alignments[i] = "-"+cluster2_alignments[i]
                current_box_i -= 1
            #if we are in gap boxes
            case _:
                if current_box_i == 0:
                    for i in range(len(cluster1_alignments)):
                        cluster1_alignments[i] = "-"+cluster1_alignments[i]
                    for i in range(len(cluster2_alignments)):
                        cluster2_alignments[i] = cluster2_list[i][current_box_j-1]+cluster2_alignments[i]
                    current_box_j -= 1
                else:
                    for i in range(len(cluster1_alignments)):
                        cluster1_alignments[i] = cluster1_list[i][current_box_i-1]+cluster1_alignments[i]
                    for i in range(len(cluster2_alignments)):
                        cluster2_alignments[i] = "-"+cluster2_alignments[i]
                    current_box_i -= 1
    return cluster1_alignments + cluster2_alignments

def compute_cluster_alignment(group, dict_multiple_alignment):
    """Compute two clusters alignement

    Parameters
    ----------
    matrix_align : the alignement matrix containing the indexes
    from where the max was get
    cluster1_list : the list of sequences from the first cluster
    cluster2_list : the list of sequences from the second cluster

    Returns
    ------
    list : all the alignments of both clusters
    """
    group_0 = eval(group[0]) if "(" in group[0] else group[0]
    group_1 = eval(group[1]) if "(" in group[1] else group[1]
    cluster1_seq_list, cluster2_seq_list = dict_multiple_alignment[group_0], dict_multiple_alignment[group_1]
    seq1, seq2 = cluster1_seq_list[0], cluster2_seq_list[0]
    matrix = needleman_wunsch_method.init_matrix(seq1, seq2)
    nlin, ncol = matrix.shape
    matrix_align = np.full(shape=(len(seq1)+1, len(seq2)+1), fill_value=-1, dtype=int)
    for i in range(1, nlin):
        for j in range(1, ncol):
            score_align, index_max = calculate_score_multiple_alignment(matrix, i, j, cluster1_seq_list, cluster2_seq_list)
            matrix[i, j] = score_align
            matrix_align[i, j] = index_max
    return compute_multiple_alignment(matrix_align, cluster1_seq_list, cluster2_seq_list)

