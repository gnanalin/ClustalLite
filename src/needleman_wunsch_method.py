"""Program which aims at reproducing the algorithm of Needleman-Wunsch.

Usage:
------
------
python needleman_wunsch_method.py
"""

__authors__ = "Stéphanie Gnanalingam"
__contact__ = "stephanie.gnanalingam@etu.u-paris.fr"
__date__ = "2024-09-06"

import pandas as pd
import numpy as np

BLOSUM_MATRIX = pd.read_csv("data/BLOSUM_MATRIX.csv", index_col=[0])
GAP = 8


def init_matrix(seq1, seq2):
    """Initialize the matrix with the correct shape and gap values.

    Parameters
    ----------
    seq1 : the first sequence
    seq2 : the second sequence

    Returns
    -------
    matrix : np.array
    """
    seq1_len, seq2_len = len(seq1), len(seq2)
    matrix = np.full(
        shape=(seq1_len + 1, seq2_len + 1), fill_value=0, dtype=int
        )
    for i in range(1, seq2_len + 1):
        matrix[0][i] = GAP * (-i)
    for i in range(1, seq1_len + 1):
        matrix[i][0] = GAP * (-i)
    return matrix


def calculate_score(matrix, current_line, current_col, current_line_aa, current_col_aa):
    """Calculate the score for matrix[current_line, current_col].

    Parameters
    ----------
    matrix : a score matrix
    current_line, current_col : the coordinates of the current box
    current_line_AA, current_col_AA : the amino acids of the current box

    Returns
    -------
    tuple : (score, index) where index will be the box which gave
    the maximal score (0: diagonal, 1: left, 2: up). If there are equals with
    the diagonal, the diagonal box will be selected
    and otherwise, it left box is selected
    """
    score_diag = (
        matrix[current_line - 1, current_col - 1]
        + BLOSUM_MATRIX.loc[current_line_aa, current_col_aa]
    )
    score_up = matrix[current_line - 1, current_col] - GAP
    score_left = matrix[current_line, current_col - 1] - GAP
    scores_array = np.array([score_diag, score_left, score_up])
    score_max = np.max(scores_array)
    index_max = np.argmax(scores_array)
    return (score_max, index_max)


def complete_matrice(matrix, seq1, seq2):
    """Calculate the score for all boxes.

    Parameters
    ----------
    matrix : a score matrix
    seq1 : the first sequence
    seq2 : the second sequence

    Returns
    -------
    score_max_align : the alignement score
    matrix_align : the matrix containing the indexes that gave the maximal
    score for each boxes
    """
    nlin, ncol = matrix.shape
    matrix_align = np.full(
        shape=(len(seq1) + 1, len(seq2) + 1), fill_value=-1, dtype=int
    )
    for i in range(1, nlin):
        for j in range(1, ncol):
            score_align, index_max = calculate_score(
                matrix, i, j, seq1[i - 1], seq2[j - 1]
            )
            matrix[i, j] = score_align
            matrix_align[i, j] = index_max
    return matrix[(nlin - 1), (ncol - 1)], matrix_align


def compute_alignement(matrix_align, seq1, seq2):
    """Compute alignment between two sequences.

    Parameters
    ----------
    matrix_align : an alignment matrix
    seq1 : the first sequence
    seq2 : the second sequence

    Returns
    -------
    tuple : (aligned seq 1, aligned seq 2)
    """
    matrix_shape = matrix_align.shape
    current_box_i, current_box_j = matrix_shape[0] - 1, matrix_shape[1] - 1
    alignement_seq1, alignement_seq2 = "", ""
    while (current_box_i >= 1) or (current_box_j >= 1):
        match matrix_align[current_box_i, current_box_j]:
            # diag
            case 0:
                alignement_seq1 = seq1[current_box_i - 1] + alignement_seq1
                alignement_seq2 = seq2[current_box_j - 1] + alignement_seq2
                current_box_i -= 1
                current_box_j -= 1
            # left
            case 1:
                alignement_seq1 = "-" + alignement_seq1
                alignement_seq2 = seq2[current_box_j - 1] + alignement_seq2
                current_box_j -= 1
            # up
            case 2:
                alignement_seq2 = "-" + alignement_seq2
                alignement_seq1 = seq1[current_box_i - 1] + alignement_seq1
                current_box_i -= 1
            # if we are in gap boxes
            case _:
                if current_box_i == 0:
                    alignement_seq1 = "-" + alignement_seq1
                    alignement_seq2 = seq2[current_box_j - 1] + alignement_seq2
                    current_box_j -= 1
                else:
                    alignement_seq2 = "-" + alignement_seq2
                    alignement_seq1 = seq1[current_box_i - 1] + alignement_seq1
                    current_box_i -= 1
    return (alignement_seq1, alignement_seq2)


if __name__ == "__main__":
    SEQ1 = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"
    SEQ2 = "MTAMEESQSDISLELPLSQETFSGLWKLLPPEDILPSPHCMDDLLLPQDVEEFFEGPSEALRVSGAPAAQDPVTETPGPVAPAPATPWPLSSFVPSQKTYQGNYGFHLGFLQSGTAKSVMCTYSPPLNKLFCQLAKTCPVQLWVSATPPAGSRVRAMAIYKKSQHMTEVVRRCPHHERCSDGDGLAPPQHLIRVEGNLYPEYLEDRQTFRHSVVVPYEPPEAGSEYTTIHYKYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRDSFEVRVCACPGRDRRTEEENFRKKEVLCPELPPGSAKRALPTCTSASPPQKKKPLDGEYFTLKIRGRKRFEMFRELNEALELKDAHATEESGDSRAHSSYLKTKKGQSTSRHKKTMVKKVGPDSD"
    matrix = init_matrix(SEQ1, SEQ2)
    score_max_align, matrix_align = complete_matrice(matrix, SEQ1, SEQ2)
    seq1_align, seq2_align = compute_alignement(matrix_align, SEQ1, SEQ2)
    print("seq1:", seq1_align)
    print("seq2:", seq2_align)
    print("score:", score_max_align)
