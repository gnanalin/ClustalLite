"""Program which aims at reproducing the algorithm of Clustal.

Usage:
------
------
python main.py
"""

__authors__ = "Stéphanie Gnanalingam"
__contact__ = "stephanie.gnanalingam@etu.u-paris.fr"
__date__ = "2024-09-04"

import argparse
import os
import copy
import numpy as np
from Bio import SeqIO
from colorama import Fore, Style
import needleman_wunsch_method
import upgma_method
import multiple_alignment


def check_file(fasta_file):
    """Check whether the path exists.

    Parameters
    ----------
    fastq_file : the path of the fastq file

    Returns
    ------
    bool : True if the path exists
    """
    return os.path.exists(fasta_file)


def parse_argument():
    """Return the path of the fastq files.

    Parameters
    ----------
    -h : for help
    fastq : the path of the fastq file

    Returns
    -------
    arg.fastq : str, the path of the fastq file
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", type=str, help="pass your fasta file")
    arg = parser.parse_args()
    return arg.fasta


def get_fastq_info(fasta_file):
    """Extract the names sequences from a fastq file.

    Parameters
    ----------
    fastq_file : the path of the fastq file

    Returns
    -------
    dict : a dictionnary in which keys are names
    and values are lists with sequence. The sequences
    are in a list so that it can ease the multiple alignment
    step, in which we suppose to have clusters
    """
    return {record.id: [str(record.seq)]
            for record in SeqIO.parse(fasta_file, "fasta")}


def compute_all_needleman_wunch(dict_sequences):
    """Compute needleman wunsch alignement for all the sequences.

    Parameters
    ----------
    dict_sequences : a dictionnary result of get_fastq_info()

    Returns
    -------
    scores_matrix : np.array, the matrix of scores
    dict_two_align : a dictionnary of all the pairwise alignements
    example : {(Seq1, Seq2): [alignment_seq1, alignement_seq2]}
    """
    dict_two_align = {}
    labels = list(dict_sequences.keys())
    number_of_sequences = len(dict_sequences)
    scores_matrix = np.full(
        shape=(number_of_sequences, number_of_sequences), fill_value=-10, dtype=int
    )
    for seq1_i in range(len(labels)):
        for seq2_i in range(seq1_i + 1, len(labels)):
            seq1_name, seq2_name = labels[seq1_i], labels[seq2_i]
            seq1, seq2 = (
                dict_sequences[labels[seq1_i]][0],
                dict_sequences[labels[seq2_i]][0],
            )
            matrix_nw = needleman_wunsch_method.init_matrix(
                seq1, seq2
            )
            score_align, matrix_align = needleman_wunsch_method.complete_matrice(
                matrix_nw, seq1, seq2
            )
            seq1_align, seq2_align = needleman_wunsch_method.compute_alignement(
                matrix_align, seq1, seq2
            )
            scores_matrix[seq1_i, seq2_i] = score_align
            dict_two_align[(seq1_name, seq2_name)] = [seq1_align, seq2_align]
            for i in range(0, len(seq1_align), 60):
                print(
                    f"{Fore.MAGENTA+Style.BRIGHT+seq1_name}: "
                    + f"{Style.RESET_ALL}{seq1_align[i:i+60]}\n"
                    + f"{Fore.MAGENTA+Style.BRIGHT+seq2_name}: "
                    + f"{Style.RESET_ALL}{seq2_align[i:i+60]}\n"
                )
            print(
                f"{Fore.RED+Style.BRIGHT} Alignement score : " +
                f"{str(score_align)}{Style.RESET_ALL}\n"
            )
    return scores_matrix, dict_two_align


def normalize_array(scores_matrix):
    """Normalize the matrix values from 0 to 1 to get distances.

    Parameters
    ----------
    scores_matrix : the matrix resulting from compute_all_needleman_wunch()

    Returns
    -------
    scores_matrix : np.array, the normalized matrix of scores
    """
    nlin = scores_matrix.shape[0]
    minimum = float("inf")
    maximum = -float("inf")
    for i in range(nlin):
        for j in range(i + 1, nlin):
            eij = scores_matrix[i, j]
            if eij < minimum:
                minimum = eij
            if eij > maximum:
                maximum = eij
    scores_matrix = (scores_matrix - minimum) / (
        maximum - minimum
    )
    return scores_matrix


def multiple_alignment_process(dict_upgma, dict_two_align, nb_seq):
    """Compute multiple alignment step.

    Parameters
    ----------
    dict_upgma : the result of from upgma algorithm to follow
    in order to align sequences
    dict_two_align : a dictinnary containing the pairwise
    alignments and the single sequences
    nb_seq : the total number of sequences to align

    Returns
    -------
    dict_multiple_alignment : dict, the alignements
    ordered_seq_names : list, sequences ordered names
    """
    dict_multiple_alignment = copy.deepcopy(dict_two_align)
    for group in dict_upgma.keys():
        cluster1 = eval(group[0]) if "(" in group[0] else group[0]
        cluster2 = eval(group[1]) if "(" in group[1] else group[1]
        if group not in dict_multiple_alignment.keys():
            dict_multiple_alignment[
                group
            ] = multiple_alignment.compute_cluster_alignment(
                group, dict_multiple_alignment
            )
            del dict_multiple_alignment[cluster1]
            del dict_multiple_alignment[cluster2]
    ordered_seq_names = upgma_method.upgma_ordered_seq_name_list(
        nb_seq, [], list(dict_upgma.keys())[-1]
    )
    return dict_multiple_alignment, ordered_seq_names


if __name__ == "__main__":
    # parsing the argument
    fasta_file = parse_argument()
    # checking the existance of the file
    if check_file(fasta_file):
        print(
            Fore.BLUE
            + Style.BRIGHT
            + "The file exists and is going to be parsed...\n"
            + Style.RESET_ALL
        )
    else:
        raise FileNotFoundError
    print(
        Fore.BLUE
        + Style.BRIGHT
        + "Here is the file's content :\n"
        + Style.RESET_ALL)
    # getting all the sequences
    all_fasta = get_fastq_info(fasta_file)
    for k, v in all_fasta.items():
        print(
            f"{Fore.MAGENTA + Style.BRIGHT + k + Style.RESET_ALL}"
            + f": {v[0]}\n"
        )
    # computing all the pairwise alignements
    print(
        Fore.BLUE
        + Style.BRIGHT
        + "Here are the Needleman-Wunsch alignements :\n"
        + Style.RESET_ALL
    )
    scores_matrix, dict_two_align = compute_all_needleman_wunch(all_fasta)
    # normalize the array to get distances instead of scores
    distance_matrix = normalize_array(scores_matrix)
    # computing UPGMA
    labels = list(all_fasta.keys())
    clusters = {label: 1 for label in labels}
    print(
        Fore.BLUE
        + Style.BRIGHT
        + "Here is the UPGMA algorithm result :\n"
        + Style.RESET_ALL
    )
    dict_upgma = upgma_method.upgma({}, labels, clusters, distance_matrix)
    dict_two_align.update(all_fasta)
    # computing multiple alignement process
    print(
        Fore.BLUE
        + Style.BRIGHT
        + "Here is the multiple alignment result :\n"
        + Style.RESET_ALL
    )
    dict_multiple_alignment, ordered_seq_names = multiple_alignment_process(
        dict_upgma, dict_two_align, len(labels)
    )
    final_alignments = dict_multiple_alignment[list(dict_multiple_alignment.keys())[-1]]
    for index_seq in range(0, len(final_alignments[0]), 40):
        for seq_name, alignement in zip(ordered_seq_names, final_alignments):
            print(
                f"{Fore.GREEN+Style.BRIGHT+seq_name+Style.RESET_ALL}"
                + f": {Style.BRIGHT+alignement[index_seq:index_seq+40]}"
            )
        print("\n")
