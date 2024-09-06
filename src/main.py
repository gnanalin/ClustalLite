"""Program which aims at reproducing the algorithm of Clustal

Usage:
------
------

"""

__authors__ = ("Stéphanie Gnanalingam")
__contact__ = ("stephanie.gnanalingam@etu.u-paris.fr")
__date__ = "2024-09-04"

import pandas as pd
import argparse
import sys, os
from Bio import SeqIO

BLOSUM_MATRIX = pd.read_csv("data/BLOSUM_MATRIX.csv", index_col=[0])


def check_file(fastq_file):
    """Check whether the path exists.

    Parameters
    ----------
    fastq_file : the path of the fastq file

    Retuns
    ------
    bool : True if the path exists
    """
    if os.path.exists(fastq_file):
        print(f"The file {fastq_file} is going to be analysed...")
        return True
    else:
        print("The file does not exists !")
        return False


def parse_argument():
    """Return the path of the fastq files.

    Parameters
    ----------
    -h : for help
    fastq : the path of the fastq file

    Returns
    -------
    str : the path of the fastq file
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq", type=str ,help="pass your fastq file")
    arg= parser.parse_args()
    return arg.fastq

def get_fastq_info(fastq_file):
    """Extract the names sequences from a fastq file

    Parameters
    ----------
    fastq_file : the path of the fastq file

    Returns
    -------
    dict : a dictionnary in which keys are names and values are sequences
    """
    return {record.id : record.seq \
        for record in SeqIO.parse(fastq_file, "fasta")}



if __name__ == "__main__":
    # parsing the argument
    fastq_file = parse_argument()
    # checking the existance of the file
    if check_file(fastq_file):
        print("ok")
    else:
        sys.exit("The file doesn't exist. Please check the path")
    for k,v in get_fastq_info(fastq_file).items():
        print(f"{k} -> {v}\n")