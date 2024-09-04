"""Programm which aims at reproducing the algorithm of UPGMA

Usage:
------
------

"""

import numpy as np

def find_lowest_coordinates(distance_matrix):
    """Find the lowest value in the matrix
    
    Parameters
    ----------
    distance_matrix : a distance matrix 
    
    Returns
    -------
    tuple : (vertical coordinate, horizontal coordinate)
    """
    min_matrix_lin = distance_matrix.min(axis = 1)
    min_matrix_col = distance_matrix.min(axis = 0)
    min_lin = min_matrix_lin.argmin(axis = 0)
    min_col = min_matrix_col.argmin(axis = 0)
    return(min_lin, min_col)

def calculate_distance(ni, dki, nj, dkj):
    """Calculate new distance between the clusters (i, j) with the minimum distance and other sequences (k)
    Parameters
    ----------
    ni : number of elements in cluster i
    dki : distance between k and i
    nj : number of elements in cluster j
    dkj : distance between k and j

    Returns
    ------
    list : return the order in which to align the sequences
    """
    return (ni*dki+nj*dkj)/ni+nj

def upgma(order_alignment_list, clusters, distance_matrix):
    """Run the upgma algorithm.

    Parameters
    ----------
    distance_matrix : a numpy array containing the distances between
    the sequences

    Returns
    ------
    list : return the order in which to align the sequences
    """
    labels = list(clusters.keys())
    if len(labels) > 1:
        lowest_coord_lin,  lowest_coord_col= find_lowest_coordinates(distance_matrix)
        seq1, seq2 = labels[lowest_coord_lin], labels[lowest_coord_col]
        if "(" not in seq1:
            order_alignment_list.append(seq1)
        if "(" not in seq2:
            order_alignment_list.append(seq2)
        new_label = "("+seq1+","+seq2+")"
        clusters[new_label] = clusters[lowest_coord_lin]+clusters[lowest_coord_col]
        
        
        
        labels.remove(seq1)
        labels.remove(seq2)
        labels.append(new_label)
        del clusters[seq1]
        del clusters[seq2]