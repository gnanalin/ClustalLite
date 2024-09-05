"""Programm which aims at reproducing the algorithm of UPGMA

Usage:
------
------
python UPGMA.py
"""

import numpy as np

def find_lowest_coordinates(distance_matrix):
    """Find the lowest value in the superior triangle of the matrix
    
    Parameters
    ----------
    distance_matrix : a distance matrix 
    
    Returns
    -------
    tuple : (vertical coordinate, horizontal coordinate)
    """
    nlin = distance_matrix.shape[0]
    minimum = float('inf')
    minimum_coord = (-1, -1)
    for i in range(nlin):
        for j in range(i+1,nlin):
            eij = distance_matrix[i, j]
            if eij < minimum:
                minimum = eij
                minimum_coord = (i, j)
    return minimum_coord

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
    return round(((ni * dki + nj * dkj) / (ni + nj)), 4)

def upgma(order_alignment_dico, labels, clusters, distance_matrix):
    """Run the upgma algorithm.

    Parameters
    ----------
    order_alignement_dico : a dictionnary which will be completed recursively
    where the keys are the two elements with the minimal score and values
    are the distances 
    labels : stores the sequences names
    clusters : dictionnary where the keys are tuples (or tuples of tuples)
    with sequences and the values are the number of sequences in them
    distance_matrix : a numpy array containing the distances between
    the sequences

    Returns
    ------
    list : return the order in which to align the sequences
    """
    # labels first contains the names of the sequences
    # in each recursion, it will assemble two sequences (seq1, seq2) in a str 
    # and add them and will delete seq1 and seq2 from the list.
    if len(labels) > 1:
        lowest_coord_lin,  lowest_coord_col = find_lowest_coordinates(distance_matrix)
        # we get the sequences which have the lowest distance
        seq1, seq2 = labels[lowest_coord_lin], labels[lowest_coord_col]
        new_label = (seq1,seq2)
        # the dictionnary will store the distance between seq1 and seq2
        order_alignment_dico[new_label] = distance_matrix[lowest_coord_lin, lowest_coord_col]
        # we have to add the number of sequences of both sequences for the new sequence name
        clusters[str(new_label)] = clusters[labels[lowest_coord_lin]]+clusters[labels[lowest_coord_col]]
        
        # in the following, the steps for modifying the distance matrix
        nlin = distance_matrix.shape[0]
        # a list which will store the new distances of the other sequences with our new assembly
        new_distance_line_l = []
        for k in range(nlin):
            if k != lowest_coord_lin and k != lowest_coord_col:
                    # from line 87 to line 85, we make sure that we get the distances from the upper triangle
                    if k > lowest_coord_lin:
                        dik = distance_matrix[lowest_coord_lin, k]
                    else:
                        dik = distance_matrix[k, lowest_coord_lin]
                    if k > lowest_coord_col:
                        dkj = distance_matrix[lowest_coord_col, k]
                    else:
                        dkj = distance_matrix[k, lowest_coord_col]
                    new_distance_line_l.append(calculate_distance(clusters[labels[lowest_coord_lin]],
                                                            dik, 
                                                            clusters[labels[lowest_coord_col]],
                                                            dkj
                                                            ))

        # delete the rown and columns corresponding to seq1 and seq2
        distance_matrix = np.delete(distance_matrix, [lowest_coord_lin, lowest_coord_col], axis = 1)
        distance_matrix = np.delete(distance_matrix, [lowest_coord_lin, lowest_coord_col], axis = 0)
        
        # in the following, we modify our matrix to add the new distances 
        # line 107 : a new column that we will add after adding the distances to get the 
        # correct format for our distance matrix
        new_distance_line_c = np.full(nlin-1, -10).reshape(-1, 1)
        distance_matrix = np.vstack(([new_distance_line_l], distance_matrix))
        distance_matrix = np.hstack((new_distance_line_c, distance_matrix))
        
        # in the labels, we don't keep our used sequences (seq1, seq2)
        # and then we add our new label
        labels.insert(0,str(new_label))
        labels.remove(seq1)
        labels.remove(seq2)
        del clusters[seq1]
        del clusters[seq2]
        
        # here is the recursion part
        return upgma(order_alignment_dico, labels, clusters, distance_matrix)
    else:
        return order_alignment_dico

if __name__ == "__main__":
    labels = ["Hu", "Ch", "Go", "Or", "Gi"]
    clusters = {label:1 for label in labels}
    distances = np.array([[-1, 15, 45, 143, 198], 
            [-10, -10, 30, 126, 179],
            [-10, -10, -10, 92, 179],
            [-10, -10, -10, -10, 179],
            [-10, -10, -10, -10, -10]])
    print(upgma({}, labels, clusters, distances))
    
    labels = ["Bsu", "Bst", "Lvi", "Amo", "Mlu"]
    clusters = {label:1 for label in labels}
    distances = np.array([[-10, 0.1715, 0.2147, 0.3091, 0.2326], 
            [-10, -10, 0.2991, 0.3399, 0.2058],
            [-10, -10, -10, 0.2795, 0.3943],
            [-10, -10, -10, -10, 0.4289],
            [-10, -10, -10, -10, -10]])
    print(upgma({}, labels, clusters, distances))
    