#!/usr/bin/env python3

## @file sparse.py
#  @brief Phylogenetic covariance matrix sparsification.
#  @author Sean Svihla
#

import sys, re
import numpy as np
import scipy.sparse


NEWICK_PATTERN = r"([,();])([^:;,()\s]*)(?:\s*:\s*([\d.e-]+)\s*)?"


## @brief Parse newick format of a tree and return sparse wavelet matrices.
#
#  This function parses a newick string representing a phylogenetic tree. It
#  is assumed that each node has an associated edge length. The sparse
#  matrices returned can be used to construct the sparse covariance matrix
#  of the tree, namely, C = Q.transpose()@B.
#
#  Parameters:
#    @param list  tokens  A list containing re.Match objects, parsed from
#                             the newick string according to NEWICK_PATTERN.
#    
#  Outputs:
#    @return scipy.sparse.csc_matrix Q
#                         A sparse matrix containing the Haar-like wavelet
#                             basis for the tree. Leaves are indexed by row
#                             and wavelets are indexed by column.
#    @return scipy.sparse.csc_matrix B
#                         The same as Q, except that the wavelet values include
#                             a factor of the trace branch length from each
#                             to the internal node associated with the wavelet.
#
def _parse_newick(newick):
    tree_size = newick.count(')') + 1
    tokens = re.finditer(NEWICK_PATTERN, newick)

    # Pre-allocate some memeory to save time. We first allocate the average
    # amount used. When more is needed, we increment by a standard devaition.
    # These values were found experimentally.
    N = 9.78                               
    M = 0.47                              
    shape = (int(N*tree_size*np.log(tree_size)),)
    pad = (0,int(M*tree_size*np.log(tree_size)))

    rows  = np.zeros(shape)
    cols  = np.zeros(shape)
    valsQ = np.zeros(shape)
    valsB = np.zeros(shape)

    tbl   = []
    stack = []
    indxs = []
    leaf_indx = 0
    indx = 0                                # index of COO matrix entry
    col  = 0
    for token in tokens:
        delim, name, length = token.groups()
        if delim == '(':
            stack.append(len(indxs))
        if length is None:                  # regex picks up stray commas,
            continue                        # just ignore them
        length = float(length)
        if delim == ')':                    # interior node
            stop = stack.pop()              # construct daughter wavelets
            a, b = indxs.pop(stop)
            while len(indxs) > stop:
                _, c = indxs.pop(stop)
                
                L0 = b - a
                L1 = c - b
                L  = L0 + L1
               
                # We use rows and cols as a "test case" for memory allocation.
                # Checking this using try/catch vs. if/else is faster since we
                # do not expect many hits.
                try:
                    rows[indx:indx+L] = np.arange(a,c)
                    cols[indx:indx+L] = np.full((L,), col, dtype=int)
                except IndexError:
                    # allocate additional memory as needed
                    np.pad(rows, pad)
                    np.pad(cols, pad)
                    np.pad(valsQ, pad)
                    np.pad(valsB, pad)
                    rows[indx:indx+L] = np.arange(a,c)
                    cols[indx:indx+L] = np.full((L,), col, dtype=int)

                # left subtree
                lval = np.sqrt(L1 / L0 / L)
                lvalvec = np.full((L0,), lval, dtype=float)
                valsQ[indx:indx+L0] = lvalvec
                valsB[indx:indx+L0] = lvalvec * tbl[a:b]

                # right subtree
                rval = -np.sqrt(L0 / L1 / L)
                rvalvec = np.full((L1,), rval, dtype=float)
                valsQ[indx+L0:indx+L] = rvalvec
                valsB[indx+L0:indx+L] = rvalvec * tbl[b:c]

                indx += L
                col += 1
                b = c
            for i in np.arange(a,b):
                tbl[i] += length * L
            indxs.append((a, b))
        else:                                         # leaf node
            indxs.append((leaf_indx, leaf_indx+1))    # store relevant info
            tbl.append(length)
            leaf_indx += 1

    # rows and cols 
    try:
        rows[indx:indx+L]  = np.arange(L)
        cols[indx:indx+L]  = np.full((L,), col, dtype=int)
    except IndexError:
        # allocate more memory as needed
        np.pad(rows, pad)
        np.pad(cols, pad)
        np.pad(valsQ, pad)
        np.pad(valsB, pad)
        rows[indx:indx+L]  = np.arange(L)
        cols[indx:indx+L]  = np.full((L,), col, dtype=int)

    # mother wavelet
    valvec = np.full((L,), 1/np.sqrt(L))
    valsQ[indx:indx+L] = valvec
    valsB[indx:indx+L] = valvec * tbl
    indx += L

    # might have over-allocated, so just take the first 'indx' entries
    Q = scipy.sparse.csc_array((valsQ[:indx], (rows[:indx], cols[:indx])),
                               shape=(L,L))
    B = scipy.sparse.csc_array((valsB[:indx], (rows[:indx], cols[:indx])),
                               shape=(L,L))
    return Q, B


def sparsify(string, filename=False):
    if filename:
        try:
            with open(string) as f:
                string = f.read()
        except FileNotFoundError:
            sys.exit(f"Cannot find the file {string}!")
    
    Q, B = _parse_newick(string)
    return Q.transpose()@B, Q
