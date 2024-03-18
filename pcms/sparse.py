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
#    @return scipy.sparse.csr_matrix Q
#                         A sparse matrix containing the Haar-like wavelet
#                             basis for the tree. Leaves are indexed by row
#                             and wavelets are indexed by column.
#    @return scipy.sparse.csr_matrix B
#                         The same as Q, except that the wavelet values include
#                             a factor of the trace branch length from each
#                             to the internal node associated with the wavelet.
#
def _parse_newick(tokens):
    rows  = []
    cols  = []
    valsQ = []
    valsB = []

    tbl   = []
    stack = []
    indxs = []
    indx = 0
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
                
                # left subtree
                lval = np.sqrt(L1 / L0 / L)
                lvalvec = np.full((L0,), lval, dtype=float)
                valsQ.extend(lvalvec)
                valsB.extend(lvalvec*tbl[a:b])

                # right subtree
                rval = -np.sqrt(L0 / L1 / L)
                rvalvec = np.full((L1,), rval, dtype=float)
                valsQ.extend(rvalvec)
                valsB.extend(rvalvec*tbl[b:c])
               
                rows.extend(np.arange(a, c))                        
                cols.extend(np.full((L,), col, dtype=int))
                
                col += 1
                b = c
            for i in np.arange(a, b):
                tbl[i] += length * L
            indxs.append((a, b))
        else:                               # leaf node
            indxs.append((indx, indx+1))    # store relevant info
            tbl.append(length)
            indx += 1

    # mother wavelet
    valvec = np.full((indx,), 1/np.sqrt(indx))
    valsQ.extend(valvec)
    valsB.extend(valvec*tbl)
    rows.extend(np.arange(indx))
    cols.extend(np.full((indx,), col, dtype=int))

    Q = scipy.sparse.csr_array((valsQ,(rows,cols)),shape=(indx,indx))
    B = scipy.sparse.csr_array((valsB,(rows,cols)),shape=(indx,indx))
    return Q, B


def sparsify(string, filename=False):
    if filename:
        try:
            with open(string) as f:
                tokens = re.finditer(NEWICK_PATTERN, f.read())
        except FileNotFoundError:
            sys.exit(f"Cannot find the file {string}!")
    else:
        tokens = re.finditer(NEWICK_PATTERN, string)
    Q, B = _parse_newick(tokens)
    return Q.transpose()@B, Q
