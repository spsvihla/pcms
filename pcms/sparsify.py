#!/usr/bin/env python3

## @file sparse.py
#  @brief Functions for sparsifying hierarchical covariance matrices. 
#  @author Sean Svihla
#

# system imports
import itertools
import numpy as np
from scipy.sparse import csr_array


## @brief Adds the n_leaves and leaf_lst features to all nodes in the tree.
#
#  The n_leaves feature is simply the number of leaves below the node.
#
#  The leaf_lst feature is a list of integers indicating the leaves beneath
#  the node. The leaves are numbered according to postorder traversal.
#
#  @param ete3.TreeNode tree    The tree
#
def add_leaf_features(tree):
    leaf_number = 0                         # postorder should visit leaves
    for node in tree.traverse("postorder"): # from left to right
        node.add_features(col_nums=np.array([], dtype=int))
        if node.is_leaf():
            node.add_features(n_leaves=1, 
                              leaf_lst=np.array([leaf_number], dtype=int))
            leaf_number += 1
            continue
        n_leaves = 0
        leaf_lst = np.ndarray((0,), dtype=int)
        for child in node.children:
            n_leaves += child.n_leaves
            leaf_lst = np.append(leaf_lst, child.leaf_lst)
        node.add_features(n_leaves=n_leaves, leaf_lst=leaf_lst) 


## @brief Computes the Haar-like wavelets associated with the tree.
#
#  @param ete3.TreeNode tree        The tree
#  @return scipy.sparse.csr_array   Array of wavelet values row-indexed by leaf
#                                       and col-indexed by wavelet
#
def compute_haar_like_wavelets(tree):
    rows = np.ndarray((0,))
    cols = np.ndarray((0,))
    vals = np.ndarray((0,))
    colno = 0
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            continue

        l_supp = node.children[0].leaf_lst   # left support
        l_size = node.children[0].n_leaves

        # handle the root separately
        if node.is_root():
            val = 1/np.sqrt(l_size)
            rows = np.append(rows, l_supp)
            cols = np.append(cols, np.full((l_size,), colno, dtype=int))
            vals = np.append(vals, np.full((l_size,), l_val))
            continue

        for child in node.children[1:]:
            r_supp = child.leaf_lst          # right support
            r_size = child.n_leaves
            n_leaves  = l_size + r_size      # whole support
            
            # handle left sub-tree
            l_val =  np.sqrt(r_size / l_size / n_leaves)
            rows = np.append(rows, l_supp)
            cols = np.append(cols, np.full((l_size,), colno, dtype=int))
            vals = np.append(vals, np.full((l_size,), l_val))
            
            # handle right sub-tree
            r_val = -np.sqrt(l_size / r_size / n_leaves)
            rows = np.append(rows, r_supp)
            cols = np.append(cols, np.full((r_size,), colno, dtype=int))
            vals = np.append(vals, np.full((r_size,), r_val))
           
            # store the column number
            node.col_nums = np.append(node.col_nums, colno)

            # increment values
            l_supp = np.append(l_supp, r_supp)
            l_size = n_leaves
            colno += 1

    return csr_array((vals,(rows,cols)), shape=(tree.n_leaves, tree.n_leaves))


## Helper function for trace_branch_length
#
_TBL_CACHE = dict()
def _tbl(u, v):
    uid = u.get_topology_id()
    vid = v.get_topology_id()
    if (uid, vid) in _TBL_CACHE:
        return _TBL_CACHE[(uid, vid)]
    if uid == vid:
        return 0
    result = u.n_leaves * u.dist + _tbl(u.up, v)
    _TBL_CACHE[(uid, vid)] = result
    return result


## @brief Computes the trace branch length between nodes.
#
#  The trace branch length between two nodes is the sum of the lengths
#  of edges connecting them, each weighted by the number of leaves which
#  descend from that edge. 
#
#  @param  ete3.TreeNode u, v    The nodes of interest
#  @return float                 Trace branch length between 'u' and 'v'
#
def trace_branch_length(u, v):
    if v in u.get_ancestors():
        return _tbl(u, v)
    if u in v.get_ancestors():
        return _tbl(v, u)
    w = u.get_common_ancestor(v)
    return _tbl(u, w) + _tbl(v, w)


## @brief Computes sparse form of the covariance matrix associated with tree.
#
#  @param ete3.TreeNode tree  The tree.
#  @return scipy.sparse.csr_matrix  The sparse covariance matrix.
#  @return scipy.sparse.csr_matrix  The sparse basis matrix.
#
def sparsify(tree):
    add_leaf_features(tree)
    phi = compute_haar_like_wavelets(tree)
    leaves = tree.get_leaves()              # leaves in postorder

    n_leaves = len(leaves)
    counter = 0

    rows = np.ndarray((0,))
    cols = np.ndarray((0,))
    vals = np.ndarray((0,))
    for node_i, node_j in itertools.product(tree.traverse("postorder"), 
                                            tree.traverse("postorder")):
        intersect = np.intersect1d(node_i.leaf_lst, node_j.leaf_lst)
        if not intersect.size > 0:
            continue 
        for i, j in itertools.product(node_i.col_nums, node_j.col_nums):
            print(f"{counter}")
            counter += 1

            rows = np.append(rows, i)
            cols = np.append(cols, j)
            val = 0
            for leaf_index in intersect:
                leaf = leaves[leaf_index]
                tbl = trace_branch_length(leaf, node_i)
                val += tbl * phi[leaf_index, i] * phi[leaf_index, j]
            vals = np.append(vals, val)

    return csr_array((vals,(rows,cols)),shape=(tree.n_leaves,tree.n_leaves)),\
           phi
