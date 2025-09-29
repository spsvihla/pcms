"""
@file haar.py
@brief A wrapper for the pcms._haar package.
@author Sean Svihla
"""

from scipy.sparse import csc_matrix

import pcms._haar 
import pcms.tree


def sparsify(input: pcms.tree.Tree | str) -> csc_matrix:
    """
    Convert a tree to sparse Haar wavelet basis and covariance matrices.

    Parameters
    ----------
    input: pcms.tree.Tree or str
        A `Tree` object or a string path to a Newick-formatted file.

    Returns
    -------
    scipy.sparse.csc_matrix
        The Haar basis matrix and corresponding covariance matrix.

    Raises
    ------
    ValueError
        If `input` is neither a `Tree` instance nor a valid filename.
    """
    tree = pcms.tree.nwk2tree(input, ensure_planted=True) if isinstance(input, str) else input
    if not tree.find_is_planted():
        raise ValueError("Tree must be planted.")
    Q, S = pcms._haar.sparsify(tree._tree)
    n_leaves = tree.find_n_leaves()
    n_wavelets = tree.find_n_wavelets()
    Q = csc_matrix(Q, shape=(n_leaves, n_wavelets))
    S = csc_matrix(S, shape=(n_wavelets, n_wavelets))
    Q.eliminate_zeros()
    S.eliminate_zeros()
    return Q, S