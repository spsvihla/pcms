import scipy

import pcms._haar
import pcms.tree


def sparsify(input: pcms.tree.Tree | str) -> scipy.sparse.csc_matrix:
    """
    Convert a tree into a sparse Haar wavelet basis and covariance matrix.

    Parameters
    ----------
    input : pcms.tree.Tree or str
        A `pcms.tree.Tree` object, or a string path to a Newick-formatted file 
        that can be parsed into a tree.

    Returns
    -------
    (scipy.sparse.csc_matrix, scipy.sparse.csc_matrix)
        A tuple containing:
        - The Haar wavelet basis matrix in sparse CSC format.
        - The corresponding covariance matrix in sparse CSC format.

    Raises
    ------
    ValueError
        If `input` is not a `pcms.tree.Tree` instance or a string.
    """
    if isinstance(input, pcms.tree.Tree):
        basis, cov = pcms._haar.sparsify(input)
        return scipy.sparse.csc_matrix(basis), scipy.sparse.csc_matrix(cov)
    elif isinstance(input, str):
        tree = pcms.tree.nwk2tree(input)
        basis, cov = pcms._haar.sparsify(tree)
        return scipy.sparse.csc_matrix(basis), scipy.sparse.csc_matrix(cov)
    else:
        raise ValueError("Expected pcms.tree.Tree or string (filename)")