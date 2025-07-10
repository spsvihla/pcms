"""
haar.py

A wrapper for the pcms._haar package.
"""

from typing import Union, Sequence, Optional

import numpy as np
from scipy.special import factorial
from scipy.stats import beta
from scipy.sparse import csc_matrix

import pcms._haar 
import pcms.tree
from pcms.tree import CriticalBetaSplittingDistribution 


def cdf_rand_basis(
    ys: Union[float, Sequence[float], np.ndarray],
    func: Union[Sequence[float], np.ndarray],
    eps: float = 0.001,
    delta: float = 0.001,
    seed: Optional[int] = None
) -> Union[float, np.ndarray]:
    """
    Evaluate the CDF of ⟨f, φ⟩ for a Haar-like wavelet φ on a critical beta-splitting tree.

    φ is associated with an interior node v of a critical beta-splitting tree T,
    and f is a non-random function on the leaves of the subtree T(v). The split 
    below φ is sampled randomly according to the critical beta-splitting 
    distribution q(n, i).

    Parameters
    ----------
    ys: float or array-like
        Points at which to evaluate the CDF.
    func: array-like
        Function values on the leaves of T(v).
    eps: float (optional, default 0.001)
        Desired precision of the estimate.
    delta: float (optional, default 0.001)
        Probability that the estimate lies within the desired precision.
    seed: int (optional, default None)
        Random seed.

    Returns
    -------
    float or np.ndarray
        Estimate(s) of the CDF at the given point(s).
    """
    ys_arr = np.atleast_1d(ys).astype(np.float64)
    func_arr = np.asarray(func, dtype=np.float64)
    n = func_arr.size
    num_iter = min(factorial(n), int(np.ceil(1.0 / (4 * eps**2 * delta))))
    pmf = CriticalBetaSplittingDistribution(n).pmf
    seed = int(seed) if seed is not None else None
    result = pcms._haar.cdf_cbst_topology(ys_arr, func_arr, pmf, int(num_iter), seed)
    return result.item() if np.isscalar(ys) else result


def cdf_rand_func(
    ys: Union[float, Sequence[float], np.ndarray],
    split_size: int,
    subtree_size: int
) -> float:
    """
    Evaluate the CDF of ⟨f, φ⟩ where φ is a Haar-like wavelet over a fixed split.

    φ is associated with an interior node that has `subtree_size` descendant leaves
    and a split of size `split_size`. In this case, the distribution of ⟨f, φ⟩ 
    is an affine transformation of a Beta(split_size, subtree_size - split_size) distribution.

    Parameters
    ----------
    ys: float or array-like
        Points at which to evaluate the CDF.
    split_size: int
        Size of one side of the split.
    subtree_size: int
        Total number of leaves in the subtree.

    Returns
    -------
    float or np.ndarray
        Value(s) of the CDF at the specified point(s).
    """
    ys = np.asarray(ys)
    a = np.sqrt((subtree_size - split_size) / split_size / subtree_size)
    b = np.sqrt(split_size / (subtree_size - split_size) / subtree_size)
    return beta.cdf((ys + b) / (a + b), split_size, subtree_size - split_size)


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
    if isinstance(input, pcms.tree.Tree):
        basis, cov = pcms._haar.sparsify(input._tree)
        return csc_matrix(basis), csc_matrix(cov)
    elif isinstance(input, str):
        tree = pcms.tree.nwk2tree(input)
        basis, cov = pcms._haar.sparsify(tree._tree)
        return csc_matrix(basis), csc_matrix(cov)
    else:
        raise ValueError("Expected pcms.tree.Tree or string (filename)")