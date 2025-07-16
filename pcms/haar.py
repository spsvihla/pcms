"""
haar.py

A wrapper for the pcms._haar package.
"""

from typing import Union, Sequence, Optional

import numpy as np
from scipy.special import factorial
from scipy.stats import dirichlet
from scipy.sparse import csc_matrix

import pcms._haar 
import pcms.tree
from pcms.tree import CriticalBetaSplittingDistribution 


def cdf_proj_cbst(
    ys: Union[float, Sequence[float], np.ndarray],
    func: Union[Sequence[float], np.ndarray],
    eps: float = 0.01,
    delta: float = 0.01,
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
    eps: float (optional, default 0.01)
        Desired precision of the estimate.
    delta: float (optional, default 0.01)
        Probability that the estimate lies within the desired precision.
    seed: int (optional, default None)
        Random seed.

    Returns
    -------
    float or np.ndarray
        Estimate(s) of the CDF at the given point(s).
    """
    ys = np.atleast_1d(ys).astype(np.float64)
    func = np.array(func, dtype=np.float64, copy=True)  # copy since f is permuted in place
    n = func.size
    num_iter = int(min(factorial(n), int(np.ceil(1.0 / (4 * eps**2 * delta)))))
    pmf = CriticalBetaSplittingDistribution(n).pmf
    seed = int(seed) if seed is not None else None
    result = pcms._haar.cdf_proj_cbst(ys, func, pmf, num_iter, seed)
    result = np.asarray(result)
    return result.item() if result.ndim == 0 or result.size == 1 else result


def cdf_proj_dirichlet_diff(
    ys: Union[float, Sequence[float], np.ndarray],
    tree: pcms.tree.Tree,
    node: int,
    eps: float = 0.01,
    delta: float = 0.01
) -> float:
    """
    Evaluate the CDF of ⟨f, φ⟩ where φ is a Haar-like wavelet over a fixed split.
    The function `f` is considered to be a difference between two dirichlet samples.

    Parameters
    ----------
    ys: float or array-like
        Points at which to evaluate the CDF.
    tree: pcms.tree.Tree
        The target tree.
    node: int
        The target node.
    eps: float (optional, default 0.01)
        Desired precision of the estimate.
    delta: float (optional, default 0.01)
        Probability that the estimate lies within the desired precision.

    Returns
    -------
    float or np.ndarray
        Value(s) of the CDF at the specified point(s).
    """
    ys = np.asarray(ys)
    n_leaves = tree.find_n_leaves()
    num_samples = int(np.ceil(1.0 / (4 * eps**2 * delta)))
    subtree_size = tree.get_subtree_size(node)
    split_size = tree.get_subtree_size(tree.get_child(node))
    m = subtree_size - split_size
    a = np.sqrt(m / split_size / subtree_size)
    b = -np.sqrt(split_size / m / subtree_size)
    alpha = np.array([split_size, m, n_leaves - m - split_size])
    rng = np.random.default_rng()
    f = rng.gamma(shape=alpha, size=(num_samples, len(alpha)))
    f /= f.sum(axis=1, keepdims=True)
    g = rng.gamma(shape=alpha, size=(num_samples, len(alpha)))
    g /= g.sum(axis=1, keepdims=True)
    coef = np.sqrt(m * split_size / subtree_size)
    samples = coef * (a * (f[:,0] - g[:,0]) + b * (f[:,1] - g[:,1]))
    samples.sort()
    result = np.searchsorted(samples, ys, side='left') / num_samples
    return result.item() if np.isscalar(ys) else result


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