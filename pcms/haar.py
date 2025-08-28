"""
haar.py

A wrapper for the pcms._haar package.
"""

from typing import Union, Optional
import random

import numpy as np
from numpy.typing import ArrayLike
from scipy.sparse import csc_matrix

import pcms._haar 
import pcms.tree


def generate_seeds(seed: int, n: int) -> list[int]:
    rng = random.Random(seed)
    return [rng.randint(0, 2**32 - 1) for _ in range(n)]


def cdf_proj_cbst(
    ys: ArrayLike,
    f: ArrayLike,
    eps: float = 0.005,
    delta: float = 0.001,
    n_samples: Optional[int] = None,
    batch_size: int = 100,
    seed: Optional[int] = None
) -> Union[float, np.ndarray]:
    """
    Evaluate the CDF of (λ_v * Δ_v^2) for a Haar-like wavelet φ on a critical 
    beta-splitting tree.

    φ is associated with an interior node v of a critical beta-splitting tree T,
    and f is a non-random function on the leaves of the subtree T(v). The split 
    below φ is sampled randomly according to the critical beta-splitting 
    distribution q(n, i).

    Parameters
    ----------
    ys: array-like
        Points at which to evaluate the CDF.
    f: array-like
        Function values on the leaves of the tree.
    eps: float (optional, default 0.005)
        Desired precision of the estimate.
    delta: float (optional, default 0.001)
        Probability that the estimate lies within the desired precision.
    n_samples: int (optional, default None)
        Number of samples to draw for the estimate (overrides eps and delta).
    seed: int (optional, default None)
        Random seed.

    Returns
    -------
    float or np.ndarray
        Estimate(s) of the CDF at the given point(s).

    Raises
    -------
    ValueError
        If n_samples < 1.
    """
    # parse inputs
    ys = np.atleast_1d(ys).astype(np.float64)
    f = np.array(f, dtype=np.float64)
    seed = int(seed) if seed is not None else None
    n_samples = int(n_samples) if n_samples is not None else int(np.ceil(1.0 / (2 * eps**2) * np.log(2 / delta)))
    if n_samples < 1:
        raise ValueError(f"n_samples must be positive, but got {n_samples}")
    batch_size = int(min(batch_size, n_samples)) if batch_size and batch_size > 0 else min(100, n_samples)

    # generate samples
    n_nodes = 2 * f.size
    samples = np.zeros((n_samples,), dtype=float)
    buffer = pcms.tree.make_buffer(n_trees=batch_size, n_nodes=n_nodes)
    n_iter = n_samples // batch_size + 1
    seeds = [None] * n_iter if seed is None else generate_seeds(seed=seed, n=n_iter)
    for i in range(0, n_samples, batch_size):
        batch_size_ = min(batch_size, n_samples - i)
        samples[i : i + batch_size_] = pcms._haar.rand_dh_component(f, batch_size_, buffer, seeds[i // batch_size])

    # compute estimate
    samples.sort()
    result = np.searchsorted(samples, ys, side='left') / n_samples
    return result.item() if np.ndim(result) == 0 or result.size == 1 else result


def cdf_proj_dod(
    ys: ArrayLike,
    tree: pcms.tree.Tree,
    node: int,
    eps: float = 0.005,
    delta: float = 0.001,
    n_samples: Optional[int] = None,
    seed: Optional[int] = None
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
    eps: float (optional, default 0.005)
        Desired precision of the estimate.
    delta: float (optional, default 0.001)
        Probability that the estimate lies within the desired precision.
    n_samples: int (optional, default None)
        Number of samples to draw for the estimate (overrides eps and delta).
    seed: int (optional, default None)
        Random seed.

    Returns
    -------
    float or np.ndarray
        Value(s) of the CDF at the specified point(s).
    """
    ys = np.asarray(ys)
    n_leaves = tree.find_n_leaves()
    n_samples = np.ceil(1.0 / (2 * eps**2) * np.log(2 / delta)) if n_samples is None else n_samples
    subtree_size = tree.get_subtree_size(node)
    split_size = tree.get_subtree_size(tree.get_child(node))
    m = subtree_size - split_size
    a = np.sqrt(m / split_size / subtree_size)
    b = -np.sqrt(split_size / m / subtree_size)
    alpha = np.array([split_size, m, n_leaves - m - split_size])
    rng = np.random.default_rng(seed)
    f = rng.gamma(shape=alpha, size=(int(n_samples), len(alpha)))
    f /= f.sum(axis=1, keepdims=True)
    g = rng.gamma(shape=alpha, size=(int(n_samples), len(alpha)))
    g /= g.sum(axis=1, keepdims=True)
    coef = np.sqrt(m * split_size / subtree_size)
    samples = coef * (a * (f[:,0] - g[:,0]) + b * (f[:,1] - g[:,1]))
    samples.sort()
    result = np.searchsorted(samples, ys, side='left') / n_samples
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
    tree = pcms.tree.nwk2tree(input, ensure_planted=True) if isinstance(input, str) else input
    if not tree.find_is_planted():
        raise ValueError("Tree must be planted.")
    Q, S = pcms._haar.sparsify(tree._tree)
    n_leaves = tree.find_n_leaves()
    n_wavelets = tree.find_n_wavelets()
    Q = csc_matrix(Q, shape=(n_leaves, n_wavelets))
    S = csc_matrix(S, shape=(n_wavelets, n_wavelets))
    return Q, S