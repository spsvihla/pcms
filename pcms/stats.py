"""
@file stats.py
@brief Functions for performing statistical tests on Haar-like coordinates and 
       Haar-like distance components.
@author Sean Svihla
"""

from typing import Optional

import numpy as np
from numpy.typing import ArrayLike

import pcms.tree


def cdf_coord_dod(
    ys: ArrayLike,
    tree: pcms.tree.Tree,
    node: int,
    eps: float = 0.005,
    delta: float = 0.001,
    n_samples: Optional[int] = None,
    seed: Optional[int] = None
) -> float:
    """
    Evaluate the c.d.f. of ⟨f, φ⟩ where φ is a Haar-like wavelet associated
    with 'node' and 'f' is a difference of Dirichlet(1,1,...,1) random 
    variables.

    Parameters
    ----------
    ys: float or array-like
        Points at which to evaluate the c.d.f.
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
    ys = np.asarray(ys).ravel().astype(np.float64)
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
    return result.item() if ys.size == 1 else result