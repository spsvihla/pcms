"""
@file haar.py
@brief A .py wrapper for the pcms._haar package.
@author Sean Svihla
"""

from typing import Union, Sequence

import numpy as np
from scipy.special import factorial
from scipy.stats import beta

import pcms._haar 
from pcms.tree import CriticalBetaSplittingDistribution 


def cdf_rand_basis(
    ys: Union[float, Sequence[float], np.ndarray],
    func: Union[Sequence[float], np.ndarray],
    eps: float,
    delta: float,
    seed: int = 0
) -> Union[float, np.ndarray]:
    """
    @brief Evaluates the CDF of <f,φ> where φ is the Haar-like wavelet associated
           with an interior node 'v' of a critical beta-splitting tree T and f 
           is a non-random function on the leaves of the sub-tree T(v). The 
           split below φ varies randomly according to the critical beta-splitting
           split distribution q(n,i).

    @param y Point(s) to evaluate the CDF.
    @param func Function on the leaves of T(v).
    @param eps Desired precision.
    @param delta Probability estimate falls within desired precision.
    @param seed A random seed.

    @return Estimate(s) of the CDF at the point(s).
    """
    ys_arr = np.atleast_1d(ys).astype(np.float64)
    func_arr = np.asarray(func, dtype=np.float64)
    n = func_arr.size
    num_iter = int(np.ceil(1.0 / (4 * eps**2 * delta)))
    if factorial(n) < num_iter:
        num_iter = n
    pmf = CriticalBetaSplittingDistribution(n).get_pmf()
    result = pcms._haar.cdf_rand_basis(ys_arr, func_arr, pmf, int(num_iter), int(seed))
    return result.item() if np.isscalar(ys) else result


def cdf_rand_func(
    ys: Union[float, Sequence[float], np.ndarray],
    split_size: int,
    subtree_size: int
) -> float:
    """
    @brief Evaluates the CDF of <f,φ> where φ us the Haar-like wavelet associated
           with an interior node with 'subtree_size' descendent leaves and a 
           split of size 'split_size'. In this case, the distribution is 
           an affine transform of a Beta(split_size, subtree_size-split_size)
           distribution

    @param ys Point(s) to evaluate the CDF.
    @param split_size The split size beneath the node.
    @param subtree_size The size of the subtree rooted at the node.

    @return Estimate(s) of the CDF at the point(s).
    """
    a = np.sqrt((subtree_size - split_size) / split_size / subtree_size)
    b = np.sqrt(split_size / (subtree_size - split_size) / subtree_size)
    return beta.cdf((ys + b) / (a + b), split_size, subtree_size - split_size)