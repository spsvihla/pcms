"""
@file stats.py
@brief Functions for performing statistical tests on Haar-like coordinates and 
       Haar-like distance components.
@author Sean Svihla
"""
from typing import Optional

import numpy as np
from numpy.typing import ArrayLike
import mvhg

import pcms.tree
import pcms.haar


def cdf_haar_coordinate_perm(
    ys: float | ArrayLike,
    tree: pcms.tree.Tree,
    nodes: int | ArrayLike,
    counts: ArrayLike,
    num_samples: Optional[int] = None,
    num_max_iter: Optional[int] = None,
    batch_size: Optional[int] = None,
    eps: Optional[float] = 0.005,
    delta: Optional[float] = 0.001,
    seed: Optional[int] = None
) -> float | ArrayLike:
    """
    Estimate the empirical cumulative distribution function of Haar-like coordinates
    via a permutation test using the multivariate hypergeometric distribution.

    Parameters
    ----------
    ys: float or ArrayLike
        Values at which to evaluate the ECDF.
    tree: pcms.tree.Tree
        Phylogenetic tree object.
    nodes: int or ArrayLike
        Internal node indices for which Haar coordinates are computed.
    counts: ArrayLike
        Leaf counts, shape (K,2), where K = number of leaves.
    num_samples: int, default = None
        Number of Monte Carlo samples to draw (overrides eps and delta). If
        None, num_samples is calculated from eps and delta according to the
        DKW inequality.
    num_iter_max: int, default = None
        Maximum number of iterations of the rejection sampler. If None, this 
        value is the default given by mvhg.multivariate_hypergeometric.
    batch_size: int, default = None
        Number of samples to collect at once. If None, no batching is used and
        all samples are computed at once.
    eps: float, default = 0.005
        Desired error tolerance for number of permutations.
    delta: float, default = 0.001
        Confidence parameter for number of permutations.
    seed: int, optional
        Random seed for reproducibility, by default None.

    Returns
    -------
    float or ArrayLike
        Empirical CDF evaluated at ys for each node.
    """
    if np.isscalar(ys):
        ys = np.array([ys], dtype=float)
    elif ys.ndim in (1, 2):
        ys = ys.astype(float, copy=False)
    else:
        raise ValueError(f"'ys' must have shape (N,) or (N, M), got shape {ys.shape}")

    if np.isscalar(nodes):
        nodes = np.array([nodes], dtype=int)
    elif nodes.ndim == 1:
        nodes = nodes.astype(int, copy=False)
    else:
        raise ValueError(f"'nodes' must have shape (N,), got shape {nodes.shape}")

    if ys.shape[0] != nodes.shape[0]:
        raise ValueError(f"Inconsistent lengths: ys.shape[0]={ys.shape[0]}, nodes.shape[0]={nodes.shape[0]}")

    intr_nodes = tree.find_interior_nodes()
    if not np.isin(nodes, intr_nodes).all():
        invalid = nodes[~np.isin(nodes, intr_nodes)]
        raise ValueError(f"Invalid node IDs: {invalid.tolist()} not found in intr_nodes")

    leaves = tree.find_leaves()
    K = leaves.size
    if counts.shape != (K, 2):
        raise ValueError(f"'counts' must have shape (K,2), got {counts.shape}")

    if not np.issubdtype(counts.dtype, np.integer):
        if np.allclose(counts, np.round(counts)):
            counts = counts.astype(int)
        else:
            raise ValueError(f"'counts' must contain integers, got dtype {counts.dtype}")

    if (counts < 0).any():
        bad = np.where(counts < 0)
        raise ValueError(f"'counts' contains negative entries at positions {bad}")

    # sampling
    if num_samples is None:
        num_samples = int(np.ceil((1 / (2 * eps**2)) * np.log(2 / delta)))
    else:
        num_samples = int(num_samples)

    Ns = np.sum(counts, axis=1)
    leaf_mask = Ns != 0
    Ns = Ns[leaf_mask]
    N = np.sum(Ns)
    Na = np.sum(counts[:,0])

    # basis matrix
    node_mask = np.isin(intr_nodes, nodes)
    Q = pcms.haar.sparsify(tree)[0][leaf_mask, :][:, node_mask]

    def collect_totals(num_samples):
        A_counts = mvhg.multivariate_hypergeometric(Ns, N, Na,
                                                    num_samples=num_samples,
                                                    num_max_iter=num_max_iter,
                                                    seed=seed)
        B_counts = counts[leaf_mask, 1] - A_counts
        A_abundances = A_counts / Na
        B_abundances = B_counts / (N - Na)
        abundance_differences = A_abundances - B_abundances

        # compute Haar-like coordinates
        coordinates = Q.T @ abundance_differences.T
        coordinates = np.abs(coordinates)

        # estimate c.d.f.
        if ys.ndim == 1:
            totals = np.count_nonzero(coordinates <= ys[:, None], axis=1)
        else:
            totals = np.count_nonzero(coordinates[:, :, None] <= ys[:, None, :], axis=1)

        return totals
    
    if batch_size is None:
        totals = collect_totals(num_samples)
        return totals / num_samples
    else:
        num_batches = num_samples // batch_size
        num_samples_remaining = num_samples - num_batches * batch_size
        totals = collect_totals(num_samples_remaining)
        for _ in range(num_batches):
            totals += collect_totals(batch_size)
        return totals / num_samples
