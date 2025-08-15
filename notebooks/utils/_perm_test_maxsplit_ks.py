# utils/_perm_test_maxsplit_ks_fast.py
import numpy as np
from typing import Optional, Tuple, Sequence
from numpy.typing import ArrayLike, NDArray


def _ks_for_m_using_sorted_order(order: NDArray[np.intp], m: int) -> float:
    """
    Compute KS statistic between seq[:m] and seq[m:] using precomputed 'order'.
    'order' is an array of original indices sorted by value (len n).
    This function uses O(n) time and O(n) memory (no n^2 allocations).
    """
    n = len(order)
    if m == 0 or m == n:
        return 0.0
    # mask_sorted[k] = True iff sorted element k belongs to left group (original idx < m)
    mask_sorted = order < m                 # boolean length n
    left_counts = np.cumsum(mask_sorted)    # length n, int
    # left_ecdf at sorted positions: left_counts / m
    # right_ecdf at sorted positions: (k+1 - left_counts) / (n-m)
    kplus1 = np.arange(1, n + 1)
    left_ecdf = left_counts / float(m)
    right_ecdf = (kplus1 - left_counts) / float(n - m)
    return float(np.max(np.abs(left_ecdf - right_ecdf)))


def max_split_ks(
    seq: ArrayLike,
    min_size: int = 5,
    grid_size: int = 1000,
    refine_radius: int = 200,
    top_k: int = 3
) -> Tuple[float, int]:
    """
    Coarse-to-fine search for the split maximizing the KS statistic.

    Parameters
    ----------
    seq : array_like
        Ordered 1D data (sequence).
    min_size : int
        Minimum points allowed in each side of split.
    grid_size : int
        Number of equally spaced candidate split points to evaluate in coarse pass.
    refine_radius : int
        Half-width of window around each top candidate to scan precisely.
    top_k : int
        Number of top coarse candidates to refine (choose small, e.g., 1..5).

    Returns
    -------
    best_val : float
        Maximum KS statistic (approx/then refined to exact within neighborhoods).
    best_m : int
        Split point index achieving the best_val (0-based, so left = seq[:best_m]).
    """
    seq = np.asarray(seq)
    n = len(seq)
    if n < 2 * min_size:
        return 0.0, 0
    order = np.argsort(seq)

    # If n small, do exact exhaustive scan (fast enough)
    if n <= 5000:
        best = 0.0
        best_m = min_size
        for m in range(min_size, n - min_size + 1):
            v = _ks_for_m_using_sorted_order(order, m)
            if v > best:
                best = v
                best_m = m
        return best, best_m

    # --- Coarse grid (equally spaced split points) ---
    # ensure grid points respect min_size and integer uniqueness
    low = min_size
    high = n - min_size
    if grid_size >= (high - low + 1):
        coarse_ms = np.arange(low, high + 1)
    else:
        coarse_ms = np.linspace(low, high, num=grid_size, dtype=int)
        coarse_ms = np.unique(coarse_ms)  # remove duplicates

    coarse_vals = np.empty(len(coarse_ms), dtype=float)
    for idx, m in enumerate(coarse_ms):
        coarse_vals[idx] = _ks_for_m_using_sorted_order(order, int(m))

    # pick top candidates and refine
    top_idx = np.argsort(coarse_vals)[-top_k:][::-1]
    best = 0.0
    best_m = min_size
    for idx in top_idx:
        m0 = int(coarse_ms[idx])
        lo = max(min_size, m0 - refine_radius)
        hi = min(n - min_size, m0 + refine_radius)
        # exact local scan
        for m in range(lo, hi + 1):
            v = _ks_for_m_using_sorted_order(order, m)
            if v > best:
                best = v
                best_m = m

    return float(best), int(best_m)


def perm_test_maxsplit_ks(
    seq: ArrayLike,
    min_size: int = 5,
    eps: float = 0.005,
    delta: float = 0.01,
    grid_size: int = 1000,
    refine_radius: int = 200,
    top_k: int = 3,
    rng: Optional[np.random.Generator] = None
) -> Tuple[float, float, NDArray[np.float64]]:
    """
    Permutation test using the coarse → refine KS scan, with adaptive number of
    permutations determined by eps (accuracy) and delta (failure probability).

    Parameters
    ----------
    seq : array_like
        Ordered sequence to test.
    min_size : int
        Minimum number of points in each segment.
    eps : float
        Desired accuracy for p-value estimate.
    delta : float
        Desired confidence (1-delta) that p-value is within eps of true p.
    grid_size : int
        Number of coarse split points to evaluate.
    refine_radius : int
        Neighborhood radius for local refinement.
    top_k : int
        Number of top coarse splits to refine.
    rng : optional numpy Generator

    Returns
    -------
    T_obs : float
        Observed max-split KS statistic.
    p_value : float
        Permutation p-value estimate.
    perm_stats : ndarray
        KS statistics from permuted sequences (null distribution).
    """
    rng = np.random.default_rng(rng)
    seq = np.asarray(seq)
    n = len(seq)

    # Number of permutations derived from Hoeffding bound
    B = int(np.ceil(np.log(2 / delta) / (2 * eps ** 2)))

    # Observed statistic
    T_obs, _ = max_split_ks(
        seq,
        min_size=min_size,
        grid_size=grid_size,
        refine_radius=refine_radius,
        top_k=top_k
    )

    Ts = np.empty(B, dtype=float)
    for b in range(B):
        perm = rng.permutation(n)
        Ts[b], _ = max_split_ks(
            seq[perm],
            min_size=min_size,
            grid_size=grid_size,
            refine_radius=refine_radius,
            top_k=top_k
        )

    p_value = (1.0 + np.sum(Ts >= T_obs)) / (B + 1.0)
    return float(T_obs), float(p_value), Ts

