import os
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.integrate import quad

from typing import Optional, Tuple
from numpy.typing import NDArray

import pcms.tree


##
# Helper functions
##
def _integrand(x, n):
    return (1 - x**n) / (1 - x)


def _harmonic_number(n):
    hn, _ = quad(_integrand, 0, 1, args=(n,))
    return hn


x_cache = {(26, 25, 5): -0.016796952511952505}

def compute_xn(n, k=25, m=5):
    for current_n in range(k + 1, n + 1):
        if (current_n, k, m) in x_cache:
            continue

        s = sum(x_cache[(current_n - i, k, m)] / i for i in range(1, current_n - k))
        s /= _harmonic_number(current_n - 1)
        x_cache[(current_n, k, m)] = s

    return x_cache[(n, k, m)]


def draw_subtree_size_counts(
    n_samples: int, n_leaves: int
) -> Tuple[NDArray, NDArray]:
    trees = pcms.tree.cbst(n_leaves=n_leaves, planted=False, n_samples=n_samples)
    counts = Counter()
    for t in trees:
        s, c = np.unique_counts(t.get_subtree_size())
        counts.update(dict(zip(s, c)))
    return np.array(list(counts.keys())), np.array(list(counts.values()))


def _worker_task(
    n_samples: int, n_leaves: int, k: int, m: int, batch_size: int, seed: Optional[int] = None
) -> Counter:
    if seed is not None:
        np.random.seed(seed)

    counts = Counter()
    for start in range(0, n_samples, batch_size):
        chunk_size = min(batch_size, n_samples - start)
        trees = pcms.tree.cbst(n_leaves=n_leaves, planted=False, n_samples=chunk_size)

        for t in trees:
            sizes = t.get_subtree_size()
            n_k = int(np.count_nonzero(sizes == k))
            n_m = n_k if k == m else int(np.count_nonzero(sizes == m))
            counts[(n_k, n_m)] += 1

    return counts


def _get_available_memory_bytes() -> int:
    try:
        import psutil
        return psutil.virtual_memory().available
    except ImportError:
        try:
            return os.sysconf("SC_AVPHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
        except (ValueError, AttributeError):
            # Conservative 4 GB baseline fallback
            return 4 * 1024 * 1024 * 1024


def _compute_safe_batch_size(
    n_leaves: int, n_jobs: int, requested_batch_size: int = -1, safety_margin: float = 0.70
) -> int:
    available_ram = _get_available_memory_bytes()
    
    # Cap total target usage to a fraction of free RAM (e.g., 70% safety headroom)
    ram_per_worker = (available_ram * safety_margin) / max(1, n_jobs)

    # A binary tree with L leaves has ~2L total nodes. 
    # Estimating ~500 bytes per leaf node (includes internal node overhead + C++/Python wrappers).
    estimated_bytes_per_tree = max(1024, n_leaves * 500)

    max_safe_batch = max(1, int(ram_per_worker / estimated_bytes_per_tree))

    if requested_batch_size <= 0:
        return max_safe_batch
    
    # Clamp explicit user-provided batch size if it exceeds safety limit
    return min(requested_batch_size, max_safe_batch)


def _dispatch_worker_batch(
    executor: ProcessPoolExecutor,
    total_to_draw: int,
    n_leaves: int,
    k: int,
    m: int,
    effective_batch_size: int,
    n_jobs: int,
    rng: np.random.Generator,
    total_counts: Counter,
) -> None:
    base_samples = total_to_draw // n_jobs
    remainder = total_to_draw % n_jobs
    worker_samples = [
        base_samples + (1 if i < remainder else 0) for i in range(n_jobs)
    ]
    seeds = rng.integers(0, 2**31 - 1, size=n_jobs)

    futures = [
        executor.submit(
            _worker_task, samples, n_leaves, k, m, effective_batch_size, seed
        )
        for samples, seed in zip(worker_samples, seeds)
        if samples > 0
    ]
    for f in futures:
        total_counts.update(f.result())


def draw_subtree_pair_size_counts(
    n_leaves: int,
    k: int,
    m: int,
    n_samples: Optional[int] = None,
    se_tol: Optional[float] = None,
    n_jobs: int = -1,
    batch_size: int = -1,
) -> Tuple[NDArray, NDArray]:
    if (n_samples is None and se_tol is None) or (
        n_samples is not None and se_tol is not None
    ):
        raise ValueError("Must specify exactly one of 'n_samples' or 'se_tol'.")

    if n_jobs == -1:
        n_jobs = os.cpu_count() or 1

    effective_batch_size = _compute_safe_batch_size(
        n_leaves=n_leaves, n_jobs=n_jobs, requested_batch_size=batch_size
    )

    total_counts = Counter()
    rng = np.random.default_rng()

    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        if n_samples is not None:
            _dispatch_worker_batch(
                executor,
                n_samples,
                n_leaves,
                k,
                m,
                effective_batch_size,
                n_jobs,
                rng,
                total_counts,
            )

    if not total_counts:
        return np.empty((0, 2), dtype=np.int64), np.empty((0,), dtype=np.int64)

    sorted_items = sorted(total_counts.items())
    unique_pairs = np.array([item[0] for item in sorted_items], dtype=np.int64)
    pair_counts = np.array([item[1] for item in sorted_items], dtype=np.int64)

    return unique_pairs, pair_counts