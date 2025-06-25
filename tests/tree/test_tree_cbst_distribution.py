# tests/test_tree_cbst_distribution.py

import pytest
import numpy as np
from pcms.tree import CriticalBetaSplittingDistribution

# ---------- Fixtures ----------

@pytest.fixture
def small_dist():
    return CriticalBetaSplittingDistribution(5)

@pytest.fixture
def larger_dist():
    return CriticalBetaSplittingDistribution(10)

# ---------- Construction & Validation ----------

@pytest.mark.parametrize("invalid_n", [0, 1, -5])
def test_invalid_n_raises(invalid_n):
    with pytest.raises(ValueError):
        CriticalBetaSplittingDistribution(invalid_n)

@pytest.mark.parametrize("invalid_type", ["ten", 3.5, None])
def test_invalid_type_raises(invalid_type):
    with pytest.raises(TypeError):
        CriticalBetaSplittingDistribution(invalid_type)

# ---------- Sampling ----------

def test_sample_type_and_range(small_dist):
    sample = small_dist.sample()
    assert isinstance(sample, int)
    assert 1 <= sample <= 4

def test_multiple_samples_have_variety(larger_dist):
    samples = [larger_dist.sample() for _ in range(1000)]
    assert all(isinstance(x, int) for x in samples)
    assert all(1 <= x <= 9 for x in samples)
    assert len(set(samples)) > 1  # samples not all identical

# ---------- Probability Mass Function (pmf) ----------

def test_pmf_length_and_nonnegative(small_dist):
    pmf = small_dist.pmf
    assert isinstance(pmf, (list, np.ndarray))
    assert len(pmf) == 4
    assert all(p >= 0 for p in pmf)
    np.testing.assert_almost_equal(sum(pmf), 1.0)

def test_pmf_exact_values():
    # Suppose n=5, known exact pmf from formula or precomputed
    expected_pmf = np.array([
        0.3,   # example values; replace with real ones you provide
        0.2,
        0.2,
        0.3
    ])
    dist = CriticalBetaSplittingDistribution(5)
    np.testing.assert_allclose(dist.pmf, expected_pmf, rtol=1e-7)

# ---------- Cumulative Distribution Function (cdf) ----------

def test_cdf_length_and_monotonic(larger_dist):
    cdf = larger_dist.cdf
    assert isinstance(cdf, (list, np.ndarray))
    assert len(cdf) == 9
    assert all(x <= y for x, y in zip(cdf, cdf[1:])), "CDF must be non-decreasing"
    assert 0 <= cdf[0] <= 1
    assert 0 <= cdf[-1] <= 1

def test_cdf_exact_values():
    # Suppose n=5, known exact cdf from formula or precomputed
    expected_cdf = np.array([
        0.3,    # example cumulative sums matching pmf above
        0.5,
        0.7,
        1.0
    ])
    dist = CriticalBetaSplittingDistribution(5)
    np.testing.assert_allclose(dist.cdf, expected_cdf, rtol=1e-7)

# ---------- Consistency Between PMF and CDF ----------

def test_cdf_consistent_with_pmf(small_dist):
    pmf = np.array(small_dist.pmf)
    cdf = np.array(small_dist.cdf)
    np.testing.assert_allclose(np.cumsum(pmf), cdf, rtol=1e-10)

# ---------- Edge Cases ----------

def test_sampling_edge_cases():
    # For a moderately large n, sample should still work
    dist = CriticalBetaSplittingDistribution(1000)
    sample = dist.sample()
    assert isinstance(sample, int)
    assert 1 <= sample <= 999

def test_pmf_and_cdf_types_and_lengths():
    n = 20
    dist = CriticalBetaSplittingDistribution(n)
    pmf = dist.pmf
    cdf = dist.cdf
    assert len(pmf) == n - 1
    assert len(cdf) == n - 1
    assert isinstance(pmf, (list, np.ndarray))
    assert isinstance(cdf, (list, np.ndarray))
