# tests/test_tree_numerics.py

import numpy as np


# --- Tests for branch length and total path length computations ---

def test_find_epl(small_structured_tree):
    epl = small_structured_tree.find_epl()
    assert epl == 4

def test_find_tbl_one_arg(small_structured_tree):
    small_structured_tree.set_edge_length(0, 1.5)
    tbl = small_structured_tree.find_tbl(0)
    assert tbl == 1.5

def test_find_tbl_two_args(small_structured_tree):
    small_structured_tree.set_edge_length(0, 1.5)
    small_structured_tree.set_edge_length(1, 2)
    tbl = small_structured_tree.find_tbl(0, 1)
    assert tbl == 3.5


# --- Tests for compute_wavelets method ---

def test_compute_wavelets_binary_daughter(small_structured_tree):
    wavelet = small_structured_tree.compute_wavelets(2)
    expected = np.array([np.sqrt(0.5), -np.sqrt(0.5)])
    assert np.array_equal(wavelet, expected)

def test_compute_wavelets_ternary_daughter(small_structured_ternary_tree):
    wavelets = small_structured_ternary_tree.compute_wavelets(3)
    expected_0 = np.array([np.sqrt(1/2), -np.sqrt(1/2)])
    expected_1 = np.array([np.sqrt(1/6), np.sqrt(1/6), -np.sqrt(2/3)])
    assert np.array_equal(wavelets[0], expected_0)
    assert np.array_equal(wavelets[1], expected_1)

def test_compute_wavelets_mother(small_structured_tree):
    wavelet = small_structured_tree.compute_wavelets(3)
    expected = np.array([1/np.sqrt(2), 1/np.sqrt(2)])
    assert np.array_equal(wavelet, expected)
