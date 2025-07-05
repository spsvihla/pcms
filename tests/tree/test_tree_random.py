# tests/test_tree_random.py

import pytest
from pcms.tree import remy, cbst, Tree


# --- Basic functionality ---

def test_remy_valid_output_type_and_leaves():
    n = 5
    tree = remy(n)
    assert isinstance(tree, Tree)
    leaves, _ = tree.find_leaves(None)
    assert len(leaves) == n

def test_cbst_valid_output_type_and_leaves():
    n = 7
    tree = cbst(n)
    assert isinstance(tree, Tree)
    leaves, _ = tree.find_leaves(None)
    assert len(leaves) == n

def test_remy_without_seed_runs():
    tree = remy(4)
    assert isinstance(tree, Tree)

def test_cbst_without_seed_runs():
    tree = cbst(4)
    assert isinstance(tree, Tree)


# --- Parameterized size tests ---

@pytest.mark.parametrize("func", [remy, cbst])
@pytest.mark.parametrize("n", [2, 3, 10, 20])
def test_tree_generation_various_sizes(func, n):
    tree = func(n)
    leaves, _ = tree.find_leaves(None)
    assert len(leaves) == n
    assert isinstance(tree, Tree)


# --- Seed reproducibility tests ---

@pytest.mark.parametrize("func", [remy, cbst])
def test_seed_reproducibility(func):
    n = 10
    seed = 12345
    tree1 = func(n, seed=seed)
    tree2 = func(n, seed=seed)
    # Trees generated with the same seed should be identical
    assert str(tree1) == str(tree2)

@pytest.mark.parametrize("func", [remy, cbst])
def test_different_seeds_produce_different_trees(func):
    n = 10
    tree1 = func(n, seed=123)
    tree2 = func(n, seed=456)
    # Trees generated with different seeds should generally differ
    assert str(tree1) != str(tree2)


# --- Invalid inputs ---

@pytest.mark.parametrize("func", [remy, cbst])
@pytest.mark.parametrize("bad_n", [0, 1, -5])
def test_invalid_n_leaves_raises_value_error(func, bad_n):
    with pytest.raises(ValueError):
        func(bad_n)
