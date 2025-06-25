# tests/test_tree_relationships.py

import numpy as np


# --- Tests for basic node manipulation ---

def test_link_nodes(small_unstructured_tree):
    small_unstructured_tree.link(0, 2)
    small_unstructured_tree.link(1, 2)
    assert small_unstructured_tree.get_parent(0) == 2
    assert small_unstructured_tree.get_sibling(1) == 0

def test_cut_nodes(small_structured_tree):
    small_structured_tree.cut(0)
    assert small_structured_tree.get_parent(0) == -1
    assert small_structured_tree.get_sibling(0) == -1

def test_swap_nodes(small_structured_tree):
    small_structured_tree.swap(0, 1)
    assert small_structured_tree.get_sibling(1) == 0


# --- Tests for parent, child, sibling queries ---

def test_get_parent_on_small_structured_tree(small_structured_tree):
    assert small_structured_tree.get_parent(0) == 2
    assert small_structured_tree.get_parent(1) == 2
    assert small_structured_tree.get_parent(3) == -1

def test_get_child_on_small_structured_tree(small_structured_tree):
    child = small_structured_tree.get_child(2)
    assert child == 0

def test_get_sibling_on_small_structured_tree(small_structured_tree):
    sibling_0 = small_structured_tree.get_sibling(0)
    sibling_1 = small_structured_tree.get_sibling(1)
    assert sibling_0 == 1
    assert sibling_1 == -1

def test_get_is_first(small_structured_tree):
    val0 = small_structured_tree.get_is_first(0)
    val1 = small_structured_tree.get_is_first(1)
    assert val0 is True
    assert val1 is False


# --- Tests for tree traversal and structure queries ---

def test_find_children(small_structured_tree):
    children = small_structured_tree.find_children(2)
    assert set(children) == {0, 1}
    assert small_structured_tree.find_children(0).size == 0

def test_find_ancestors(small_structured_tree):
    assert np.array_equal(small_structured_tree.find_ancestors(0), np.array([2, 3]))
    assert np.array_equal(small_structured_tree.find_ancestors(2), np.array([3], dtype=int))
    assert np.array_equal(small_structured_tree.find_ancestors(3), np.array([], dtype=int))

def test_find_path(large_structured_tree):
    there, back = large_structured_tree.find_path(0, 3)
    assert len(there) == 1
    assert np.array_equal(there, np.array([2]))
    assert np.array_equal(back, np.array([], dtype=int))

def test_find_root(small_structured_tree):
    root = small_structured_tree.find_root()
    assert root == 3

def test_find_leaves(small_structured_tree):
    leaves, depths = small_structured_tree.find_leaves(None)
    assert set(leaves) == {0, 1}
    assert all(isinstance(d, (int, np.integer)) for d in depths)

def test_find_subtree_start_indices(small_structured_tree):
    start_indices = small_structured_tree.find_subtree_start_indices()
    assert len(start_indices) == 4
    assert all(0 <= i < 2 for i in start_indices)
