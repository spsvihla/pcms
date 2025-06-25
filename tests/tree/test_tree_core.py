# tests/test_tree_core.py

import pytest
import pcms.tree
import numpy as np


# --- Tree Initialization Tests ---

def test_tree_init_with_invalid_str():
    with pytest.raises(TypeError):
        pcms.tree.Tree("3")

def test_tree_init_with_none():
    with pytest.raises(TypeError):
        pcms.tree.Tree(None)

def test_tree_init_with_valid_int(singleton_tree):
    assert isinstance(singleton_tree, pcms.tree.Tree)

def test_del_method():
    t = pcms.tree.Tree(1)
    del t


# --- Tree Creation from C++ Object ---

def test_from_cpp_tree_creates_valid_instance(small_structured_tree):
    t2 = pcms.tree.Tree._from_cpp_tree(small_structured_tree._tree)
    assert isinstance(t2, pcms.tree.Tree)

def test_from_cpp_tree_initializes_attributes(small_structured_tree):
    t2 = pcms.tree.Tree._from_cpp_tree(small_structured_tree._tree)
    assert hasattr(t2, "_tree")
    assert t2.print_label == "none"


# --- Basic Tree Properties ---

def test_tree_node_count():
    t = pcms.tree.Tree(5)
    assert t.n_nodes == 5

def test_print_label_default(singleton_tree):
    assert singleton_tree.print_label == "none"

@pytest.mark.parametrize("label", ["index", "name", "none"])
def test_print_label_setter_valid_input(singleton_tree, label):
    singleton_tree.print_label = label
    assert singleton_tree.print_label == label

def test_print_label_setter_invalid_input(singleton_tree):
    with pytest.raises(ValueError):
        singleton_tree.print_label = "invalid_label"


# --- String Representation ---

def test_tree_str_representation_structured(small_structured_tree):
    s = str(small_structured_tree)
    assert isinstance(s, str)

def test_tree_str_representation_unstructured(small_unstructured_tree):
    s = str(small_unstructured_tree)
    assert isinstance(s, str)


# --- Edge Length Tests ---

def test_set_edge_length(small_structured_tree):
    for u in range(small_structured_tree.n_nodes):
        small_structured_tree.set_edge_length(u, 3.14)
        val = small_structured_tree.get_edge_length(u)
        assert val == 3.14

def test_get_edge_length_returns_float(small_structured_tree):
    for u in range(small_structured_tree.n_nodes):
        val = small_structured_tree.get_edge_length(u)
        assert isinstance(val, float)

def test_get_edge_length_vector_returns_vector(small_structured_tree):
    lengths = small_structured_tree.get_edge_length(None)
    assert hasattr(lengths, "__len__")

def test_set_edge_length_invalid_index(small_unstructured_tree):
    with pytest.raises(IndexError):
        small_unstructured_tree.set_edge_length(-1, 1.0)
    with pytest.raises(IndexError):
        small_unstructured_tree.set_edge_length(5, 1.0)


# --- Subtree Size Tests ---

def test_get_subtree_size_returns_int(small_structured_tree):
    for u in range(small_structured_tree.n_nodes):
        size = small_structured_tree.get_subtree_size(u)
        assert isinstance(size, int)

def test_get_subtree_size_vector_returns_vector(small_structured_tree):
    full_sizes = small_structured_tree.get_subtree_size(None)
    assert hasattr(full_sizes, "__len__")


# --- Boolean and Array Return Tests ---

def test_get_is_first_returns_bool(small_structured_tree):
    val = small_structured_tree.get_is_first(0)
    assert isinstance(val, bool)

def test_find_children_returns_array(small_structured_tree):
    for u in range(small_structured_tree.n_nodes):
        children = small_structured_tree.find_children(u)
        assert hasattr(children, "__len__")

def test_find_ancestors_returns_array(small_structured_tree):
    for u in range(small_structured_tree.n_nodes):
        ancestors = small_structured_tree.find_ancestors(u)
        assert hasattr(ancestors, "__len__")


# --- Path, Root, Leaves, and Indices ---

def test_find_path_returns_tuple(small_structured_tree):
    if small_structured_tree.n_nodes > 1:
        path = small_structured_tree.find_path(0, 1)
        assert isinstance(path, tuple)
        assert len(path) == 2

def test_find_root_returns_int(small_structured_tree):
    root = small_structured_tree.find_root()
    assert isinstance(root, int)

def test_find_leaves_returns_tuple(small_structured_tree):
    leaves, depths = small_structured_tree.find_leaves(None)
    assert hasattr(leaves, "__len__")
    assert hasattr(depths, "__len__")

def test_find_subtree_start_indices_returns_array(small_structured_tree):
    indices = small_structured_tree.find_subtree_start_indices()
    assert hasattr(indices, "__len__")

def test_find_leaves_returns_ints(small_structured_tree):
    leaves, depths = small_structured_tree.find_leaves(None)
    assert all(isinstance(l, (int, np.integer)) for l in leaves)
    assert all(isinstance(d, (int, np.integer)) for d in depths)


# --- Numeric Return Type Checks ---

def test_find_epl_returns_numeric(small_structured_tree):
    epl = small_structured_tree.find_epl()
    assert isinstance(epl, (int, float))

def test_find_tbl_no_args_returns_vector(small_structured_tree):
    val = small_structured_tree.find_tbl(None, None)
    assert hasattr(val, "__len__")

def test_find_tbl_one_arg_returns_float(small_structured_tree):
    val = small_structured_tree.find_tbl(0, None)
    assert isinstance(val, float)

def test_find_tbl_two_args_returns_float(small_structured_tree):
    if small_structured_tree.n_nodes > 1:
        val = small_structured_tree.find_tbl(0, 1)
        assert isinstance(val, float)


# --- Compute Wavelets Return Types ---

def test_compute_wavelets_binary_daughter_returns_vector(small_structured_tree):
    wavelet = small_structured_tree.compute_wavelets(2)
    assert isinstance(wavelet, np.ndarray)

def test_compute_wavelets_ternary_daughter_returns_tuple(small_structured_ternary_tree):
    wavelets = small_structured_ternary_tree.compute_wavelets(3)
    assert isinstance(wavelets, list)
    assert all(isinstance(w, np.ndarray) for w in wavelets)

def test_compute_wavelets_mother_returns_vector(small_structured_tree):
    wavelet = small_structured_tree.compute_wavelets(3)
    assert isinstance(wavelet, np.ndarray)


# --- Invalid Index Error Tests ---

@pytest.mark.parametrize("method,args", [
    ("get_parent", [-1]),
    ("get_parent", [10]),
    ("get_child", [-1]),
    ("get_child", [10]),
    ("get_sibling", [-1]),
    ("get_sibling", [5]),
    ("set_edge_length", [-1, 1.0]),
    ("set_edge_length", [5, 1.0]),
    ("get_name", [-1]),
    ("get_name", [5]),
    ("set_name", [-1, "bad"]),
    ("set_name", [5, "bad"]),
    ("link", [0, 10]),
    ("link", [-1, 0]),
    ("cut", [-1]),
    ("cut", [5]),
    ("swap", [-1, 1]),
    ("swap", [0, 5]),
    ("find_children", [-1]),
    ("find_children", [100]),
    ("find_ancestors", [-1]),
    ("find_ancestors", [100]),
    ("find_path", [-1, 1]),
    ("find_path", [0, 100]),
    ("compute_wavelets", [-1]),
    ("compute_wavelets", [10]),
])
def test_methods_raise_index_error(small_unstructured_tree, method, args):
    func = getattr(small_unstructured_tree, method)
    with pytest.raises(IndexError):
        func(*args)


# --- RuntimeError Tests ---

def test_link_runtime_error(small_structured_tree):
    with pytest.raises(RuntimeError):
        small_structured_tree.link(0, 1)


# --- Compute Wavelets ValueError Test ---

def test_compute_wavelets_value_error_leaf_node(small_structured_tree):
    with pytest.raises(ValueError):
        small_structured_tree.compute_wavelets(0)
