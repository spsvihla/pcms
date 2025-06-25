# tests/test_tree_nwk2tree.py

import pytest
import pcms.tree


# --- Helper to write and parse a Newick string ---

def write_and_parse(tmp_path, newick):
    file = tmp_path / "tree.nwk"
    file.write_text(newick)
    return pcms.tree.nwk2tree(file)


# --- Valid Newick string tests ---

def test_nwk2tree_single_leaf(tmp_path):
    tree = write_and_parse(tmp_path, "A;")
    assert tree.n_nodes == 1
    leaves, _ = tree.find_leaves(None)
    assert len(leaves) == 1

def test_nwk2tree_basic_binary(tmp_path):
    tree = write_and_parse(tmp_path, "(A,B);")
    assert tree.n_nodes == 3
    leaves, _ = tree.find_leaves(None)
    assert set(tree.get_name(i) for i in leaves) == {"A", "B"}

def test_nwk2tree_larger_binary(tmp_path):
    tree = write_and_parse(tmp_path, "((A,B),(C,D));")
    assert tree.n_nodes == 7
    leaves, _ = tree.find_leaves(None)
    assert set(tree.get_name(i) for i in leaves) == {"A", "B", "C", "D"}

def test_nwk2tree_named_internal_nodes(tmp_path):
    tree = write_and_parse(tmp_path, "((A,B)X,(C,D)Y)Root;")
    names = [tree.get_name(i) for i in range(tree.n_nodes)]
    assert "X" in names and "Y" in names and "Root" in names

def test_nwk2tree_with_branch_lengths(tmp_path):
    tree = write_and_parse(tmp_path, "(A:0.1,B:0.2):0.3;")
    assert tree.n_nodes == 3
    lengths = [tree.get_edge_length(i) for i in range(tree.n_nodes)]
    assert all(isinstance(l, (float, int)) for l in lengths)

def test_nwk2tree_nested_structure(tmp_path):
    tree = write_and_parse(tmp_path, "(((A,B),C),D);")
    assert tree.n_nodes == 7
    leaves, _ = tree.find_leaves(None)
    assert set(tree.get_name(i) for i in leaves) == {"A", "B", "C", "D"}

def test_nwk2tree_unlabeled_internal_nodes(tmp_path):
    tree = write_and_parse(tmp_path, "((A,B),C);")
    assert tree.n_nodes == 5
    leaves, _ = tree.find_leaves(None)
    assert set(tree.get_name(i) for i in leaves) == {"A", "B", "C"}

def test_nwk2tree_polytomy(tmp_path):
    tree = write_and_parse(tmp_path, "(A,B,C,D);")
    leaves, _ = tree.find_leaves(None)
    assert set(tree.get_name(i) for i in leaves) == {"A", "B", "C", "D"}


# --- Invalid Newick string tests ---

def test_nwk2tree_empty_string(tmp_path):
    file = tmp_path / "empty.nwk"
    file.write_text("")
    with pytest.raises(RuntimeError):
        pcms.tree.nwk2tree(file)

def test_nwk2tree_unbalanced_parens(tmp_path):
    file = tmp_path / "bad_parens.nwk"
    file.write_text("(A,(B,C);")
    with pytest.raises(RuntimeError):
        pcms.tree.nwk2tree(file)

def test_nwk2tree_missing_semicolon(tmp_path):
    file = tmp_path / "no_semicolon.nwk"
    file.write_text("(A,B)")
    with pytest.raises(RuntimeError):
        pcms.tree.nwk2tree(file)

def test_nwk2tree_malformed_branch_lengths(tmp_path):
    file = tmp_path / "bad_lengths.nwk"
    file.write_text("(A:abc,B:0.1);")
    with pytest.raises(RuntimeError):
        pcms.tree.nwk2tree(file)
