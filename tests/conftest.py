# tests/conftest.py

import pytest
import pcms.tree


@pytest.fixture
def singleton_tree():
    return pcms.tree.Tree(1)

@pytest.fixture
def small_structured_tree():
    tree = pcms.tree.Tree(4)
    tree.link(1, 2)
    tree.link(0, 2)
    tree.link(2, 3)
    return tree

@pytest.fixture
def small_unstructured_tree():
    return pcms.tree.Tree(3)

@pytest.fixture
def large_structured_tree():
    tree = pcms.tree.Tree(5)
    tree.link(1, 2)
    tree.link(0, 2)
    tree.link(3, 4)
    tree.link(2, 4)
    return tree

@pytest.fixture
def small_structured_ternary_tree():
    tree = pcms.tree.Tree(5)
    tree.link(2, 3)
    tree.link(1, 3)
    tree.link(0, 3)
    tree.link(3, 4)
    return tree

@pytest.fixture
def uniform_random_tree():
    return pcms.tree.remy(5, seed=42)

@pytest.fixture
def cbst_random_tree():
    return pcms.tree.cbst(5, seed=123)
