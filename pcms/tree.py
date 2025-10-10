"""
@file tree.py
@brief A wrapper for the pcms._tree package.
@author Sean Svihla
"""

from typing import Optional, Union, List, Tuple
from functools import wraps

import numpy as np
from numpy.typing import NDArray, ArrayLike

import pcms._tree


def check_bounds(default=None, arg_names=("u",), constraints=None):
    """
    Decorator that checks whether one or more node indices are in bounds.

    Parameters
    ----------
    default: callable 
        Returns a default value derived from self.
    arg_names: tuple
        Argument names to bounds-check.
    constraints: dict 
        Map from arg_name to allowed extra values (e.g., -1 for v).
    """
    constraints = constraints or {}
    
    def decorator(method):
        @wraps(method)
        def wrapper(self, *args, **kwargs):
            args = list(args)
            for i, name in enumerate(arg_names):
                # Try positional first
                if i < len(args):
                    value = args[i]
                else:
                    value = kwargs.get(name, None)

                # Check for None and set default 
                if value is None:
                    if default is None:
                        raise ValueError(f"Node index '{name}' cannot be None")
                    value = default(self)
                    if i < len(args):
                        args[i] = value
                    else:
                        kwargs[name] = value
                # Check bounds
                else:
                    allowed = constraints.get(name, [])
                    if not (value in allowed or 0 <= value < self.n_nodes):
                        raise IndexError(f"Node index '{name}' = {value} out of bounds.")

            return method(self, *args, **kwargs)
        return wrapper
    return decorator


class Tree:
    """
    Wrapper class for `pcms._tree.Tree`.

    This class provides a Python interface to the underlying C++ implementation
    of a tree, enabling access to tree operations and properties from Python.
    """
    def __init__(self, n_nodes: int, print_label: Optional[str] = None) -> None:
        """
        Initialize a tree with the specified number of nodes.

        Parameters
        ----------
        n_nodes: int 
            The number of nodes in the tree.
        """
        if not isinstance(n_nodes, int):
            raise TypeError("Expected 'n_nodes' to have type 'int'")
        self._tree = pcms._tree.Tree(n_nodes)
        if print_label is None:
            self.print_label = "none"
        else:
            self.print_label = print_label

    @classmethod
    def _from_cpp_tree(cls, tree: pcms._tree.Tree) -> "Tree":
        obj = cls.__new__(cls)  # bypass __init__
        obj._tree = tree
        obj._print_label = "none"
        return obj

    def __del__(self) -> None:
        """Default destructor."""
        pass

    def __str__(self):
        return self._tree.to_string(self._print_label)

    def copy(self):
        """
        Returns
        -------
        pcms.tree.Tree
            A copy of the current tree.
        """
        return Tree._from_cpp_tree(pcms._tree.Tree(self._tree))

    def prune(self, keep: ArrayLike, is_leaf_idx: bool = True) -> "Tree":
        """
        Prunes leaves of the tree and returns a copy of the new tree. The
        original tree is not modified.

        Parameters
        ----------
        keep: NDArray
            The leaves to keep in the new tree.

        Returns
        -------
        pcms.tree.Tree
            The pruned tree.
        """
        keep = np.asarray(keep)
        leaves = self.find_leaves()
        if is_leaf_idx:
            keep = leaves[keep]
        if not np.all(np.isin(keep, leaves)):
            raise ValueError("Some values in 'keep' are not leaves.")
        keep = np.sort(keep)
        pruned_tree, inverse_postorder = self._tree.prune(keep)
        return Tree._from_cpp_tree(pruned_tree), inverse_postorder

    @property
    def print_label(self) -> str:
        """Return the label format used when printing the tree."""
        return self._print_label

    @print_label.setter
    def print_label(self, value: str) -> None:
        """
        Set the label format used when printing the tree.

        Parameters
        ----------
        value: str
            One of 'index', 'name', or 'none'.

        Raises
        -------
        ValueError
            If value is not one of 'index', 'name', or 'none'.
        """
        if value not in ["index", "name", "none"]:
            raise ValueError(f"Print label must be 'index', 'name', or 'none', not {value}")
        self._print_label = value

    @property
    def n_nodes(self) -> None:
        """
        Gets the number of nodes in the tree.

        Returns
        -------
        int
            The number of nodes in the tree.
        """
        return self._tree.get_n_nodes()

    @check_bounds() 
    def get_parent(self, u: int) -> int:
        """
        Gets the parent of a specified node.

        Parameters
        ----------
        u: int 
            Target node index.

        Returns
        -------
        int
            The index of the parent node.

        Raises
        -------
        IndexError
            If the node index is out of bounds.
        """
        return self._tree.get_parent(u)

    @check_bounds()
    def get_child(self, u: int) -> int:
        """
        Gets the child of a specified node.

        Parameters
        ----------
        u: int
            Target node index.

        Returns
        -------
        int
            The index of the child node.

        Raises
        -------
        IndexError
            If the node index is out of bounds.
        """
        return self._tree.get_child(u)

    @check_bounds()
    def get_sibling(self, u: int) -> int:
        """
        Gets the sibling of a specified node.

        Parameters
        ----------
        u: int
            Target node index.

        Returns
        -------
        int
            The index of the sibling node.

        Raises
        -------
        IndexError
            If the node index is out of bounds.
        """
        return self._tree.get_sibling(u)

    @check_bounds()
    def get_is_first(self, u: int) -> bool:
        """
        Gets whether the node is first among its siblings.

        Parameters
        ----------
        u: int
            Target node index.

        Returns
        -------
        bool
            The is_first parameter of the node.

        Raises
        -------
        IndexError
            If the node index is out of bounds.
        """
        return self._tree.get_is_first(u)

    @check_bounds(default=lambda self: None)
    def get_subtree_size(self, u: Optional[int] = None) -> Union[int, NDArray]:
        """
        Gets the size of the subtree rooted at a specified node.

        Parameters
        ----------
        u: int (optional, default None)
            Target node index.

        Returns
        -------
        int or np.ndarray
            The size of the subtree. If None is given, returns the entire
            subree_size vector.

        Raises
        -------
        IndexError
            If the node index is out of bounds.
        """
        return self._tree.get_subtree_size(u)

    @check_bounds(default=lambda self: None)
    def get_edge_length(self, u: Optional[int] = None) -> Union[float, NDArray]:
        """
        Gets the edge length of the edge connecting a node to its parent.

        Parameters
        ----------
        u: int (optional, default None)
            Target node index.

        Returns
        -------
        int or np.ndarray
            The edge length to the parent. If None is given, returns entire 
            edge_length vector.

        Raises
        -------
        IndexError
            If the node index is out of bounds.
        """
        return self._tree.get_edge_length(u)

    @check_bounds()
    def set_edge_length(self, u: int, value: float) -> None:
        """
        Sets the edge length of the edge connecting a node to its parent.

        Parameters
        ----------
        u: int
            Target node index.
        value: float
            The new edge length to set.

        Raises
        -------
        IndexError
            If the node index is out of bounds.
        """
        self._tree.set_edge_length(u, value)

    @check_bounds()
    def get_name(self, u: int) -> str:
        """
        Gets the name of the node.

        Parameters
        ----------
        u: int
            Target node index.

        Returns
        -------
        str
            The name of the node.

        Raises
        -------
        IndexError
            If the node index is out of bounds.
        """
        return self._tree.get_name(u)

    @check_bounds()
    def set_name(self, u: int, value: str) -> None:
        """
        Sets the name of the node.

        Parameters
        ----------
        u: int
            Target node index.
        value: str
            The new node name.

        Raises
        -------
        IndexError
            If the node index is out of bounds.
        """
        self._tree.set_name(u, value)

    @check_bounds(arg_names=("u", "v"), constraints={"v": [-1]})
    def link(self, u: int, v: int) -> None:
        """
        Links two nodes by making one the child of the other.

        Parameters
        ----------
        u: int
            Target node index.
        v: int
            Parent node.

        Raises
        -------
        IndexError
            If the node index is out of bounds.
        """
        if self.get_parent(u) != -1:
            raise RuntimeError("Cannot link a node before cutting.")
        self._tree.link(u, v)

    @check_bounds()
    def cut(self, u: int) -> None:
        """
        Cuts the link between a node and its parent.

        Parameters
        ----------
        u: int
            Target node index.

        Raises
        -------
        IndexError
            If the node index is out of bounds.
        """
        self._tree.cut(u)

    @check_bounds(arg_names=("u", "v"))
    def swap(self, u: int, v: int) -> None:
        """
        Swaps the positions of two nodes in the tree.

        Parameters
        ----------
        u: int
            First target node index.
        v: int
            Second target node index.

        Raises
        -------
        IndexError
            If the node index is out of bounds.
        """
        self._tree.swap(u, v)

    @check_bounds(arg_names=("u", "v"))
    def splice(self, u: int, v: int) -> None:
        """
        Splices (cuts and links) a node onto another.

        Parameters
        ----------
        u: int
            First target node index.
        v: int
            Second target node index.

        Raises
        -------
        IndexError
            If the node index is out of bounds.
        """
        self._tree.splice(u, v)

    @check_bounds(default=lambda self: self.find_root())
    def find_postorder(self, u: int) -> NDArray:
        """
        Parameters
        ----------
        u: int 
            The target node.

        Returns
        -------
        NDArray
            Node indices in postorder.
        """
        return self._tree.find_postorder(u)

    @check_bounds(default=lambda self: self.find_root())
    def find_mirror_postorder(self, u: int) -> NDArray:
        """
        Parameters
        ----------
        u: int 
            The target node.

        Returns
        -------
        NDArray
            Node indices in mirror postorder.
        """
        return self._tree.find_mirror_postorder(u)

    @check_bounds()
    def find_children(self, u: int) -> NDArray:
        """
        Finds the children of a specified node.

        Parameters
        ----------
        u: int
            Target node index.

        Returns
        -------
        np.ndarray
            A vector of child node indices.

        Raises
        -------
        IndexError
            If the node index is out of bounds.
        """
        return self._tree.find_children(u)

    @check_bounds()
    def find_ancestors(self, u: int) -> NDArray:
        """
        Finds the ancestors of a specified node.

        Parameters
        ----------
        u: int
            Target node index.

        Returns
        -------
        np.ndarray
            A vector of ancestor node indices.

        Raises
        -------
        IndexError
            If the node index is out of bounds.
        """
        return self._tree.find_ancestors(u)

    @check_bounds(arg_names=("u", "v"), constraints={"v": [-1]})
    def find_path(self, u: int, v: int) -> Tuple[NDArray, NDArray]:
        """
        Finds the path from one node to another.

        Parameters
        ----------
        u: int
            First target node index.
        v: int
            Second target node index.

        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            A pair of vectors: the first contains the ancestors of the first
            node, and the second contains the ancestors of the second node.

        Raises
        -------
        IndexError
            If the node index is out of bounds.
        """
        return self._tree.find_path(u, v)

    def find_root(self) -> int:
        """
        Finds the root of the tree.

        Returns
        -------
        int
            The index of the root node.
        """
        return self._tree.find_root()

    def find_is_planted(self) -> bool:
        """
        Finds whether the tree is planted.

        Returns
        -------
        bool
            Whether the tree is planted.
        """
        return self._tree.find_is_planted()

    @check_bounds(lambda self: self.find_root())
    def find_leaves(self, u: Optional[int] = None, return_depths = False) -> Union[NDArray, Tuple[NDArray, NDArray]]:
        """
        Finds the leaf nodes and their depths beneath a specified node.

        Parameters
        ----------
        u: int (optional, default None)
            The index of the node.
        return_depths: bool (optional, default False)
            Whether to return depths of leaves.

        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            A pair of vectors: the first contains the leaf node indices,
            and the second contains their corresponding depths. If None is
            given, the leaves of the entire tree are returned.

        Raises
        -------
        IndexError
            If the node index is out of bounds.
        """
        leaves, depths = self._tree.find_leaves(u)
        return (leaves, depths) if return_depths else leaves

    @check_bounds(lambda self: self.find_root())
    def find_interior_nodes(self, u: Optional[int] = None, return_depths = False) -> Union[NDArray, Tuple[NDArray, NDArray]]:
        """
        Finds the interior nodes beneath a specified node.

        Parameters
        ----------
        u: int (optional, default None)
            The index of the node.
        return_depths: bool (optional, default False)
            Whether to return depths of leaves.

        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            A pair of vectors: the first contains the interior node indices,
            and the second contains their corresponding depths. If None is
            given, the interior of the entire tree are returned.

        Raises
        -------
        IndexError
            If the node index if out of bounds.
        """
        interior_nodes, depths = self._tree.find_interior_nodes(u)
        return (interior_nodes, depths) if return_depths else interior_nodes

    def find_subtree_start_indices(self) -> NDArray:
        """
        Finds the leaf-indices (i.e., leaves indexed by postorder) of the 
        leftmost leaf in the subtree rooted at each node. Note that the 
        node-index (i.e., the indexing of all nodes in the tree) can be 
        obtained by referencing the output of find_leaves.

        Returns
        -------
        np.ndarray
            The subtree start index of each node in the tree.
        """
        return self._tree.find_subtree_start_indices()

    def find_n_leaves(self) -> int:
        """
        Find the number of leaves in the tree.

        Returns 
        -------
        int
            The number of leaves in the tree.
        """
        return self._tree.find_n_leaves()

    @check_bounds()
    def find_n_children(self, u: int) -> int:
        """
        Find the number of children of a node.

        Parameters
        ----------
        u: int 
            Target node index.

        Returns
        -------
        int
            The number of children of the node.
        """
        return self._tree.find_n_children(u)
    
    def find_n_wavelets(self) -> int:
        """
        Find the number of wavelets in the tree's Haar-like basis.

        Returns
        -------
        int 
            The number of wavelets in the basis
        """
        return self._tree.find_n_wavelets()

    def find_epl(self) -> int:
        """
        Finds the expected path length (EPL) of the tree.

        Returns
        -------
        int 
            The expected path length.

        Raises
        -------
        IndexError
            If the node index is out of bounds.
        """
        return self._tree.find_epl()

    @check_bounds(default=lambda self: None, arg_names=("u", "v"))
    def find_tbl(self, u: Optional[int] = None, v: Optional[int] = None) -> Union[float, NDArray]:
        """
        Finds the trace length of the edge above a node, or the unique
        shortest path between two nodes, or all nodes.

        Parameters
        ----------
        u: int (optional, default None)
            Target node index or first node.
        v: int (optional, default None)
            Optional second node index.

        Returns
        -------
        float or np.ndarray
            If one argument: the trace length of the edge to the parent.
            If two arguments: the trace length of the path between u and v.
            If none: a vector of the trace length of each node.

        Raises
        -------
        IndexError
            If the node index is out of bounds.
        ValueError
            If parameter 'v' given when 'u' is None.
        """
        if u is None and v is not None:
            u = v
        return self._tree.find_tbl(u, v)

    @check_bounds()
    def compute_wavelets(self, u: int) -> Union[NDArray, List[NDArray]]:
        """
        Constructs wavelets associated to an interior node.

        Parameters
        ----------
        u: int
            Target node index.

        Returns
        -------
        List[np.ndarray]
            A list of NumPy arrays, each containing the values of a wavelet
            associated with the target node index.

        Raises
        -------
        IndexError
            If the node index is out of bounds.
        """
        if self.get_child(u) == -1:
            raise ValueError("Node must be an interior node, not a leaf.")
        wavelets = self._tree.compute_wavelets(u)
        return wavelets if len(wavelets) > 1 else wavelets[0]

    @check_bounds()
    def compute_supports(self, u: int) -> Union[NDArray, List[NDArray]]:
        """
        Finds supports of wavelets associated to an interior node.

        Note:
            The first returned vector contains a list of the indices of
            all nodes in the support of all wavelets associated with the
            target node. The second vector contains the boundary indices
            of a partition into each support.

            For example, the ith wavelet associated with 'u' takes on its 
            left-value on leaves partition[i] to partition[i+1].

        Parameters
        ----------
        u: int
            Target node index.

        Returns
        -------
        List[np.ndarray]
            A pair of vectors: the first containing the indices of leaf
            nodes in the support, the second the boundaries between 
            partitions of the support.

        Raises
        -------
        IndexError
            If the node index is out of bounds.
        ValueError
            If the input node is a leaf node.
        """
        if self.get_child(u) == -1:
            raise ValueError("Input 'u' cannot be a leaf node")
        supports = self._tree.compute_supports(u)
        return supports if len(supports) > 1 else supports[0]


def nwk2tree(filename: str, ensure_planted: bool = False) -> Tree:
    """
    Construct a Tree object from a Newick format string.

    Parameters
    ----------
    newick_string: str
        The Newick string.
    ensure_planted: bool (optional, default False)
        Whether to ensure that the tree is planted (if not already planted in
        the Newick string).

    Returns
    -------
    Tree
        A Tree object corresponding to the Newick string.
    """
    with open(filename) as f:
        nwk_str = f.read()
    nwk_str = nwk_str.replace("\n","")
    if nwk_str == "":
        raise RuntimeError(f"Newick string in {filename} is empty")
    if not nwk_str.endswith(';'):
        raise RuntimeError(f"Newick string in {filename} must end with ';'")
    return Tree._from_cpp_tree(pcms._tree.nwk2tree(nwk_str, ensure_planted))


def make_buffer(n_trees: int, n_nodes: int):
    return [pcms._tree.Tree(n_nodes) for _ in range(n_trees)]


def _random_tree_generator(
    sampler: callable,
    batched_sampler: callable,
    n_leaves: int,
    planted: bool, 
    n_samples: int, 
    seed: Optional[int] = None,
    buffer: Optional[List[Tree]] = None,
):
    if n_leaves < 2:
        raise ValueError("Required n_leaves >= 2.")
    n_nodes = 2 * n_leaves if planted else 2 * n_leaves - 1
    n_nodes = int(n_nodes)
    if n_samples > 1:
        buffer = make_buffer(n_samples, n_nodes) if buffer is None else buffer
        batched_sampler(buffer, planted, n_samples, seed)
        return [Tree._from_cpp_tree(t_) for t_ in buffer]
    elif n_samples == 1:
        t = pcms.tree.Tree(n_nodes)
        sampler(t._tree, planted, seed)
        return t
    else:
        raise ValueError("n_samples should be at least one.")


def remy(
    n_leaves: int,
    planted: bool = True, 
    n_samples: int = 1, 
    seed: Optional[int] = None,
    buffer: Optional[List[Tree]] = None
) -> Tree:
    """
    Construct a random Tree object from the uniform distribution.

    Parameters
    ----------
    n_leaves: int 
        Number of leaves in the tree.
    planted: bool
        Whether the tree should be planted.
    n_samples: int
        Number of independent samples to generate.
    seed: int
        A seed for random generation.
    buffer: List[Tree] (optional, default None)
        A buffer for placing generated trees.

    Returns
    -------
    Tree | List[Tree]
        Randomly generated Tree object(s)

    Raises
    -------
    ValueError
        If n_leaves < 2.
    """
    return _random_tree_generator(pcms._tree.remy, pcms._tree.remy_batched,
                                  n_leaves, planted, n_samples, seed, buffer)


def cbst(
    n_leaves: int, 
    planted: bool = True, 
    n_samples: int = 1, 
    seed: Optional[int] = None,
    buffer: Optional[List[Tree]] = None
) -> Tree:
    """
    Construct a random Tree object from the critical beta-splitting
    distribution.

    Parameters
    ----------
    n_leaves: int 
        Number of leaves in the tree.
    planted: bool
        Whether the tree should be planted.
    n_samples: int
        Number of independent samples to generate.
    seed: int
        A seed for random generation.
    buffer: List[Tree] (optional, default None)
        A buffer for placing generated trees.

    Returns
    -------
    Tree | List[Tree]
        Randomly generated Tree object(s)

    Raises
    -------
    ValueError
        If n_leaves < 2.
    """
    return _random_tree_generator(pcms._tree.cbst, pcms._tree.cbst_batched,
                                  n_leaves, planted, n_samples, seed, buffer)