/**
 * @file tree.hpp
 * @brief Header file for the Tree class and related operations.
 * 
 * This file defines the `TreeNode` structure and the `Tree` class, which represents
 * a tree data structure supporting various operations such as linking, cutting, 
 * swapping subtrees, and finding paths, ancestors, and children. It also provides
 * functions for managing edge lengths and subtree sizes, and for printing the tree structure.
 * 
 * The `Tree` class is used to model hierarchical relationships, where each node can 
 * have a parent, children, and siblings. The tree is implemented as a static structure, 
 * with a fixed number of nodes.
 * 
 * Functions are provided to manipulate the tree structure, including linking nodes, 
 * cutting edges, and swapping subtrees. Additionally, utility functions for traversing
 * the tree, finding leaves, and computing the effective path length (EPL) are included.
 * 
 * @author Sean Svihla
 */
#ifndef TREE_HPP
#define TREE_HPP

// standard library includes
#include <memory>
#include <string>
#include <vector>

// Pybind11 includes
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

/**
 * @struct TreeTopology
 * @brief A structure to represent parent-child relationships in the tree.
 */
struct TreeTopology {
    int parent = -1;         ///< Parent node index (-1 if no parent)
    int child = -1;          ///< Child node index (-1 if no child)
    int sibling = -1;        ///< Sibling node index (-1 if no sibling)
    bool is_first = false;   ///< Whether node is the first (leftmost) sibling
};

/**
 * @struct TreeNumerics
 * @brief A structure to represent numerical tree attributes.
 */
struct TreeNumerics {
    py::array_t<int> subtree_size;       ///< Number of children in the subtree
    py::array_t<double> edge_length;     ///< Edge length to parent

    TreeNumerics(size_t n_nodes) 
    : subtree_size(n_nodes), edge_length(n_nodes) 
    {
        auto subtree_size_ = subtree_size.mutable_unchecked<1>();
        for(py::ssize_t i = 0; i < static_cast<py::ssize_t>(n_nodes); ++i)
        {
            subtree_size_(i) = 1;
        }
        auto edge_length_ = edge_length.mutable_unchecked<1>();
        for(py::ssize_t i = 0; i < static_cast<py::ssize_t>(n_nodes); ++i)
        {
            edge_length_(i) = 0.0;
        }
    } 
};

/**
 * @class Tree
 * @brief A class to represent a tree data structure.
 *
 * The tree is implemented as a set of nodes, where each node contains
 * information about its parent, child, sibling, and subtree size.
 * The class provides various methods for manipulating and querying the tree.
 */
class Tree {
public:
    /**
     * @brief Constructs a tree with the specified size.
     * 
     * Initializes the tree with a given number of nodes. The tree is
     * initialized without parent-child-sibling relationships.
     * 
     * @param n_nodes The number of nodes in the tree.
     */
    Tree(int n_nodes);

    /**
     * @brief Default destructor.
     */
    ~Tree();

    /**
     * @brief Gets the number of nodes in the tree.
     * @return The number of nodes in the tree.
     */
    int get_n_nodes() const;

    /**
     * @brief Gets the parent of a specified node.
     * @param u Target node index.
     * @return The index of the parent node.
     */
    int get_parent(int u) const;

    /**
     * @brief Gets the child of a specified node.
     * @param u Target node index.
     * @return The index of the child node.
     */
    int get_child(int u) const;

    /**
     * @brief Gets the sibling of a specified node.
     * @param u Target node index.
     * @return The index of the sibling node.
     */
    int get_sibling(int u) const;

    /**
     * @brief Gets whether the node is first among its siblings.
     * @param u Target node index.
     * @return The is_first parameter of the node.
     */
    bool get_is_first(int u) const;

    /**
     * @brief Gets the size of the subtree rooted at a specified node.
     * @param u Target node index.
     * @return The size of the subtree.
     */
    int get_subtree_size(int u) const;

    /**
     * @brief Gets the entire subtree size vector for the tree.
     * @return The subtree size vector.
     */
    const py::array_t<int>& get_subtree_size() const;

    /**
     * @brief Gets the edge length of the edge connecting a node to its parent.
     * @param u Target node index.
     * @return The edge length to the parent.
     */
    double get_edge_length(int u) const;

    /**
     * @brief Gets the entire edge length vector of the tree.alignas
     * @return The edge length vector.
     */
    const py::array_t<double>& get_edge_length() const;

    /**
     * @brief Sets the edge length of the edge connecting a node to its parent.
     * @param u Target node index.
     * @param value The new edge length to set.
     */
    void set_edge_length(int u, double value);

    /**
     * @brief Gets the name of the node.
     * @param u Target node index.
     * @return The namae of the node. 
     */
    std::string get_name(int u) const;

    /**
     * @brief Sets the name of the node.
     * @param u Target node index.
     * @param value The new node name.
     */
    void set_name(int u, const std::string& value);

    /**
     * @brief Links two nodes by making one the child of the other.
     * @param u Target node index.
     * @param v Parent node.
     */
    void link(int u, int v);

    /**
     * @brief Cuts the link between a node and its parent.
     * @param u Target node index.
     */
    void cut(int u);

    /**
     * @brief Swaps the positions of two nodes in the tree.
     * @param u First target node index.
     * @param v Second target node index.
     */
    void swap(int u, int v);

    /**
     * @brief Finds the children of a specified node.
     * @param u Target node index.
     * @return A vector of child node indices.
     * @throws py::index_error If the node index is out of bounds.
     */
    py::array_t<int> find_children(int u) const;

    /**
     * @brief Finds the ancestors of a specified node.
     * @param u Target node index.
     * @return A vector of ancestor node indices.
     */
    py::array_t<int> find_ancestors(int u) const;

    /**
     * @brief Finds the path from one node to another.
     * @param u First target node index.
     * @param v Second target node index.
     * @return A pair of vectors: the first contains the ancestors of the first
     *         node, and the second contains the ancestors of the second node.
     */
    std::pair<py::array_t<int>, py::array_t<int>> find_path(int u, int v) const;

    /**
     * @brief Finds the root of the tree.
     * @return The index of the root node.
     */
    int find_root() const;

    /**
     * @brief Finds the leafs nodes and their depths beneath a specified node.
     * @param u The index of the node.
     * @return A pair of vectors: the first contains the leaf node indices,
     *         and the second contains their corresponding depths.
     */
    std::pair<py::array_t<int>, py::array_t<int>> find_leaves(int u) const;

    /**
     * @brief Computes the number of leaves in the tree. 
     * @return The number of leaves in the tree.
     */
    int find_n_leaves() const;

    /**
     * @brief Finds the expected path length (EPL) of the tree.
     * @return The expected path length.
     */
    int find_epl() const;

    /**
     * @brief Finds the trace length of the edge above the node.
     * @param u Target node index.
     * @return The trace length of the edge to the parent. 
     */
    double find_tbl(int u) const;

    /**
     * @brief Finds the trace length of the unique shortest path between
     *        two nodes.
     * @param u First taget node index.
     * @param v Second target node index.
     * @return The trace length of the path between u and v.
     */
    double find_tbl(int u, int v) const;

    /**
     * @brief Finds the trace length of all nodes in the tree.
     * @return A vector of the trace length of each node.
     */
    py::array_t<double> find_tbl() const;

    /**
     * @brief Constructs wavelets associated to an interior node.
     * @param u Target node index.
     * @return A list of NumPy arrays, each containing the values of a wavelet
     *         associated with the Target node index.
     */
    py::list compute_wavelets(int u) const;

    /**
     * @brief 
     * @return 
     */
    int find_nnz_max() const;

    /**
     * @brief Finds supports of wavelets associated to an interior node. 
     * 
     * Note that the first returned vector contains a list of the indices of 
     * all nodes in the support of all wavelets associated with the target 
     * node. The second vector contains the boundary indices of a partition
     * into each support.
     * 
     * For example, the the ith wavelet associated with 'u' takes on its 
     * left-value on leaves partition[i] to partition[i+1].
     * 
     * @param u Target node index.
     * @return A pair of vectors: the first containing the indices of leaf
     *         nodes in the support, the second the boundaries between 
     *         partitions of the support.
     */
    py::list compute_supports(int u) const;

    /**
     * @brief Prints the tree in an ASCII-art format.
     * 
     * This function starts from the root and prints the tree structure using
     * ASCII characters to represent the hierarchy of nodes.
     * 
     * @param label The labelling style (none, index, or name).
     */
    std::string to_string(const std::string& label) const;

private:
    int n_nodes;                        ///< Number of nodes in the tree.
    std::vector<TreeTopology> topology; ///< Parent-child-sibling relationships
    TreeNumerics numerics;              ///< Edge length and subtree sizes
    std::vector<std::string> names;     ///< Array of node names.

    /**
     * @brief Private parent node setter. 
     * @param u Target node index.
     * @param v Parent node index.
     */
    void set_parent(int u, int v);

    /**
     * @brief Private child node setter. 
     * @param u Target node index.
     * @param v Child node index.
     */
    void set_child(int u, int v);

    /**
     * @brief Private sibling node setter. 
     * @param u Target node index.
     * @param v Sibling node index.
     */
    void set_sibling(int u, int v);

    /**
     * @brief Private is_first setter.  
     * @param u Target node index. 
     * @param value
     */
    void set_is_first(int u, bool value);

    /**
     * @brief Private subtree_size setter. 
     * @param u Target node index.  
     * @param value
     */
    void set_subtree_size(int u, int value);

    /**
     * @brief Private find ancestors method.  
     * @param u Target node index. 
     * @return A vector of ancestor node indices.
     */
    std::vector<int> find_ancestors_(int u) const;

    /**
     * @brief Helper function to print a specific node and its children.
     * @param node The index of the node to print.
     * @param prefix The prefix string used to format the tree output.
     * @param is_last Whether this node is the last child of its parent.
     */
    std::string to_string_(int node, const std::string& prefix, bool is_last, 
                           const std::string& label) const;
};

// inlined class methods
inline int 
Tree::get_n_nodes() const
{
    return n_nodes;
}

inline int 
Tree::get_parent(int u) const
{
    return topology[u].parent;
}

inline void
Tree::set_parent(int u, int v)
{
    topology[u].parent = v;
}

inline int 
Tree::get_child(int u) const
{
    return topology[u].child;
}

inline void
Tree::set_child(int u, int v)
{
    topology[u].child = v;
}

inline int 
Tree::get_sibling(int u) const
{
    return topology[u].sibling;
}

inline void
Tree::set_sibling(int u, int v)
{
    topology[u].sibling = v;
}

inline bool
Tree::get_is_first(int u) const
{
    return topology[u].is_first;
}

inline void
Tree::set_is_first(int u, bool value)
{
    topology[u].is_first = value;
}

inline int 
Tree::get_subtree_size(int u) const
{
    auto subtree_size_ = numerics.subtree_size.unchecked<1>();
    return subtree_size_[u];
}

inline const py::array_t<int>&
Tree::get_subtree_size() const
{
    return numerics.subtree_size;
}

inline void
Tree::set_subtree_size(int u, int value)
{
    auto subtree_size_ = numerics.subtree_size.mutable_unchecked<1>();
    subtree_size_(u) = value;
}

inline double 
Tree::get_edge_length(int u) const
{
    auto edge_length_ = numerics.edge_length.unchecked<1>();
    return edge_length_[u];
}

inline const py::array_t<double>&
Tree::get_edge_length() const 
{
    return numerics.edge_length;
}

inline void 
Tree::set_edge_length(int u, double value)
{
    auto edge_length_ = numerics.edge_length.mutable_unchecked<1>();
    edge_length_(u) = value;
}

inline std::string
Tree::get_name(int u) const
{
    return names[u];
}

inline void
Tree::set_name(int u, const std::string& value)
{
    names[u] = value;
}

inline double 
Tree::find_tbl(int u) const 
{
    return get_subtree_size(u) * get_edge_length(u);
}

/**
 * @brief Construct a Tree object from a Newick format string.
 * @param newick_string The Newick string.
 * @return A Tree object corresponding to the Newick stirng.
 */
Tree *nwk2tree(const std::string& newick_string);

#endif // TREE_HPP