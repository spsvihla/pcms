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
    std::vector<int> subtree_size;       ///< Number of children in the subtree
    std::vector<double> edge_length;     ///< Edge length to parent

    TreeNumerics(size_t n_nodes) 
    : subtree_size(n_nodes, 1), edge_length(n_nodes, 0.0) {} 
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
     * Initializes the tree with a given number of nodes, setting the initial
     * relationships between the nodes to default values.
     * 
     * @param size The number of nodes in the tree.
     */
    Tree(int size);

    /**
     * @brief Default destructor.
     */
    ~Tree();

    /**
     * @brief Gets the size of the tree.
     * 
     * @return The number of nodes in the tree.
     */
    int get_size() const;

    /**
     * @brief Gets the parent of a specified node.
     * 
     * @param u The index of the node.
     * @return The index of the parent node.
     * @throws py::index_error If the node index is out of bounds.
     */
    int get_parent(int u) const;

    /**
     * @brief Gets the child of a specified node.
     * 
     * @param u The index of the node.
     * @return The index of the child node.
     * @throws py::index_error If the node index is out of bounds.
     */
    int get_child(int u) const;

    /**
     * @brief Gets the sibling of a specified node.
     * 
     * @param u The index of the node.
     * @return The index of the sibling node.
     * @throws py::index_error If the node index is out of bounds.
     */
    int get_sibling(int u) const;

    /**
     * @brief Gets whether the node is first among its siblings.
     * 
     * @param u The index of the node.
     * @return The is_first parameter of the node.
     * @throws py::index_error If the node index is out of bounds.
     */
    bool is_first(int u) const;

    /**
     * @brief Gets the size of the subtree rooted at a specified node.
     * 
     * @param u The index of the node.
     * @return The size of the subtree.
     * @throws py::index_error If the node index is out of bounds.
     */
    int get_subtree_size(int u) const;

    /**
     * @brief Gets the entire subtree size vector for the tree.
     * 
     * @return The subtree size vector.
     */
    const std::vector<int>& get_subtree_size() const;

    /**
     * @brief Gets the edge length of the edge connecting a node to its parent.
     * 
     * @param u The index of the node.
     * @return The edge length to the parent.
     * @throws py::index_error If the node index is out of bounds.
     */
    double get_edge_length(int u) const;

    /**
     * @brief Gets the entire edge length vector of the tree.alignas
     * 
     * @return The edge length vector.
     */
    const std::vector<double>& get_edge_length() const;

    /**
     * @brief Sets the edge length of the edge connecting a node to its parent.
     * 
     * @param u The index of the node.
     * @param value The new edge length to set.
     * @throws py::index_error If the node index is out of bounds.
     */
    void set_edge_length(int u, double value);

    /**
     * @brief Gets the name of the node.
     * 
     * @param u The index of the node.
     * @return The namae of the node. 
     * @throws py::index_error If the node index is out of bounds.
     */
    std::string get_name(int u) const;

    /**
     * @brief Sets the name of the node.
     * 
     * @param u The index of the node.
     * @param value The new node name.
     * @throws py::index_error If the node index is out of bounds.
     */
    void set_name(int u, const std::string& value);

    /**
     * @brief Links two nodes by making one the child of the other.
     * 
     * @param u The child node index.
     * @param v The parent node index.
     * @throws py::value_error If the node indices are out of bounds or if
     *         the link cannot be made.
     */
    void link(int u, int v);

    /**
     * @brief Cuts the link between a node and its parent.
     * 
     * @param u The index of the node to cut from its parent.
     * @throws py::index_error If the node index is out of bounds.
     */
    void cut(int u);

    /**
     * @brief Swaps the positions of two nodes in the tree.
     * 
     * @param u The index of the first node.
     * @param v The index of the second node.
     */
    void swap(int u, int v);

    /**
     * @brief Finds the children of a specified node.
     * 
     * @param u The index of the node.
     * @return A vector of child node indices.
     * @throws py::index_error If the node index is out of bounds.
     */
    std::vector<int> find_children(int u) const;

    /**
     * @brief Finds the ancestors of a specified node.
     * 
     * @param u The index of the node.
     * @return A vector of ancestor node indices.
     * @throws py::index_error If the node index is out of bounds.
     */
    std::vector<int> find_ancestors(int u) const;

    /**
     * @brief Finds the path from one node to another.
     * 
     * @param u The index of the first node.
     * @param v The index of the second node.
     * @return A pair of vectors: the first contains the ancestors of the first
     *         node, and the second contains the ancestors of the second node.
     * @throws py::index_error If any of the node indices are out of bounds.
     */
    std::pair<std::vector<int>, std::vector<int>> find_path(int u, int v) const;

    /**
     * @brief Finds the root of the tree.
     * 
     * @return The index of the root node.
     */
    int find_root() const;

    /**
     * @brief Finds the leafs nodes and their depths beneath a specified node.
     * 
     * @param u The index of the node.
     * @return A pair of vectors: the first contains the leaf node indices,
     *         and the second contains their corresponding depths.
     * @throws py::index_error If the node index is out of bounds.
     */
    std::pair<std::vector<int>, std::vector<int>> find_leaves(int u) const;

    /**
     * @brief Finds the leaves of the tree.
     * 
     * @return A pair of vectors: the first contains the leaf node indices,
     *         and the second contains their corresponding depths.
     */
    std::pair<std::vector<int>, std::vector<int>> find_leaves() const;

    /**
     * @brief Finds the expected path length (EPL) of the tree.
     * 
     * @return The expected path length.
     */
    int find_epl() const;

    /**
     * @brief Finds the trace length of the edge above the node.
     * 
     * @param u
     * @return The trace length of the edge to the parent. 
     * @throws py::index_error If the node index is out of bounds.
     */
    double find_tbl(int u) const;

    /**
     * @brief Finds the trace length of the unique shortest path between
     *        two nodes.
     * 
     * @param u
     * @param v
     * @return The trace length of the path between u and v.
     * @throws py::index_error If the node index is out of bounds.
     */
    double find_tbl(int u, int v) const;

    /**
     * @brief Finds the trace length of all nodes in the tree.
     * 
     * @return A vector of the trace length of each node.
     */
    std::vector<double> find_tbl() const;

    /**
     * @brief Constructs wavelets associated to an interior node.
     * 
     * Note that the wavelets are only constructed for their values. The
     * support of each wavelet is not computed, but may be computed by 
     * calling Tree::find_supports.
     * 
     * @param u
     * @return 
     * @throws py::index_error If the node index is out of bounds or a leaf.
     */
    std::vector<std::vector<double>> find_wavelets(int u) const;

    /**
     * @brief Finds supports of wavelets associated to an interior node. 
     * 
     * @param u
     * @return 
     * @throws py::index_error If the node index is out of bounds or a leaf.
     */
    std::pair<std::vector<int>, std::vector<int>> find_supports(int u) const;

    /**
     * @brief
     * 
     * @return
     */
    int find_nnz() const;

    /**
     * @brief Prints the tree in an ASCII-art format.
     * 
     * This function starts from the root and prints the tree structure using
     * ASCII characters to represent the hierarchy of nodes.
     * 
     * @param label The labelling style (none, index, or name).
     */
    void print(const std::string& label) const;

private:
    int n_nodes;                        ///< Number of nodes in the tree.
    std::vector<TreeTopology> topology; ///< Parent-child-sibling relationships
    TreeNumerics numerics;              ///< Edge length and subtree sizes
    std::vector<std::string> names;     ///< Array of node names.

    /**
     * @brief Helper function to print a specific node and its children.
     * 
     * @param node The index of the node to print.
     * @param prefix The prefix string used to format the tree output.
     * @param is_last Whether this node is the last child of its parent.
     */
    void print_node(int node, const std::string& prefix, bool is_last, 
                    const std::string& label) const;
};

Tree *nwk2tree(const std::string& filename);

#endif // TREE_HPP