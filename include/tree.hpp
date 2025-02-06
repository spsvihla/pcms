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
 * The `remy` function generates a random tree with a specified number of leaves.
 * 
 * This header is part of a larger project aimed at managing and manipulating tree data structures.
 * 
 * @author Sean Svihla
 */
#ifndef TREE_HPP
#define TREE_HPP

// Standard library includes
#include <vector>

/**
 * @struct TreeNode
 * @brief A structure to represent a node in the tree.
 *
 * Each node stores information about its parent, child, sibling,
 * its subtree size, and the edge length to its parent.
 */
struct TreeNode {
    int parent;        ///< Parent node index (-1 if no parent)
    int child;         ///< Child node index (-1 if no child)
    int sibling;       ///< Sibling node index (-1 if no sibling)
    int subtree_size;  ///< Size of the subtree rooted at this node
    float edge_length; ///< Length of the edge to the parent node
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
     * @throws py::value_error If the node index is out of bounds.
     */
    int get_parent(int u) const;

    /**
     * @brief Gets the child of a specified node.
     * 
     * @param u The index of the node.
     * @return The index of the child node.
     * @throws py::value_error If the node index is out of bounds.
     */
    int get_child(int u) const;

    /**
     * @brief Gets the sibling of a specified node.
     * 
     * @param u The index of the node.
     * @return The index of the sibling node.
     * @throws py::value_error If the node index is out of bounds.
     */
    int get_sibling(int u) const;

    /**
     * @brief Gets the size of the subtree rooted at a specified node.
     * 
     * @param u The index of the node.
     * @return The size of the subtree.
     * @throws py::value_error If the node index is out of bounds.
     */
    int get_subtree_size(int u) const;

    /**
     * @brief Gets the edge length of the edge connecting a node to its parent.
     * 
     * @param u The index of the node.
     * @return The edge length to the parent.
     * @throws py::value_error If the node index is out of bounds.
     */
    float get_edge_length(int u) const;

    /**
     * @brief Sets the edge length of the edge connecting a node to its parent.
     * 
     * @param u The index of the node.
     * @param value The new edge length to set.
     * @throws py::value_error If the node index is out of bounds.
     */
    void set_edge_length(int u, float value);

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
     * @throws py::value_error If the node index is out of bounds.
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
     * @throws py::value_error If the node index is out of bounds.
     */
    std::vector<int> find_children(int u) const;

    /**
     * @brief Finds the ancestors of a specified node.
     * 
     * @param u The index of the node.
     * @return A vector of ancestor node indices.
     * @throws py::value_error If the node index is out of bounds.
     */
    std::vector<int> find_ancestors(int u) const;

    /**
     * @brief Finds the support nodes and their depths for a specified node.
     * 
     * @param u The index of the node.
     * @return A pair of vectors: the first contains support node indices,
     *         and the second contains their corresponding depths.
     * @throws py::value_error If the node index is out of bounds.
     */
    std::pair<std::vector<int>, std::vector<int>> find_support(int u) const;

    /**
     * @brief Finds the path from one node to another.
     * 
     * @param u The index of the first node.
     * @param v The index of the second node.
     * @return A pair of vectors: the first contains the ancestors of the first
     *         node, and the second contains the ancestors of the second node.
     * @throws py::value_error If any of the node indices are out of bounds.
     */
    std::pair<std::vector<int>, std::vector<int>> find_path(int u, int v) const;

    /**
     * @brief Finds the root of the tree.
     * 
     * @return The index of the root node.
     */
    int find_root() const;

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
     * @brief Prints the tree in an ASCII-art format.
     * 
     * This function starts from the root and prints the tree structure using
     * ASCII characters to represent the hierarchy of nodes.
     */
    void print() const;

private:
    int n_nodes;     ///< Number of nodes in the tree.
    std::vector<TreeNode> nodes; ///< Array of tree nodes.

    /**
     * @brief Helper function to print a specific node and its children.
     * 
     * @param node The index of the node to print.
     * @param prefix The prefix string used to format the tree output.
     * @param is_last Whether this node is the last child of its parent.
     */
    void print_node(int node, const std::string& prefix, bool is_last) const;
};

#endif // TREE_HPP