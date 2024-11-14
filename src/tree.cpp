/**
 * @file tree_node.c
 * 
 * Contains the implementation of the TreeNode class and related functions.
 * This file handles tree node operations such as swapping subtrees and computing 
 * Lowest Common Ancestors (LCAs) in a tree structure.
 * 
 * @author Sean Svihla
 */

// Standard library includes
#include <iostream>
#include <numeric>
#include <random>
#include <stack>
#include <vector>

// Project-specific includes
#include "tree.hpp"

// Pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;


// constructor
Tree::Tree(int n_nodes)
: n_nodes(n_nodes), nodes(n_nodes)
{
    for(int i = 0; i < n_nodes; i++)
    {
        nodes.at(i).parent = -1;
        nodes.at(i).child = -1;
        nodes.at(i).sibling = -1;
        nodes.at(i).subtree_size = 1;
        nodes.at(i).edge_length = 1.0f;
    }
}

int 
Tree::get_size() const
{
    return n_nodes;
}

int 
Tree::get_parent(int u) const
{
    return nodes.at(u).parent;
}

int 
Tree::get_child(int u) const
{
    return nodes.at(u).child;
}

int 
Tree::get_sibling(int u) const
{
    return nodes.at(u).sibling;
}

int 
Tree::get_subtree_size(int u) const
{
    return nodes.at(u).subtree_size;
}

float 
Tree::get_edge_length(int u) const
{
    return nodes.at(u).edge_length;
}

void 
Tree::set_edge_length(int u, float value)
{
    nodes.at(u).edge_length = value;
}

void 
Tree::link(int u, int v)
{
    assert(nodes.at(u).parent == -1);

    if(v == -1)
    {
        return;
    }

    int is_leaf = nodes.at(v).child == -1 ? 1 : 0;

    // update pointers
    nodes.at(u).parent = v;
    nodes.at(u).sibling = nodes.at(v).child;
    nodes.at(v).child = u;

    // update subtree_size
    int diff = nodes.at(u).subtree_size - is_leaf;
    for(; v != -1; v = nodes.at(v).parent)
    {
        nodes.at(v).subtree_size += diff;
    }
}

void 
Tree::cut(int u)
{
    int p = nodes.at(u).parent;
    if(p == -1)
    {
        // node is an orphan
        return;
    }

    // update pointers
    if(nodes.at(p).child == u)
    {
        nodes.at(p).child = nodes.at(u).sibling;
    }
    else
    {
        int c = nodes.at(p).child;
        while(nodes.at(c).sibling != u)
        {
            c = nodes.at(c).sibling;
        }
        nodes.at(c).sibling = nodes.at(u).sibling;
    }
    int is_leaf = nodes.at(p).child == -1 ? 1 : 0;

    // update subtree_size
    int diff = is_leaf - nodes.at(u).subtree_size;
    for(; p != -1; p = nodes.at(p).parent)
    {
        nodes.at(p).subtree_size += diff;
    }

    nodes.at(u).parent = -1;
    nodes.at(u).sibling = -1;
}

void 
Tree::swap(int u, int v)
{
    int pu = nodes.at(u).parent;
    int pv = nodes.at(v).parent;
    cut(u);
    cut(v);
    link(u, pv);
    link(v, pu);
}

std::vector<int> 
Tree::find_children(int u) const
{
    std::vector<int> children;
    for(int c = nodes.at(u).child; c != -1; c = nodes.at(c).sibling)
    {
        children.push_back(c);
    }
    return children;
}

std::vector<int> 
Tree::find_ancestors(int u) const
{
    std::vector<int> ancestors;
    for(u = nodes.at(u).parent; u != -1; u = nodes.at(u).parent)
    {
        ancestors.push_back(u);
    }
    return ancestors;
}

std::pair<std::vector<int>, std::vector<int>> 
Tree::find_support(int u) const
{
    std::vector<int> support;
    std::vector<int> depths;
    std::stack<std::pair<int, int>> stack;
    stack.push(std::make_pair(u, 0));

    while(!stack.empty())
    {
        int v = stack.top().first;
        int depth = stack.top().second;
        stack.pop();

        if(nodes.at(v).child == -1)
        {
            support.push_back(v);     // leaf
            depths.push_back(depth);  // depth
        }

        for(int c = nodes.at(v).child; c != -1; c = nodes.at(c).sibling)
        {
            stack.push({c, depth + 1});
        }
    }
    return std::make_pair(support, depths);
}

std::pair<std::vector<int>,std::vector<int>> 
Tree::find_path(int u, int v) const
{
    std::vector<int> u_path = find_ancestors(u);
    std::vector<int> v_path = find_ancestors(v);
    while(!u_path.empty() && !v_path.empty()) 
    {
        if(u_path.back() != v_path.back())
        {
            break;
        }
        u_path.pop_back();
        v_path.pop_back();
    }
    return make_pair(u_path, v_path);
}

int
Tree::find_root() const
{
    return find_ancestors(0).back();
}

std::pair<std::vector<int>,std::vector<int>>
Tree::find_leaves() const
{
    return find_support(find_root());
}

int
Tree::find_epl() const
{
    auto [leaves, depths] = find_leaves();
    return std::accumulate(depths.begin(), depths.end(), 0.0);
}

void
Tree::print() const
{
    int root = find_root();
    print_node(root, "", true);
}

void
Tree::print_node(int node, const std::string& prefix, bool is_last) const
{
    std::cout << prefix;

    if (is_last) {
        std::cout << "└── X" << std::endl;
    } else {
        std::cout << "├── X" << std::endl;
    }

    std::vector<int> children = find_children(node);
    for (size_t i = 0; i < children.size(); ++i) {
        print_node(children[i], prefix + (is_last ? "    " : "│   "), i == children.size() - 1);
    }
}