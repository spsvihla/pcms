/**
 * @file tree_node.cpp
 * 
 * Contains the implementation of the TreeNode class and related functions.
 * This file handles tree node operations such as swapping subtrees and computing 
 * Lowest Common Ancestors (LCAs) in a tree structure.
 * 
 * @author Sean Svihla
 */

// standard library includes
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <stack>
#include <string>
#include <string_view>
#include <vector>

// project-specific includes
#include "tree.hpp"

// pybind11 includes
#include <pybind11/pybind11.h>

namespace py = pybind11;


// constructor
Tree::Tree(int n_nodes)
: n_nodes(n_nodes), topology(n_nodes), numerics(n_nodes), names(n_nodes) {}

// destructor
Tree::~Tree() {}

int 
Tree::get_size() const
{
    return n_nodes;
}

int 
Tree::get_parent(int u) const
{
    return topology.at(u).parent;
}

int 
Tree::get_child(int u) const
{
    return topology.at(u).child;
}

int 
Tree::get_sibling(int u) const
{
    return topology.at(u).sibling;
}

int 
Tree::get_subtree_size(int u) const
{
    return numerics.subtree_size.at(u);
}

double 
Tree::get_edge_length(int u) const
{
    return numerics.edge_length.at(u);
}

void 
Tree::set_edge_length(int u, double value)
{
    numerics.edge_length.at(u) = value;
}

std::string
Tree::get_name(int u) const
{
    return names.at(u);
}

void
Tree::set_name(int u, const std::string& value)
{
    names.at(u) = value;
}

void 
Tree::link(int u, int v)
{
    // Check bounds
    if(u < 0 || u >= n_nodes)
    {
        throw std::out_of_range("Node indices out of bounds: u = " + std::to_string(u));
    }
    if(v < -1 || v >= n_nodes)
    {
        throw std::out_of_range("Node indices out of bounds: v = " + std::to_string(v));
    }

    // Check if u is already linked
    if(topology[u].parent != -1) 
    {
        throw std::runtime_error("Cannot link a node before cutting");
    }

    if(v == -1) return;

    int is_leaf = (topology[v].child == -1);

    // update pointers
    topology[u].parent = v;
    topology[u].sibling = topology[v].child;
    topology[v].child = u;

    // update subtree_size
    int diff = numerics.subtree_size[u] - is_leaf;
    while(v >= 0)
    {
        numerics.subtree_size[v] += diff;
        v = topology[v].parent;
    }
}

void 
Tree::cut(int u)
{
    // Check bounds
    if(u < 0 || u >= n_nodes)
    {
        throw std::out_of_range("Node indices out of bounds");
    }

    int p = topology[u].parent;
    if(p == -1)
    {
        // node is an orphan
        return;
    }

    // update pointers
    if(topology[p].child == u)
    {
        topology[p].child = topology[u].sibling;
    }
    else
    {
        int c = topology[p].child;
        while(topology[c].sibling != u)
        {
            c = topology[c].sibling;
        }
        topology[c].sibling = topology[u].sibling;
    }
    int is_leaf = (topology[p].child == -1);

    // update subtree_size
    int diff = is_leaf - numerics.subtree_size[u];
    for(; p != -1; p = topology[p].parent)
    {
        numerics.subtree_size[p] += diff;
    }

    topology[u].parent = -1;
    topology[u].sibling = -1;
}

void 
Tree::swap(int u, int v)
{
    int pu = topology.at(u).parent;
    int pv = topology.at(v).parent;
    cut(u);
    cut(v);
    link(u, pv);
    link(v, pu);
}

std::vector<int> 
Tree::find_children(int u) const
{
    // Check bounds
    if(u < 0 || u >= n_nodes)
    {
        throw std::out_of_range("Node indices out of bounds");
    }

    std::vector<int> children;
    for(int c = topology[u].child; c != -1; c = topology[c].sibling)
    {
        children.push_back(c);
    }
    return children;
}

std::vector<int> 
Tree::find_ancestors(int u) const
{
    // Check bounds
    if(u < 0 || u >= n_nodes)
    {
        throw std::out_of_range("Node indices out of bounds");
    }

    std::vector<int> ancestors;
    for(u = topology[u].parent; u != -1; u = topology[u].parent)
    {
        ancestors.push_back(u);
    }
    return ancestors;
}

std::pair<std::vector<int>, std::vector<int>> 
Tree::find_support(int u) const
{
    // Check bounds
    if(u < 0 || u >= n_nodes)
    {
        throw std::out_of_range("Node indices out of bounds");
    }

    std::vector<int> support;
    std::vector<int> depths;
    std::stack<std::pair<int, int>> stack;
    stack.push(std::make_pair(u, 0));

    while(!stack.empty())
    {
        int v = stack.top().first;
        int depth = stack.top().second;
        stack.pop();

        if(topology[v].child == -1)
        {
            support.push_back(v);     // leaf
            depths.push_back(depth);  // depth
        }

        for(int c = topology[v].child; c != -1; c = topology[c].sibling)
        {
            stack.push(std::make_pair(c, depth + 1));
        }
    }
    return std::make_pair(support, depths);
}

std::pair<std::vector<int>,std::vector<int>> 
Tree::find_path(int u, int v) const
{
    // Check bounds
    if(u < 0 || u >= n_nodes)
    {
        throw std::out_of_range("Node indices out of bounds: u = " + std::to_string(u));
    }
    if(v < 0 || v >= n_nodes)
    {
        throw std::out_of_range("Node indices out of bounds: v = " + std::to_string(v));
    }

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
    std::vector<int> ancestors = find_ancestors(0);
    if(!ancestors.empty())
    {
        return ancestors.back();
    }
    else
    {
        // no ancestors, so node 0 is root
        return 0;
    }
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
Tree::print(const std::string& label) const
{
    int root = find_root();
    print_node(root, "", true, label);
}

void
Tree::print_node(int node, const std::string& prefix, bool is_last, 
                 const std::string& label) const
{
    std::cout << prefix;

    // add branches
    if (is_last) 
    {
        std::cout << "└── ";
    } 
    else 
    {
        std::cout << "├── ";
    }

    // add labels or X marker
    if(label == "index") 
    {
        std::cout << node << std::endl;
    } 
    else if(label == "name")
    {
        std::cout << names[node] << std::endl;
    }
    else if(label == "none")
    {
        std::cout << "X" << std::endl;
    }
    else
    {
        throw py::value_error("Unrecognized labelling: " + label);
    }

    std::vector<int> children = find_children(node);
    for(size_t i = 0; i < children.size(); ++i) 
    {
        print_node(children[i], prefix + (is_last ? "    " : "│   "), 
                   i == children.size() - 1, label);
    }
}

std::string
read_nwk(std::string filename)
{
    std::ifstream file(filename);
    if(!file)
    {
        throw std::runtime_error("Could not open file " + filename);
    }

    std::string nwk_str;
    std::getline(file, nwk_str); // read the entire Newick string from the file

    file.close();
    
    return nwk_str;
}

/* first pass: count the number of nodes in the tree */
int 
count_nodes(std::string_view nwk_str)
{
    int count = 1;
    for(char c: nwk_str)
    {
        if(c == '(' || c == ',')
        {
            count++;
        }
    }
    return count;
}

void
flush_buffer(int& curr_node, Tree* tree, std::vector<char>& char_buf, 
             std::stack<int>& stack, bool is_name)
{
    if (char_buf.empty()) return;

    std::string str(char_buf.begin(), char_buf.end());  
    char_buf.clear();

    if(is_name)
    {
        tree->set_name(curr_node, str);
    }
    else
    {
        try
        {
            tree->set_edge_length(curr_node, std::stod(str));
        }
        catch(const std::invalid_argument& e)
        {
            throw py::value_error("Invalid edge length: " + str);
        }
    }
}

inline void
parse_newick(char c, int& curr_node, Tree* tree, std::vector<char>& char_buf, 
             std::stack<int>& stack, bool& is_name, bool is_valid_char[256])
{
    switch(c) 
    {
        case ';':   // root
            flush_buffer(curr_node, tree, char_buf, stack, is_name); 
            break;
        case '(':   // subtree begin
            stack.push(-1); 
            break;
        case ',':   // sibling 
            flush_buffer(curr_node, tree, char_buf, stack, is_name);
            stack.push(curr_node++);
            is_name = true;
            break;
        case ')':   // subtree end
        {
            flush_buffer(curr_node, tree, char_buf, stack, is_name);
            stack.push(curr_node++);
            while(!stack.empty() && stack.top() != -1) 
            {
                int child = stack.top();
                tree->link(child, curr_node);
                stack.pop();
            }
            if(!stack.empty()) stack.pop(); // pop '-1'
            is_name = true;
            break;
        }
        case ':':   // edge length
            flush_buffer(curr_node, tree, char_buf, stack, is_name);
            is_name = false;
            break;
        default:
            if(is_valid_char[(unsigned char)c]) char_buf.push_back(c);
    }
}

/* second pass: assign parent-child relationships */
Tree
*nwk2tree(const std::string& filename)
{
    std::string nwk_str = read_nwk(filename);

    // allocate a new Tree
    int n_nodes = count_nodes(nwk_str);
    Tree *tree = new Tree(n_nodes);

    // create a stack for storing subtree parents
    std::stack<int> stack;
    int curr_node = 0;

    // create a lookup table for alpha-numeric characters
    bool is_valid_char[256] = {false};

    // initialize valid characters
    for(char c = '0'; c <= '9'; c++) is_valid_char[(unsigned char)c] = true;
    for(char c = 'A'; c <= 'Z'; c++) is_valid_char[(unsigned char)c] = true;
    for(char c = 'a'; c <= 'z'; c++) is_valid_char[(unsigned char)c] = true;
    is_valid_char[(unsigned char)'_'] = true;
    is_valid_char[(unsigned char)'.'] = true;

    // create a character buffer for temporary name/edge length storage
    bool is_name = true;
    std::vector<char> char_buf;

    // all characters allowed if quoted
    bool is_quoted = false;

    for(int i = 0; i < static_cast<int>(nwk_str.size()); i++)
    {
        char c = nwk_str[i];

        if(c == '\'' || c == '\"') 
        {
            is_quoted ^= 1; // toggle quoting
        } 
        else if(is_quoted) 
        {
            char_buf.push_back(c);
        } 
        else 
        {
            parse_newick(c, curr_node, tree, char_buf, stack, is_name, 
                         is_valid_char);
        }
    }

    return tree;
}