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
#include <queue>
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


inline py::array_t<int>
std_vec2py_array_t_int(std::vector<int> arr)
{
    py::array_t<int> result(arr.size());
    auto result_ = result.mutable_unchecked<1>();
    for (ssize_t i = 0; i < static_cast<ssize_t>(arr.size()); ++i)
        result_(i) = arr[i];
    return result;
}


// constructor
Tree::Tree(int n_nodes)
: n_nodes(n_nodes), topology(n_nodes), numerics(n_nodes), names(n_nodes) {}

// destructor
Tree::~Tree() {}

void 
Tree::link(int u, int v)
{
    if(v == -1) 
    {
        return;
    }

    int is_leaf = (get_child(v) == -1);

    // update pointers
    set_parent(u, v);
    set_sibling(u, get_child(v));
    set_child(v, u);

    // update is_first
    set_is_first(u, true);
    set_is_first(get_sibling(u), false);

    // update subtree_size
    int diff = get_subtree_size(u) - is_leaf;
    for(; v >= 0; v = get_parent(v))
    {
        set_subtree_size(v, get_subtree_size(v) + diff);
    }
}

void 
Tree::cut(int u)
{
    int p = get_parent(u);
    
    if(p == -1)
    {
        return;
    }

    // update pointers
    if(get_child(p) == u)
    {
        set_child(p, get_sibling(u));

        // update is_first
        if(int s = get_sibling(u); s != -1)
        {
            set_is_first(u, false);
            set_is_first(s, true);
        }
    }
    else
    {
        int c = get_child(p);
        while(get_sibling(c) != u)
        {
            c = get_sibling(c);
        }
        set_sibling(c, get_sibling(u));
    }
    int is_leaf = (get_child(p) == -1);

    // update subtree_size
    int diff = is_leaf - get_subtree_size(u);
    for(; p != -1; p = get_parent(p))
    {
        set_subtree_size(p, get_subtree_size(p) + diff);
    }

    set_parent(u, -1);
    set_sibling(u, -1);
}

void 
Tree::swap(int u, int v)
{
    int pu = get_parent(u);
    int pv = get_parent(v);
    cut(u);
    cut(v);
    link(u, pv);
    link(v, pu);
}

py::array_t<int> 
Tree::find_children(int u) const
{
    std::vector<int> children;
    for(int c = get_child(u); c != -1; c = get_sibling(c))
    {
        children.push_back(c);
    }
    return std_vec2py_array_t_int(children);
}

std::vector<int>
Tree::find_ancestors_(int u) const
{
    std::vector<int> ancestors;
    for(int p = get_parent(u); p != -1; p = get_parent(p))
    {
        ancestors.push_back(p);
    }
    return ancestors;
}

py::array_t<int> 
Tree::find_ancestors(int u) const
{
    std::vector<int> ancestors = find_ancestors_(u);
    if(ancestors.empty())
    {
        return py::array_t<int>();
    }
    return std_vec2py_array_t_int(ancestors);
}

std::pair<py::array_t<int>,py::array_t<int>> 
Tree::find_path(int u, int v) const
{
    std::vector<int> u_path = find_ancestors_(u);
    std::vector<int> v_path = find_ancestors_(v);
    while(!u_path.empty() && !v_path.empty() && u_path.back() == v_path.back()) 
    {
        u_path.pop_back();
        v_path.pop_back();
    }
    return std::make_pair(
        std_vec2py_array_t_int(u_path), 
        std_vec2py_array_t_int(v_path)
    );
}

int
Tree::find_root() const
{
    std::vector<int> ancestors = find_ancestors_(0);
    if(!ancestors.empty())
    {
        return ancestors.back();
    }
    else
    {
        return 0;
    }
}

bool
Tree::find_is_planted() const
{
    return (find_children(find_root()).size() == 1);
}

std::pair<py::array_t<int>, py::array_t<int>> 
Tree::find_leaves(int u) const
{
    py::array_t<int> support(get_subtree_size(u));
    py::array_t<int> depths(get_subtree_size(u));

    auto support_ = support.mutable_unchecked<1>();
    auto depths_ = depths.mutable_unchecked<1>();
    
    // depth-first search
    int idx = 0;
    std::stack<std::pair<int, int>> stack;
    stack.push(std::make_pair(u, 0));
    while(!stack.empty())
    {
        int v = stack.top().first;
        int depth = stack.top().second;
        stack.pop();

        if(topology[v].child == -1)
        {
            support_(idx) = v;
            depths_(idx++) = depth;
        }

        for(int c = topology[v].child; c != -1; c = topology[c].sibling)
        {
            stack.push(std::make_pair(c, depth + 1));
        }
    }

    return std::make_pair(support, depths);
}

py::array_t<int>
Tree::find_subtree_start_indices() const
{
    py::array_t<int> subtree_starts(get_n_nodes());
    auto subtree_starts_ = subtree_starts.mutable_unchecked<1>();

    // breadth-first search
    std::queue<std::pair<int, int>> q;
    q.push({find_root(), 0});
    while(!q.empty())
    {
        auto [u, start] = q.front();
        q.pop();

        subtree_starts_(u) = start;

        // enqueue children
        int offset = 0;
        for(int c = get_child(u); c != -1; c = get_sibling(c))
        {
            q.push({c, start + offset});
            offset += get_subtree_size(c);
        }
    }

    return subtree_starts;
}

int 
Tree::find_n_leaves() const
{
    int counter = 0;
    for(int i = 0; i < get_n_nodes(); ++i)
    {
        counter += static_cast<int>(get_child(i) == -1);
    }
    return counter;
}

int
Tree::find_epl() const
{
    auto [leaves, depths] = find_leaves(find_root());
    auto depths_ = depths.mutable_unchecked<1>();
    int epl = 0;
    for(py::ssize_t i = 0; i < depths.size(); ++i)
    {
        epl += depths_[i];
    }
    return epl;
}

double
Tree::find_tbl(int u, int v) const
{
    auto [u_path, v_path] = find_path(u, v);
    
    auto u_path_ = u_path.unchecked<1>();
    auto v_path_ = v_path.unchecked<1>();
    
    double u_sum = find_tbl(u);
    double v_sum = find_tbl(v);
    for(py::ssize_t i = 0; i < u_path.size(); ++i)
    {
        u_sum += find_tbl(u_path_[i]);
    }
    for(py::ssize_t i = 0; i < v_path.size(); ++i)
    {
        v_sum += find_tbl(v_path_[i]);
    }

    return u_sum + v_sum;
}

py::array_t<double>
Tree::find_tbl() const
{
    py::array_t<double> tbl(n_nodes);
    auto tbl_ = tbl.mutable_unchecked<1>();
    for(py::ssize_t i = 0; i < static_cast<py::ssize_t>(n_nodes); ++i)
    {
        tbl_(i) = find_tbl(static_cast<int>(i));
    }
    return tbl;
}

py::list
Tree::compute_wavelets(int u) const
{
    py::array_t<int> children = find_children(u);
    auto children_ = children.unchecked<1>();

    py::list wavelets;

    int l0 = get_subtree_size(children_[0]);

    // mother wavelet
    if(children.size() == 1)
    {
        double val = 1.0 / sqrt(l0);
        
        py::array_t<double> wavelet(static_cast<py::ssize_t>(l0));
        auto wavelet_ = wavelet.mutable_unchecked<1>();
        for(py::ssize_t j = 0; j < static_cast<py::ssize_t>(l0); ++j)
        {
            wavelet_(j) = val;
        }

        wavelets.append(wavelet);
        return wavelets;
    }

    // daughter wavelet
    for(py::ssize_t i = 1; i < children.size(); ++i)
    {
        int l1 = get_subtree_size(children_[i]);
        int s = l0 + l1;
        double val0 = sqrt(static_cast<double>(l1) / (l0 * s));
        double val1 = -1 * sqrt(static_cast<double>(l0) / (l1 * s));

        py::array_t<double> wavelet(static_cast<py::ssize_t>(s));
        auto wavelet_ = wavelet.mutable_unchecked<1>();
        for(py::ssize_t j = 0; j < static_cast<py::ssize_t>(l0); ++j)
        {
            wavelet_(j) = val0;
        }
        for(py::ssize_t j = l0; j < static_cast<py::ssize_t>(s); ++j)
        {
            wavelet_(j) = val1;
        }

        wavelets.append(wavelet);
        l0 += l1;
    }
    return wavelets;
}

py::list
Tree::compute_supports(int u) const
{
    auto [leaves, depths] = find_leaves(u);
    auto leaves_ = leaves.unchecked<1>();

    py::array_t<int> children = find_children(u);
    auto children_= children.unchecked<1>();

    py::list supports;

    int s = get_subtree_size(children_[0]);

    // mother wavelet
    if(children.size() == 1)
    {
        py::array_t<int> support(s);
        auto support_ = support.mutable_unchecked<1>();
        for(py::ssize_t j = 0; j < static_cast<py::ssize_t>(s); ++j)
        {
            support_(j) = leaves_[j];
        }

        supports.append(support);
        return supports;
    }

    // daughter wavelet
    for(py::ssize_t i = 1; i < children.size(); ++i)
    {
        s += get_subtree_size(children_[i]);

        py::array_t<int> support(s);
        auto support_ = support.mutable_unchecked<1>();
        for(py::ssize_t j = 0; j < static_cast<py::ssize_t>(s); ++j)
        {
            support_(j) = leaves_[j];
        }

        supports.append(support);
    }

    return supports;
}

std::string
Tree::to_string(const std::string& label) const
{
    if (!(label == "index" || label == "name" || label == "none"))
    {
        throw py::value_error("Unrecognized labelling: " + label);
    }
    return to_string_(find_root(), "", true, label);
}

std::string
Tree::to_string_(int node, const std::string& prefix, bool is_last,
                 const std::string& label) const
{
    std::string result = prefix;

    // add branches
    result += (is_last ? "└── " : "├── ");

    // add labels or X marker
    if (label == "index")
    {
        result += std::to_string(node) + "\n";
    }
    else if (label == "name")
    {
        result += names[node] + "\n";
    }
    else if (label == "none")
    {
        result += "X\n";
    }

    py::array_t<int> children = find_children(node);
    auto children_ = children.unchecked<1>();
    for (int i = 0; i < static_cast<int>(children.size()); ++i)
    {
        result += to_string_(children_[i], prefix + (is_last ? "    " : "│   "),
                             i == static_cast<int>(children.size()) - 1, label);
    }
    return result;
}

// first pass: count the number of nodes in the tree
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
            throw std::runtime_error("Invalid edge length: " + str);
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

// second pass: assign parent-child relationships
Tree
*nwk2tree(const std::string& nwk_str)
{
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

    // count parentheses
    int parens = 0;

    for(int i = 0; i < static_cast<int>(nwk_str.size()); ++i)
    {
        char c = nwk_str[i];

        parens += (c == '(');
        parens -= (c == ')');
        if(parens < 0) 
        {
            throw std::runtime_error("Unbalanced parentheses: extra ')'");
        }

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

    if(parens > 0) 
    {
        throw std::runtime_error("Unbalanced parentheses: extra '('");
    }

    return tree;
}