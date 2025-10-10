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
#include <algorithm>    // for std::max and std::binary_search
#include <cmath>        // for sqrt
#include <numeric>      // for std::accumulate and std::transform_reduce
#include <queue>        // for std::queue
#include <random>       // for RNG
#include <stack>        // for std::stack
#include <stdexcept>    // for std::runtime_error, std::invalid_argument
#include <string>       // for std::string
#include <string_view>  // for std::string_view
#include <vector>       // for std::vector

// project-specific includes
#include "tree.hpp"
#include "utils.hpp"

// pybind11 includes
#include <pybind11/pybind11.h>

namespace py = pybind11;


void
Tree::reset()
{
    for(int i = 0; i < n_nodes; ++i)
    {
        set_edge_length(i, 0.0);
        set_subtree_size(i, 1);
        set_parent(i, -1);
        set_child(i, -1);
        set_sibling(i, -1);
        set_is_first(i, true);
        set_name(i, "");
    }
}

// NOTE: These getters can live here since they will never be called from 
//       C++ code, so will never have the chance to be inlined.

py::array_t<int>
Tree::get_subtree_size() const
{
    return std_vec2py_array_t(subtree_size);
}

py::array_t<double>
Tree::get_edge_length() const 
{
    return std_vec2py_array_t(edge_length);
}

void 
Tree::link_(int u, int v)
{
    if(v == -1) 
    {
        return;
    }

    // update pointers
    set_parent(u, v);
    set_sibling(u, get_child(v));
    set_child(v, u);

    // update is_first
    set_is_first(u, true);
    set_is_first(get_sibling(u), false);
}

void
Tree::link(int u, int v)
{
    bool is_leaf = (get_child(v) == -1);

    link_(u, v);

    // update subtree_size
    int diff = get_subtree_size(u) - (is_leaf ? 1 : 0);
    for(; v >= 0; v = get_parent(v))
    {
        set_subtree_size(v, get_subtree_size(v) + diff);
    }
}

void 
Tree::cut_(int u)
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

    set_parent(u, -1);
    set_sibling(u, -1);
}

void 
Tree::cut(int u)
{
    int p = get_parent(u);
    if(p == -1)
    {
        return;
    }
    
    bool was_internal = (get_child(p) != -1);

    cut_(u);

    bool is_leaf = (get_child(p) == -1);

    // update subtree_size
    int diff = (is_leaf && was_internal ? 1 : 0) - get_subtree_size(u);
    for(; p != -1; p = get_parent(p))
    {
        set_subtree_size(p, get_subtree_size(p) + diff);
    }
}

void
Tree::swap_(int u, int v)
{
    int pu = get_parent(u);
    int pv = get_parent(v);
    cut_(u);
    cut_(v);
    link_(u, pv);
    link_(v, pu);
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

void 
Tree::splice_(int u, int v)
{
    cut_(u);
    link_(u, v);
}

void 
Tree::splice(int u, int v)
{
    cut(u);
    link(u, v);
}

std::vector<int>
Tree::find_postorder_(int u) const
{
    std::vector<int> postorder;
    postorder.reserve(n_nodes);
    std::function<void(int)> depth_first_search = [&](int v)
    {
        for(int c = get_child(v); c != -1; c = get_sibling(c))
        {
            depth_first_search(c);
        }
        postorder.emplace_back(v);
    };
    depth_first_search(u);
    return postorder;
}

py::array_t<int>
Tree::find_postorder(int u) const
{
    return std_vec2py_array_t(find_postorder_(u));
}

std::vector<int>
Tree::find_mirror_postorder_(int u) const
{
    std::vector<int> mirror_postorder;
    mirror_postorder.reserve(n_nodes);
    std::function<void(int)> mirror_depth_first_search = [&](int v)
    {
        std::vector<int> children = find_children_(v);
        for(std::size_t i = children.size(); i-- > 0;)
        {
            mirror_depth_first_search(children[i]);
        }
        mirror_postorder.emplace_back(v);
    };
    mirror_depth_first_search(u);
    return mirror_postorder;
}

py::array_t<int>
Tree::find_mirror_postorder(int u) const
{
    return std_vec2py_array_t(find_mirror_postorder_(u));
}

Tree*
Tree::prune(std::vector<int> keep) const
{
    Tree tree_copy(*this);  // heap members are handled properly

    int root = tree_copy.find_root();
    auto [leaves, leaf_depths] = tree_copy.find_leaves_(root);
    auto [intr_nodes, intr_node_depths] = tree_copy.find_interior_nodes_(root);

    for(int leaf : leaves)
    {
        if(!std::binary_search(keep.begin(), keep.end(), leaf))
        {
            tree_copy.cut_(leaf);
        }
    }

    for(int node : intr_nodes)
    {
        if(node == root)
        {
            // do not prune the root
            continue;
        }

        int child = tree_copy.get_child(node);
        if(child == -1)
        {
            tree_copy.cut_(node);
            continue;
        }

        int sibling = tree_copy.get_sibling(child);
        if(sibling == -1)
        {
            int parent = tree_copy.get_parent(node);
            tree_copy.cut_(node);
            tree_copy.cut_(child);
            double edge_length = tree_copy.get_edge_length(child) 
                                    + tree_copy.get_edge_length(node);
            tree_copy.set_edge_length(child, edge_length);
            tree_copy.link_(child, parent);
        }
    }

    std::vector<int> postorder = tree_copy.find_postorder_(root);
    std::unordered_map<int, int> inverse_postorder;
    for(std::size_t i = 0; i < postorder.size(); ++i) 
    {
        inverse_postorder[postorder[i]] = static_cast<int>(i);
    }

    int pruned_n_nodes = static_cast<int>(postorder.size());
    Tree* pruned_tree = new Tree(pruned_n_nodes);
    for(int old_idx : tree_copy.find_mirror_postorder_(root))
    {
        int new_idx = inverse_postorder[old_idx];
        int old_idx_parent = tree_copy.get_parent(old_idx);
        int new_idx_parent = old_idx_parent == -1 ? -1 : inverse_postorder[old_idx_parent];
        pruned_tree->set_edge_length(new_idx, tree_copy.get_edge_length(old_idx));
        pruned_tree->link_(new_idx, new_idx_parent);
    }

    pruned_tree->update_subtree_size();

    return pruned_tree;
}

void
Tree::update_subtree_size()
{
    std::vector<int> postorder = find_postorder_(find_root());
    for(int i = 0; i < n_nodes; ++i)
    {
        int node = postorder[i];
        if(get_child(node) == -1)
        {
            set_subtree_size(node, 1);
            continue;
        }
        std::vector<int> children = find_children_(node);
        int size = std::transform_reduce(
            children.begin(), children.end(),
            0, std::plus<>(),
            [&](int child) { return get_subtree_size(child); }
        );
        set_subtree_size(node, size);
    }
}

std::vector<int>
Tree::find_children_(int u) const
{
    std::vector<int> children;
    for(int c = get_child(u); c != -1; c = get_sibling(c))
    {
        children.push_back(c);
    }
    return children;
}

py::array_t<int> 
Tree::find_children(int u) const
{
    return std_vec2py_array_t(find_children_(u));
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
    return std_vec2py_array_t(find_ancestors_(u));
}

std::pair<std::vector<int>,std::vector<int>>
Tree::find_path_(int u, int v) const 
{
    std::vector<int> u_path = find_ancestors_(u);
    std::vector<int> v_path = find_ancestors_(v);
    while(!u_path.empty() && !v_path.empty() && u_path.back() == v_path.back()) 
    {
        u_path.pop_back();
        v_path.pop_back();
    }
    return std::make_pair(u_path, v_path);
}

py::tuple
Tree::find_path(int u, int v) const
{
    auto [u_path, v_path] = find_path_(u, v);
    return py::make_tuple(
        std_vec2py_array_t(u_path), 
        std_vec2py_array_t(v_path)
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
    return (find_children_(find_root()).size() == 1);
}

std::pair<std::vector<int>, std::vector<int>>
Tree::find_leaves_(int u) const
{
    std::vector<int> leaves;
    std::vector<int> depths;

    leaves.reserve(get_subtree_size(u));
    depths.reserve(get_subtree_size(u));

    std::function<void(int, int)> depth_first_search = [&](int v, int depth)
    {
        for(int c = get_child(v); c != -1; c = get_sibling(c))
        {
            depth_first_search(c, depth + 1);
        }
        if(get_child(v) == -1)
        {
            leaves.emplace_back(v);
            depths.emplace_back(depth);
        }
    };

    depth_first_search(u, 0);

    return std::make_pair(leaves, depths);
}

py::tuple
Tree::find_leaves(int u) const
{
    auto [leaves, depths] = find_leaves_(u);
    return py::make_tuple(
        std_vec2py_array_t(leaves), 
        std_vec2py_array_t(depths)
    );
}

std::pair<std::vector<int>, std::vector<int>>
Tree::find_interior_nodes_(int u) const 
{
    std::vector<int> intr_nodes;
    std::vector<int> depths;

    intr_nodes.reserve(get_subtree_size(u));
    depths.reserve(get_subtree_size(u));

    std::function<void(int, int)> depth_first_search = [&](int v, int depth)
    {
        for(int c = get_child(v); c != -1; c = get_sibling(c))
        {
            depth_first_search(c, depth + 1);
        }
        if(get_child(v) != -1)
        {
            intr_nodes.emplace_back(v);
            depths.emplace_back(depth);
        }
    };

    depth_first_search(u, 0);

    return std::make_pair(intr_nodes, depths);
}

py::tuple
Tree::find_interior_nodes(int u) const
{
    auto [interior_nodes, depths] = find_interior_nodes_(u);
    return py::make_tuple(
        std_vec2py_array_t(interior_nodes), 
        std_vec2py_array_t(depths)
    );
}

std::vector<int>
Tree::find_subtree_start_indices_() const
{
    std::vector<int> subtree_starts(get_n_nodes());

    // breadth-first search
    std::queue<std::pair<int, int>> q;
    q.push({find_root(), 0});
    while(!q.empty())
    {
        auto [u, start] = q.front();
        q.pop();

        subtree_starts[u] = start;

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

py::array_t<int>
Tree::find_subtree_start_indices() const
{
    std::vector<int> subtree_starts = find_subtree_start_indices_();
    return std_vec2py_array_t(subtree_starts);
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
Tree::find_n_children(int u) const
{
    int counter = 0;
    for(int c = nodes[u].child; c != -1; c = nodes[c].sibling)
    {
        counter += 1;
    }
    return counter;
}

int 
Tree::find_n_wavelets() const
{
    int counter = 0;
    for(int i = 0; i < get_n_nodes(); ++i)
    {
        int n_children = find_n_children(i);
        counter += std::max(1, n_children - 1) * (n_children != 0);
    }
    return counter;
}

int
Tree::find_epl() const
{
    auto [leaves, depths] = find_leaves_(find_root());
    return std::accumulate(depths.begin(), depths.end(), 0);
}

double
Tree::find_tbl(int u, int v) const
{
    auto [u_path, v_path] = find_path_(u, v);
    
    double u_sum = find_tbl(u);
    double v_sum = find_tbl(v);
    for(int i = 0; i < static_cast<int>(u_path.size()); ++i)
    {
        u_sum += find_tbl(u_path[i]);
    }
    for(int i = 0; i < static_cast<int>(v_path.size()); ++i)
    {
        v_sum += find_tbl(v_path[i]);
    }

    return u_sum + v_sum;
}

std::vector<double>
Tree::find_tbl_() const 
{
    std::vector<double> tbl(get_n_nodes());
    for(int i = 0; i < get_n_nodes(); ++i)
    {
        tbl[i] = find_tbl(i);
    }
    return tbl;
}

py::array_t<double>
Tree::find_tbl() const
{
    return std_vec2py_array_t(find_tbl_());
}

py::list
Tree::compute_wavelets(int u) const
{
    std::vector<int> children = find_children_(u);

    py::list wavelets;

    int l0 = get_subtree_size(children[0]);

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
    for(py::ssize_t i = 1; i < static_cast<py::ssize_t>(children.size()); ++i)
    {
        int l1 = get_subtree_size(children[i]);
        int s = l0 + l1;

        double l0_ = static_cast<double>(l0);
        double l1_ = static_cast<double>(l1);
        double s_  = static_cast<double>(s);

        double val0 =  sqrt(l1_ / (l0_ * s_));
        double val1 = -sqrt(l0_ / (l1_ * s_));

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
    auto [leaves, depths] = find_leaves_(u);

    std::vector<int> children = find_children_(u);

    py::list supports;

    int size = get_subtree_size(children[0]);

    // mother wavelet
    if(children.size() == 1)
    {
        py::array_t<int> support(size);
        auto support_ = support.mutable_unchecked<1>();
        for(py::ssize_t j = 0; j < static_cast<py::ssize_t>(size); ++j)
        {
            support_(j) = leaves[j];
        }

        supports.append(support);
        return supports;
    }

    // daughter wavelet
    for(py::ssize_t i = 1; i < static_cast<py::ssize_t>(children.size()); ++i)
    {
        size += get_subtree_size(children[i]);

        py::array_t<int> support(size);
        auto support_ = support.mutable_unchecked<1>();
        for(py::ssize_t j = 0; j < static_cast<py::ssize_t>(size); ++j)
        {
            support_(j) = leaves[j];
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

    std::vector<int> children = find_children_(node);
    for (int i = 0; i < static_cast<int>(children.size()); ++i)
    {
        result += to_string_(children[i], prefix + (is_last ? "    " : "│   "),
                             i == static_cast<int>(children.size()) - 1, label);
    }
    return result;
}

// first pass: count the number of nodes in the tree
std::pair<int, bool>
count_nodes(std::string_view nwk_str)
{
    int count = 1;
    int depth = 0;
    int top_children = 1;
    for(char c: nwk_str)
    {
        switch(c)
        {
            case '(':
            {
                count++;
                depth++;
                break;
            }
            case ')':
            {
                depth--;
                break;
            }
            case ',':
            {
                count++;
                top_children += (depth == 0);
                break;
            }
            default:
            {
                break;
            }
        }
    }
    bool is_planted = (top_children > 1);
    return std::make_pair(count, is_planted);
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
        {
            flush_buffer(curr_node, tree, char_buf, stack, is_name); 
            break;
        }
        case '(':   // subtree begin
        {
            stack.push(-1); 
            break;
        }
        case ',':   // sibling 
        {
            flush_buffer(curr_node, tree, char_buf, stack, is_name);
            stack.push(curr_node++);
            is_name = true;
            break;
        }
        case ')':   // subtree end
        {
            flush_buffer(curr_node, tree, char_buf, stack, is_name);
            stack.push(curr_node++);
            while(!stack.empty() && stack.top() != -1) 
            {
                int child = stack.top();
                tree->link_(child, curr_node);
                stack.pop();
            }
            if(!stack.empty()) stack.pop(); // pop '-1'
            is_name = true;
            break;
        }
        case ':':   // edge length
        {
            flush_buffer(curr_node, tree, char_buf, stack, is_name);
            is_name = false;
            break;
        }
        default:
        {
            if(is_valid_char[(unsigned char)c]) char_buf.push_back(c);
        }
    }
}

// second pass: assign parent-child relationships
Tree*
nwk2tree(const std::string& nwk_str, bool ensure_planted)
{
    // allocate a new Tree
    auto [n_nodes, is_planted] = count_nodes(nwk_str);
    bool do_plant = (ensure_planted && !is_planted);
    Tree *tree = new Tree(n_nodes + do_plant);

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
    bool is_quote_char[256] = {false};
    is_quote_char[(unsigned char)'\''] = true;
    is_quote_char[(unsigned char)'\"'] = true;
    char curr_quote_char;

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

        if((!is_quoted && is_quote_char[(unsigned char)c]) || (is_quoted && c == curr_quote_char)) 
        {
            curr_quote_char = c; // set quote character
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

    // plant tree
    if(do_plant)
    {
        tree->link_(tree->get_n_nodes() - 2, tree->get_n_nodes() - 1);
    }

    tree->update_subtree_size();

    if(parens > 0) 
    {
        throw std::runtime_error("Unbalanced parentheses: extra '('");
    }

    return tree;
}