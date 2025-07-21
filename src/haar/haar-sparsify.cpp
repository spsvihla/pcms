/**
 * @file haar-sparsify.cpp
 * 
 * This file implements the sparsify function and helper functions.
 * 
 * @author Sean Svila
 */

// standard library includes
#include <algorithm>
#include <numeric>
#include <stack>
#include <vector>

// project-specific includes
#include "tree.hpp"
#include "haar-sparsify.hpp"

// pybind11 and numpy includes
#include <pybind11/pybind11.h>

namespace py = pybind11;


// compute number of non-zero entries in wavelet basis matrix
int
compute_nnz_wavelets(Tree* tree)
{
    int nnz = 0;
    for(int i = 0; i < tree->get_n_nodes(); ++i)
    {
        if(tree->get_child(i) == -1L)
        {
            continue;
        }
        py::array_t<int> children = tree->find_children(i);
        auto children_ = children.unchecked<1>();
        int n_children = static_cast<int>(children.size());
        nnz += tree->get_subtree_size(children_[0]) * std::max(1, n_children - 1);
        for(int j = 1; j < n_children; ++j)
        {
            nnz += tree->get_subtree_size(children_[j]) * (n_children - j);
        }
    }
    return nnz;
}

py::tuple
sparsify(Tree* tree)
{
    int nnz = compute_nnz_wavelets(tree);

    int n_leaves = tree->find_n_leaves();

    py::array_t<int> subtree_start = tree->find_subtree_start_indices();
    auto subtree_start_ = subtree_start.unchecked<1>();

    py::array_t<double> Q_values(static_cast<py::ssize_t>(nnz));
    py::array_t<double> S_values(static_cast<py::ssize_t>(nnz));
    py::array_t<double> trace_length(static_cast<py::ssize_t>(n_leaves));
    py::array_t<py::ssize_t> indices(static_cast<py::ssize_t>(nnz));
    py::array_t<py::ssize_t> indptr(static_cast<py::ssize_t>(n_leaves+1));

    double* Q_values_ = static_cast<double*>(Q_values.request().ptr);
    double* S_values_ = static_cast<double*>(S_values.request().ptr);
    double* trace_length_ = static_cast<double*>(trace_length.request().ptr);
    py::ssize_t* indices_ = static_cast<py::ssize_t*>(indices.request().ptr);
    py::ssize_t* indptr_ = static_cast<py::ssize_t*>(indptr.request().ptr);

    indptr_[0] = 0;

    for(int i = 0; i < n_leaves; ++i)
    {
        trace_length_[i] = 0.0;
    }

    // build wavelet basis and trace length sparse matrices
    std::size_t values_idx = 0;
    std::size_t indptr_idx = 0;

    std::stack<int> subtree_start_stack;
    std::stack<int> subtree_size_stack;

    for(int i = 0; i < tree->get_n_nodes() - 2; ++i)
    {
        // accumulate trace length
        int start = subtree_start_[static_cast<py::ssize_t>(i)];
        int size = tree->get_subtree_size(i);
        double tbl = tree->find_tbl(i);
        for(int j = 0; j < size; ++j)
        {
            trace_length_[static_cast<std::size_t>(start + j)] += tbl;
        }

        // construct wavelets
        if(tree->get_is_first(i))           // first child node
        {
            subtree_start_stack.push(start);
            subtree_size_stack.push(size);
        }
        else                                // daughter wavelet
        {
            start = subtree_start_stack.top();

            int rsize = size;
            int lsize = subtree_size_stack.top();
            int sum = lsize + rsize;

            // use floating-point arithmetic to avoid overflow with large integers
            double lval =  sqrt(static_cast<double>(rsize) / (static_cast<double>(lsize) * static_cast<double>(sum)));
            double rval = -sqrt(static_cast<double>(lsize) / (static_cast<double>(rsize) * static_cast<double>(sum)));

            // left subtree
            for(std::size_t j = 0; j < static_cast<std::size_t>(lsize); ++j)
            {
                S_values_[values_idx] = trace_length_[start + j] * lval;
                Q_values_[values_idx] = lval;
                indices_[values_idx] = static_cast<py::ssize_t>(start + j);
                values_idx++;
            }

            // right subtree
            for(std::size_t j = lsize; j < static_cast<std::size_t>(sum); ++j)
            {
                S_values_[values_idx] = trace_length_[start + j] * rval;
                Q_values_[values_idx] = rval;
                indices_[values_idx] = static_cast<py::ssize_t>(start + j);
                values_idx++;
            }

            indptr_[++indptr_idx] = values_idx;
            subtree_size_stack.top() = sum;
        }
        if(tree->get_sibling(i) == -1L)      // last child node
        {
            subtree_start_stack.pop();
            subtree_size_stack.pop();
        }
    }

    // mother wavelet
    int root = tree->find_root();
    int size = tree->get_subtree_size(root);

    double tbl = tree->find_tbl(root - 1);
    for(int j = 0; j < size; ++j)
    {
        trace_length_[j] += tbl;
    }

    double val = sqrt(1.0 / size);
    for(std::size_t j = 0; j < static_cast<std::size_t>(size); ++j)
    {
        S_values_[values_idx] = trace_length_[j] * val;
        Q_values_[values_idx] = val;
        indices_[values_idx] = static_cast<py::ssize_t>(j);
        values_idx++;
    }

    indptr_[++indptr_idx] = values_idx;

    py::tuple Q = py::make_tuple(Q_values, indices, indptr);
    py::tuple S = py::make_tuple(S_values, indices.attr("copy")(), indptr.attr("copy")());

    return py::make_tuple(Q, S);
}