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
py::ssize_t
compute_nnz_wavelets(Tree* tree)
{
    py::ssize_t nnz = 0;
    for(int i = 0; i < tree->get_n_nodes(); ++i)
    {
        if(tree->get_child(i) == -1)
        {
            continue;
        }
        py::array_t<int> children = tree->find_children(i);
        auto children_ = children.unchecked<1>();
        py::ssize_t n_children = children.size();
        nnz += tree->get_subtree_size(children_[0]) * std::max(1L, n_children - 1);
        for(py::ssize_t j = 1; j < n_children; ++j)
        {
            nnz += tree->get_subtree_size(children_[j]) * (n_children - j);
        }
    }
    return nnz;
}

// compute upper bound on number of non-zero entries in sparsified covariance matrix
// TODO: generalize bound to non-binary trees
py::ssize_t
compute_nnz_max_cov(Tree* tree)
{
    int n = tree->find_n_leaves();
    int epl = tree->find_epl();
    double frac = 2.0 * (epl + 1) / (n * n) - 3.0 / n;
    return static_cast<py::ssize_t>(ceil(frac * n * n));
}

// double precision dot product
inline double 
dsdot(
    const py::detail::unchecked_reference<double, 1L>& A_values_,
    const py::detail::unchecked_reference<py::ssize_t, 1L>& A_indices_,
    py::ssize_t a_start, py::ssize_t a_end,
    const py::detail::unchecked_reference<double, 1L>& B_values_,
    const py::detail::unchecked_reference<py::ssize_t, 1L>& B_indices_,
    py::ssize_t b_start, py::ssize_t b_end)
{
    double result = 0.0;
    py::ssize_t i = a_start;
    py::ssize_t j = b_start;
    while(i < a_end && j < b_end)
    {
        py::ssize_t a_row = A_indices_[i];
        py::ssize_t b_row = B_indices_[j];
        if(a_row == b_row)
        {
            result += A_values_[i] * B_values_[j];
            ++i;
            ++j;
        }
        else if(a_row < b_row)
        {
            ++i;
        }
        else
        {
            ++j;
        }
    }
    return result;
}

// double precision sparse matrix-transpose-matrix multiply
// C <- A^T * B
// NOTE: This implementation takes advantage of the unique structure of
//       Haar-like basis matrices; it is not a functional general sparse
//       matrix multiply.
inline py::tuple
dspmtm(const py::detail::unchecked_reference<double, 1L>& A_values_, 
       const py::detail::unchecked_reference<py::ssize_t, 1L>& A_indices_, 
       const py::detail::unchecked_reference<py::ssize_t, 1L>& A_indptr_, 
       const py::detail::unchecked_reference<double, 1L>& B_values_,
       const py::detail::unchecked_reference<py::ssize_t, 1L>& B_indices_,
       const py::detail::unchecked_reference<py::ssize_t, 1L>& B_indptr_,
       py::ssize_t nnz_max)
{
    // allocate memory for new matrix
    py::array_t<double> values(nnz_max);
    py::array_t<py::ssize_t> indices(nnz_max);
    py::array_t<py::ssize_t> indptr(B_indptr_.size());

    auto values_ = values.mutable_unchecked<1>();
    auto indices_ = indices.mutable_unchecked<1>();
    auto indptr_ = indptr.mutable_unchecked<1>();

    indptr_(0) = 0;

    py::ssize_t values_idx = 0;
    py::ssize_t indptr_idx = 0;

    for(py::ssize_t i = 0; i < B_indptr_.size() - 1; ++i)
    {
        py::ssize_t bi0 = B_indptr_[i];             // start index
        py::ssize_t bi1 = B_indptr_[i+1];           // end index
        py::ssize_t br0 = B_indices_[bi0];          // start row
        py::ssize_t br1 = B_indices_[bi1-1];        // end row

        for(py::ssize_t j = 0; j < A_indptr_.size() - 1; ++j)
        {
            py::ssize_t ai0 = A_indptr_[j];         // start index
            py::ssize_t ai1 = A_indptr_[j+1];       // end index
            py::ssize_t ar0 = A_indices_[ai0];      // start row
            py::ssize_t ar1 = A_indices_[ai1-1];    // edn row

            // skip if rows to not overlap
            if(br1 < ar0 || ar1 < br0)
            {
                continue;
            }

            // compute dot product
            double val = dsdot(
                A_values_, A_indices_, ai0, ai1,
                B_values_, B_indices_, bi0, bi1
            );

            if(val != 0) 
            {
                values_(values_idx) = val;
                indices_(values_idx) = j;
                ++values_idx;
            }
        }

        indptr_(++indptr_idx) = values_idx;
    }

    values.resize({values_idx});
    indices.resize({values_idx});

    return py::make_tuple(values, indices, indptr);
}

py::tuple
sparsify(Tree* tree)
{
    py::ssize_t nnz = compute_nnz_wavelets(tree);

    int n_leaves = tree->find_n_leaves();

    py::array_t<int> subtree_starts = tree->find_subtree_start_indices();
    auto subtree_starts_ = subtree_starts.unchecked<1>();

    py::array_t<double> Q_values(nnz);
    py::array_t<double> S_values(nnz);
    py::array_t<double> trace_length(n_leaves);
    py::array_t<py::ssize_t> indices(nnz);
    py::array_t<py::ssize_t> indptr(n_leaves+1);

    auto Q_values_ = Q_values.mutable_unchecked<1>();
    auto S_values_ = S_values.mutable_unchecked<1>();
    auto trace_length_ = trace_length.mutable_unchecked<1>();
    auto indices_ = indices.mutable_unchecked<1>();
    auto indptr_ = indptr.mutable_unchecked<1>();

    indptr_(0) = 0;

    for(py::ssize_t i = 0; i < n_leaves; ++i)
    {
        trace_length_(i) = 0.0;
    }

    // build wavelet basis and trace length sparse matrices
    py::ssize_t values_idx = 0;
    py::ssize_t indptr_idx = 0;

    std::stack<int> subtree_start_stack;
    std::stack<int> subtree_size_stack;
    for(int i = 0; i < tree->get_n_nodes() - 2; ++i)
    {
        // accumulate trace length
        int start = subtree_starts_[static_cast<py::ssize_t>(i)];
        int size = tree->get_subtree_size(i);
        double tbl = tree->find_tbl(i);
        for(int j = 0; j < size; ++j)
        {
            trace_length_(start + j) += tbl;
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
            for(int j = 0; j < lsize; ++j)
            {
                S_values_(values_idx) = trace_length_[start + j];
                Q_values_(values_idx) = lval;
                indices_(values_idx) = static_cast<py::ssize_t>(start + j);
                values_idx++;
            }

            // right subtree
            for(int j = lsize; j < sum; ++j)
            {
                S_values_(values_idx) = trace_length_[start + j];
                Q_values_(values_idx) = rval;
                indices_(values_idx) = static_cast<py::ssize_t>(start + j);
                values_idx++;
            }

            indptr_(++indptr_idx) = values_idx;
            subtree_size_stack.top() = sum;
        }
        if(tree->get_sibling(i) == -1)      // last child node
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
        trace_length_(j) += tbl;
    }

    double val = sqrt(1.0 / static_cast<double>(size));
    for(int j = 0; j < size; ++j)
    {
        S_values_(values_idx) = trace_length_[j];
        Q_values_(values_idx) = val;
        indices_(values_idx) = static_cast<py::ssize_t>(j);
        values_idx++;
    }
    indptr_(++indptr_idx) = values_idx;

    // element-wise multiply trace length by wavelets
    for(py::ssize_t i = 0; i < nnz; ++i)
    {
        S_values_(i) = S_values_[i] * Q_values_[i];
    }

    // sparse matrix-matrix multiply
    py::ssize_t nnz_max = compute_nnz_max_cov(tree);
    py::tuple S = dspmtm(Q_values_, indices_, indptr_, S_values_, indices_, indptr_, nnz_max);

    py::tuple Q = py::make_tuple(Q_values, indices, indptr);

    return py::make_tuple(Q, S);
}