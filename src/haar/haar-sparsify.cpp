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
    return static_cast<py::ssize_t>(ceil(frac * n * (n - 1)));
}

// double precision sparse dot product
inline double 
dsdot(
    const py::detail::unchecked_reference<double, 1L>& A_values_,
    const py::detail::unchecked_reference<py::ssize_t, 1L>& A_rowidx_,
    py::ssize_t a_start, py::ssize_t a_end,
    const py::detail::unchecked_reference<double, 1L>& B_values_,
    const py::detail::unchecked_reference<py::ssize_t, 1L>& B_rowidx_,
    py::ssize_t b_start, py::ssize_t b_end)
{
    double result = 0.0;
    py::ssize_t i = a_start, j = b_start;
    while(i < a_end && j < b_end)
    {
        py::ssize_t a_row = A_rowidx_[i];
        py::ssize_t b_row = B_rowidx_[j];
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
// TODO: modify B <- A^T * B in place
inline py::tuple
dspmtm(const py::detail::unchecked_reference<double, 1L>& A_values_, 
       const py::detail::unchecked_reference<py::ssize_t, 1L>& A_rowidx_, 
       const py::detail::unchecked_reference<py::ssize_t, 1L>& A_colptr_, 
       const py::detail::unchecked_reference<double, 1L>& B_values_,
       const py::detail::unchecked_reference<py::ssize_t, 1L>& B_rowidx_,
       const py::detail::unchecked_reference<py::ssize_t, 1L>& B_colptr_,
       py::ssize_t nnz_max)
{
    // allocate memory for new matrix
    py::array_t<double> values(nnz_max);
    py::array_t<py::ssize_t> rowidx(nnz_max);
    py::array_t<py::ssize_t> colptr(B_colptr_.size());

    auto values_ = values.mutable_unchecked<1>();
    auto rowidx_ = rowidx.mutable_unchecked<1>();
    auto colptr_ = colptr.mutable_unchecked<1>();

    // sparse matrix multiplication
    py::ssize_t idx = 0;
    py::ssize_t colptr_idx = 0;
    for(py::ssize_t i = 0; i < B_colptr_.size() - 1; ++i)
    {
        colptr_(colptr_idx++) = idx;

        py::ssize_t b_start_idx = B_colptr_[i];
        py::ssize_t b_end_idx = B_colptr_[i+1];
        py::ssize_t b_start_row = B_rowidx_[b_start_idx];
        py::ssize_t b_end_row = B_rowidx_[b_end_idx];

        for(py::ssize_t j = 0; j < A_colptr_.size() - 1; ++j)
        {
            py::ssize_t a_start_idx = A_colptr_[j];
            py::ssize_t a_end_idx = A_colptr_[j+1];
            py::ssize_t a_start_row = A_rowidx_[a_start_idx];
            py::ssize_t a_end_row = A_rowidx_[a_end_idx];

            // check for disjoint support
            if(b_end_row < a_start_row || a_end_row < b_start_row)
            {
                continue;
            }

            // compute dot product
            double val = dsdot(A_values_, A_rowidx_, a_start_idx, a_end_idx,
                               B_values_, B_rowidx_, b_start_idx, b_end_idx);
            if(val != 0) 
            {
                if(idx >= nnz_max)
                {
                    throw std::runtime_error(
                        "Exceeded allocated nnz_max in sparse multiply: idx = " 
                        + std::to_string(idx) + ", nnz_max = " 
                        + std::to_string(nnz_max)
                    );
                }
                values_(idx) = val;
                rowidx_(idx) = j;
                idx++;
            }
        }
    }
    colptr_(colptr_idx) = idx;

    values.resize({idx});
    rowidx.resize({idx});

    return py::make_tuple(values, rowidx, colptr);
}

py::tuple
sparsify(Tree* tree)
{
    py::ssize_t nnz = compute_nnz_wavelets(tree);

    py::array_t<int> subtree_starts = tree->find_subtree_start_indices();
    auto subtree_starts_ = subtree_starts.unchecked<1>();

    py::array_t<double> wavelets(nnz);
    py::array_t<double> trlength(nnz);
    py::array_t<py::ssize_t> col_ptrs(tree->find_n_leaves()+1);
    py::array_t<py::ssize_t> row_idxs(nnz);

    auto wavelets_ = wavelets.mutable_unchecked<1>();
    auto trlength_ = trlength.mutable_unchecked<1>();
    auto col_ptrs_ = col_ptrs.mutable_unchecked<1>();
    auto row_idxs_ = row_idxs.mutable_unchecked<1>();

    for(py::ssize_t i = 0; i < nnz; ++i)
    {
        trlength_(i) = 0.0;
    }

    col_ptrs_(0) = 0;

    // build wavelet basis and trace length sparse matrices
    py::ssize_t n_idx = 0;
    py::ssize_t m_idx = 0;

    std::stack<int> subtree_start_stack;
    std::stack<int> subtree_size_stack;
    for(int i = 0; i < tree->get_n_nodes(); ++i)
    {
        if(tree->get_is_first(i))           // first child node
        {
            subtree_start_stack.push(subtree_starts_[static_cast<py::ssize_t>(i)]);
            subtree_size_stack.push(tree->get_subtree_size(i));
        }
        else                                // daughter wavelet
        {
            int start = subtree_start_stack.top();

            int rsize = tree->get_subtree_size(i);
            int lsize = subtree_size_stack.top();
            int sum = lsize + rsize;

            double tbl = tree->find_tbl(i);

            // use floating-point arithmetic to avoid overflow with large integers
            double lval = sqrt(static_cast<double>(rsize) / (static_cast<double>(lsize) * static_cast<double>(sum)));
            double rval = -1 * sqrt(static_cast<double>(lsize) / (static_cast<double>(lsize) * static_cast<double>(sum)));

            // left subtree
            for(int j = 0; j < lsize; ++j)
            {
                trlength_(n_idx) = trlength_[n_idx] + tbl;
                row_idxs_(n_idx) = static_cast<py::ssize_t>(start + j);
                wavelets_(n_idx) = lval;
                n_idx++;
            }

            // right subtree
            for(int j = lsize; j < sum; ++j)
            {
                trlength_(n_idx) = trlength_[n_idx] + tbl;
                row_idxs_(n_idx) = static_cast<py::ssize_t>(start + j);
                wavelets_(n_idx) = rval;
                n_idx++;
            }

            col_ptrs_(++m_idx) = n_idx;
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
    double val = sqrt(1.0 / static_cast<double>(size));
    for(int j = 0; j < size; ++j)
    {
        trlength_(n_idx) = 0; // BUG
        row_idxs_(n_idx) = static_cast<py::ssize_t>(j);
        wavelets_(n_idx) = val;
        n_idx++;
    }
    col_ptrs_(++m_idx) = n_idx;

    // element-wise multiply trace length by wavelets
    for(py::ssize_t i = 0; i < nnz; ++i)
    {
        trlength_(i) = trlength_[i] * wavelets_[i];
    }

    // sparse matrix-matrix multiply
    py::ssize_t nnz_max = compute_nnz_max_cov(tree);
    py::tuple S = dspmtm(wavelets_, row_idxs_, col_ptrs_, trlength_, row_idxs_, col_ptrs_, nnz_max);

    py::tuple Q = py::make_tuple(wavelets, row_idxs, col_ptrs);
    
    return py::make_tuple(Q, S);
}