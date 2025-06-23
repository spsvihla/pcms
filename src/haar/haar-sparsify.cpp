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
        nnz += tree->get_subtree_size(children_[0]) * (n_children - 1);
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
    int n = tree->get_n_nodes();
    int epl = tree->find_epl();
    return 3 / n - 2 * (epl + 1) / (n * n);
}

// double precision dot product
inline double
ddot(const py::detail::unchecked_reference<double, 1L>& arr1_,
     const py::detail::unchecked_reference<double, 1L>& arr2_,
     py::ssize_t size)
{
    double sum = 0.0;
    for(py::ssize_t i = 0; i < size; ++i)
    {
        sum += arr1_[i] * arr2_[i];
    }
    return sum;
}

// double precision sparse matrix multiply
// C <- A^T * B
// TODO: modify B <- A^T * B in place
inline py::tuple
dspmm_transp(const py::detail::unchecked_reference<double, 1L>& A_values_, 
             const py::detail::unchecked_reference<py::ssize_t, 1L>& A_rowidx_, 
             const py::detail::unchecked_reference<py::ssize_t, 1L>& A_colptr_, 
             const py::detail::unchecked_reference<double, 1L>& B_values_,
             const py::detail::unchecked_reference<py::ssize_t, 1L>& B_rowidx_,
             const py::detail::unchecked_reference<py::ssize_t, 1L>& B_colptr_,
             int nnz_max)
{
    int n = A_colptr_.size() - 1;  // A has n columns
    int m = B_colptr_.size() - 1;  // B has m columns

    // allocate memory for new matrix
    py::array_t<double> values(nnz_max);
    py::array_t<py::ssize_t> rowidx(nnz_max);
    py::array_t<py::ssize_t> colptr(m);

    auto values_ = values.mutable_unchecked<1>();
    auto rowidx_ = rowidx.mutable_unchecked<1>();
    auto colptr_ = colptr.mutable_unchecked<1>();

    // sparse matrix multiplication
    int idx = 0;
    int colptr_idx = 0;
    int b_start = 0;
    for(int i = 1; i < m; ++i)
    {
        int b_end = B_colptr_[i];
        int a_start = 0;
        for(int j = 1; j < n; ++j)
        {
            int a_end = A_colptr_[j];

            if(b_end < a_start || a_end < b_start)
            {
                a_start = a_end;
                continue;
            }
            
            int start = std::max(a_start, b_start);
            int size = std::min(a_end, b_end) - start;

            double val = ddot(A_values_, B_values_, size);

            if(val != 0) 
            {
                values_(idx) = val;
                rowidx_(idx) = j;
                idx++;
            }
            a_start = a_end;
        }
        colptr_(colptr_idx++) = idx;
        b_start = b_end;
    }

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
    py::array_t<py::ssize_t> col_ptrs(tree->find_n_leaves());
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
        else
        {
            int start = subtree_start_stack.top();
            int rsize = tree->get_subtree_size(i);
            int lsize = subtree_size_stack.top();
            int sum = lsize + rsize;
            double tbl = tree->find_tbl(i);
            double lval = sqrt(rsize / lsize / sum);
            double rval = sqrt(lsize / rsize / sum);
            for(int j = 0; j < lsize; ++j)
            {
                trlength_(n_idx) = trlength_[n_idx] + tbl;
                row_idxs_(n_idx) = static_cast<py::ssize_t>(start + j);
                wavelets_(n_idx) = lval;
                n_idx++;
            }
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
        if(tree->get_sibling(i) == -1)      // leaf node
        {
            subtree_start_stack.pop();
            subtree_size_stack.pop();
        }
    }

    // element-wise multiply trace length by wavelets
    for(py::ssize_t i = 0; i < nnz; ++i)
    {
        trlength_(i) = trlength_[i] * wavelets_[i];
    }

    // sparse matrix-matrix multiply
    py::ssize_t nnz_max = compute_nnz_max_cov(tree);
    py::tuple S = dspmm_transp(wavelets_, row_idxs_, col_ptrs_, trlength_, row_idxs_, col_ptrs_, nnz_max);

    py::tuple Q = py::make_tuple(wavelets, row_idxs, col_ptrs);
    
    return py::make_tuple(Q, S);
}