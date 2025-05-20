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
#include <queue>
#include <stack>
#include <vector>

// project-specific includes
#include "tree.hpp"
#include "haar-sparsify.hpp"

// pybind11 and numpy includes
#include <pybind11/pybind11.h>

namespace py = pybind11;


// TODO: optimize linear algebra functions

// double precision element-wise vector-plus-constant
// a <- a + c
inline void
dewvpc(double* ptr, int size, double value)
{
    for(int i = 0; i < size; ++i)
    {
        ptr[i] += value;
    }
}

// double precision element-wise vector-multuply-vector
// b <- a * b 
inline void
dewvmv(double* a, double* b, int size)
{
    for(int i = 0; i < size; ++i)
    {
        b[i] *= a[i];
    }
}

// double precision dot product
inline double
ddot(double* a, double* b, int size)
{
    double sum = 0.0;
    for(int i = 0; i < size; ++i)
    {
        sum += a[i] * b[i];
    }
    return sum;
}

// TODO: modify B <- A^T * B in place

// double precision sparse matrix multiply
// C <- A^T * B
// inline SparseMatrix<double>
// dspmm_transp(SparseMatrix<double>& A, SparseMatrix<double>& B, int n, int m, int nnz_max)
// {
//     std::vector<double> values;
//     std::vector<int> rowidx;
//     std::vector<int> colptr;

//     values.reserve(nnz_max);
//     rowidx.reserve(nnz_max);
//     colptr.reserve(m);

//     int idx = 0;
//     int b_start = 0;
//     for(int i = 1; i < m; ++i)
//     {
//         int b_end = B.colptr[i];
//         int a_start = 0;
//         for(int j = 1; j < n; ++j)
//         {
//             int a_end = A.colptr[j];

//             if(b_end < a_start || a_end < b_start)
//             {
//                 a_start = a_end;
//                 continue;
//             }
            
//             int start = std::max(a_start, b_start);
//             int size = std::min(a_end, b_end) - start;

//             double val = ddot(&A.values[start], &B.values[start], size);

//             if(val != 0) 
//             {
//                 values.push_back(val);
//                 rowidx.push_back(j);
//                 idx++;
//             }
//             a_start = a_end;
//         }
//         colptr.push_back(idx);
//         b_start = b_end;
//     }

//     return SparseMatrix<double>{values, rowidx, colptr};
// }

py::ssize_t
compute_nnz(Tree* tree)
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

py::array_t<int>
compute_subtree_starts(Tree* tree) 
{
    py::array_t<int> subtree_starts(tree->get_n_nodes());
    auto subtree_starts_ = subtree_starts.mutable_unchecked<1>();

    // breadth-first search
    std::queue<std::pair<int, int>> q;
    q.push({tree->find_root(), 0});
    while(!q.empty())
    {
        auto [u, start] = q.front();
        q.pop();

        subtree_starts_(u) = start;

        // enqueue children
        int offset = 0;
        for(int c = tree->get_child(u); c != -1; c = tree->get_sibling(c))
        {
            q.push({c, start + offset});
            offset += tree->get_subtree_size(c);
        }
    }

    return subtree_starts;
}

std::pair<SparseMatrix<double>, SparseMatrix<double>> 
sparsify(Tree* tree)
{
    int nnz = compute_nnz(tree);

    py::array_t<int> subtree_starts = compute_subtree_starts(tree);
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

    py::ssize_t n_idx = 0;
    py::ssize_t m_idx = 0;

    std::stack<int> subtree_start_stack;
    std::stack<int> subtree_size_stack;
    for(int i = 0; i < tree->get_n_nodes(); ++i)
    {
        if(tree->get_is_first(i))
        {
            subtree_start_stack.push(subtree_starts_[i]);
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
            for(int j = lsize; j < rsize; ++j)
            {
                trlength_(n_idx) = trlength_[n_idx] + tbl;
                row_idxs_(n_idx) = static_cast<py::ssize_t>(start + j);
                wavelets_(n_idx) = rval;
                n_idx++;
            }
            col_ptrs_(++m_idx) = n_idx;
            subtree_size_stack.top() = sum;
        }
        if(tree->get_sibling(i) == -1)
        {
            subtree_start_stack.pop();
            subtree_size_stack.pop();
        }
    }

    // TODO: multiple wavelet by trace length

    // TODO: sparse matrix-matrix multiply
}