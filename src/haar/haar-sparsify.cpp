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

// compute upper bound on number of non-zero entries in sparsified covariance matrix
// TODO: generalize bound to non-binary trees
int
compute_nnz_max_cov(Tree* tree)
{
    int n = tree->find_n_leaves();
    int epl = tree->find_epl();
    double frac = 2.0 * (epl + 1) / (n * n) - 3.0 / n;
    return static_cast<int>(ceil(frac * n * n));
}

// double precision sparse matrix-transpose-matrix multiply
// C <- A^T * B
// NOTE: This implementation takes advantage of the unique structure of
//       Haar-like basis matrices; it is not a functional general sparse
//       matrix multiply.
inline py::tuple
dsmtm(const double* A_values, const py::ssize_t* A_indices, 
      const py::ssize_t* A_indptr, const double* B_values,
      const py::ssize_t* B_indices, const py::ssize_t* B_indptr,
      std::size_t nnz_max, std::size_t num_cols)
{
    // allocate output CSR arrays
    py::array_t<double> values(static_cast<py::ssize_t>(nnz_max));
    py::array_t<py::ssize_t> indices(static_cast<py::ssize_t>(nnz_max));
    py::array_t<py::ssize_t> indptr(static_cast<py::ssize_t>(num_cols + 1));

    auto* values_  = static_cast<double*>(values.request().ptr);
    auto* indices_ = static_cast<py::ssize_t*>(indices.request().ptr);
    auto* indptr_  = static_cast<py::ssize_t*>(indptr.request().ptr);

    indptr_[0] = 0;

    std::size_t values_idx = 0;

    std::vector<std::size_t> cols(num_cols);

    // for each column in B (row in B^T)
    for(std::size_t i = 0; i < num_cols; ++i)
    {
        std::size_t cols_size = 0;

        // B's block range
        std::size_t bi0 = B_indptr[i];
        std::size_t bi1 = B_indptr[i + 1];
        std::size_t br0 = B_indices[bi0];
        std::size_t br1 = B_indices[bi1 - 1];

        // find overlapping A blocks
        // TODO: find a faster way to determine these
        for(std::size_t j = 0; j < num_cols; ++j)
        {
            std::size_t ai0 = A_indptr[j];
            std::size_t ai1 = A_indptr[j + 1];
            std::size_t ar0 = A_indices[ai0];
            std::size_t ar1 = A_indices[ai1 - 1];

            bool overlap = (ar0 <= br1) && (br0 <= ar1);

            cols[cols_size] = j;
            cols_size += overlap;
        }

        // for each overlapping pair: dense dot product
        for(std::size_t idx = 0; idx < cols_size; ++idx)
        {
            std::size_t j = cols[idx];

            std::size_t ai0 = A_indptr[j];
            std::size_t ai1 = A_indptr[j + 1];
            std::size_t ar0 = A_indices[ai0];
            std::size_t ar1 = A_indices[ai1 - 1];

            std::size_t start_row = std::max(ar0, br0);
            std::size_t end_row = std::min(ar1, br1);
            std::size_t overlap_width = end_row - start_row + 1;

            std::size_t offset_a = start_row - ar0;
            std::size_t offset_b = start_row - br0;

            double val = 0.0;
            for(std::size_t k = 0; k < overlap_width; ++k)
            {
                val += A_values[ai0 + offset_a + k] * B_values[bi0 + offset_b + k];
            }

            if(val != 0.0)
            {
                values_[values_idx] = val;
                indices_[values_idx] = j;
                ++values_idx;
            }
        }

        indptr_[i + 1] = values_idx;
    }

    // shrink output arrays
    values.resize({static_cast<py::ssize_t>(values_idx)});
    indices.resize({static_cast<py::ssize_t>(values_idx)});

    return py::make_tuple(values, indices, indptr);
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

    // sparse matrix-matrix multiply
    int nnz_max = compute_nnz_max_cov(tree);
    py::tuple S = dsmtm(Q_values_, indices_, indptr_, S_values_, indices_, indptr_, nnz_max, n_leaves);

    py::tuple Q = py::make_tuple(Q_values, indices, indptr);

    return py::make_tuple(Q, S);
}