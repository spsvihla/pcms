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

// external includes
#include <mkl.h>

// project-specific includes
#include "tree.hpp"
#include "haar-sparsify.hpp"

// pybind11 and numpy includes
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;


// compute number of non-zero entries in wavelet basis matrix
int
find_nnz(Tree* tree)
{
    int nnz = 0;
    for(int i = 0; i < tree->get_n_nodes(); ++i)
    {
        if(tree->get_child(i) == -1L)
        {
            continue;
        }
        std::vector<int> children = tree->find_children_(i);
        int n_children = static_cast<int>(children.size());
        nnz += tree->get_subtree_size(children[0]) * std::max(1, n_children - 1);
        for(int j = 1; j < n_children; ++j)
        {
            nnz += tree->get_subtree_size(children[j]) * (n_children - j);
        }
    }
    return nnz;
}

py::tuple mkl2py_csc(sparse_matrix_t A) 
{
    sparse_index_base_t indexing;
    MKL_INT rows, cols;
    MKL_INT *col_start = nullptr;
    MKL_INT *col_end = nullptr;
    MKL_INT *row_indx = nullptr;
    double *values = nullptr;

    // Export CSC structure from MKL
    mkl_sparse_d_export_csc(
        A, 
        &indexing, 
        &rows, 
        &cols, 
        &col_start, 
        &col_end, 
        &row_indx, 
        &values
    );

    MKL_INT nnz = col_end[cols - 1];

    // Build indptr array from col_start and last col_end element
    std::vector<MKL_INT> indptr(cols + 1);
    #pragma omp simd
    for(MKL_INT i = 0; i < cols; ++i) 
    {
        indptr[i] = col_start[i];
    }
    indptr[cols] = col_end[cols - 1];

    // Deep copy the values and indices into new Python-owned arrays
    auto values_array = py::array_t<double>(nnz);
    std::memcpy(values_array.mutable_data(), values, sizeof(double) * nnz);

    auto indices_array = py::array_t<MKL_INT>(nnz);
    std::memcpy(indices_array.mutable_data(), row_indx, sizeof(MKL_INT) * nnz);

    // indptr is already built in a std::vector, so we copy it too
    auto indptr_array = py::array_t<MKL_INT>(indptr.size());
    std::memcpy(indptr_array.mutable_data(), indptr.data(), sizeof(MKL_INT) * indptr.size());

    return py::make_tuple(values_array, indices_array, indptr_array);
}

py::tuple
sparsify(Tree* tree)
{
    int nnz = find_nnz(tree);
    int n_leaves = tree->find_n_leaves();
    int n_wavelets = tree->find_n_wavelets();
    int n_nodes = tree->get_n_nodes();

    std::vector<int> subtree_starts = tree->find_subtree_start_indices_();

    MKL_INT rows = n_leaves;
    MKL_INT cols = n_wavelets;

    std::vector<double> Q_values(nnz);
    std::vector<MKL_INT> Q_indices(nnz);
    std::vector<MKL_INT> Q_indptr(cols + 1);

    std::vector<double> S_values(nnz);
    std::vector<MKL_INT> S_indices(nnz);
    std::vector<MKL_INT> S_indptr(cols + 1);

    std::vector<double> trace_length(n_leaves);
    #pragma omp simd
    for(int i = 0; i < n_leaves; ++i)
    {
        trace_length[i] = 0.0;
    }

    int idx = 0;
    int col = 0;

    std::stack<int> subtree_start_stack;
    std::stack<int> subtree_size_stack;

    Q_indptr[0] = 0;
    S_indptr[0] = 0;

    for(int i = 0; i < n_nodes - 1; ++i)
    {
        // accumulate trace branch length
        int start = subtree_starts[i];
        int size = tree->get_subtree_size(i);
        double tbl = tree->find_tbl(i);
        #pragma omp simd
        for(int j = 0; j < size; ++j)
        {
            trace_length[start + j] += tbl;
        }

        // construct wavelet
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

            double rsize_ = static_cast<double>(rsize);
            double lsize_ = static_cast<double>(lsize);
            double sum_ = static_cast<double>(sum);

            double lval =  sqrt(rsize_ / (lsize_ * sum_));
            double rval = -sqrt(lsize_ / (rsize_ * sum_));

            // left subtree
            #pragma omp simd
            for(int j = 0; j < lsize; ++j)
            {
                Q_values[idx] = lval;
                S_values[idx] = lval * trace_length[start + j];
                Q_indices[idx] = start + j;
                S_indices[idx] = start + j;
                idx++;
            }

            // right subtree
            #pragma omp simd
            for(int j = lsize; j < sum; ++j)
            {
                Q_values[idx] = rval;
                S_values[idx] = rval * trace_length[start + j];
                Q_indices[idx] = start + j;
                S_indices[idx] = start + j;
                idx++;
            }

            col++;
            Q_indptr[col] = idx;
            S_indptr[col] = idx;

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

    double val = sqrt(1.0 / size);
    #pragma omp simd
    for(int i = 0; i < size; ++i)
    {
        Q_values[idx] = val;
        S_values[idx] = val * trace_length[i];
        Q_indices[idx] = i;
        S_indices[idx] = i;
        idx++;
    }

    col++;
    Q_indptr[col] = idx;
    S_indptr[col] = idx;

    sparse_matrix_t Q;
    mkl_sparse_d_create_csc(
        &Q,
        SPARSE_INDEX_BASE_ZERO,
        rows,
        cols,
        Q_indptr.data(),
        Q_indptr.data() + 1,
        Q_indices.data(),
        Q_values.data()
    );

    sparse_matrix_t R;
    mkl_sparse_d_create_csc(
        &R,
        SPARSE_INDEX_BASE_ZERO,
        rows,
        cols,
        S_indptr.data(),
        S_indptr.data() + 1,
        S_indices.data(),
        S_values.data()
    );

    sparse_matrix_t S;
    mkl_sparse_spmm(SPARSE_OPERATION_TRANSPOSE, Q, R, &S);

    auto Q_py = mkl2py_csc(Q);
    auto S_py = mkl2py_csc(S);

    mkl_sparse_destroy(Q);
    mkl_sparse_destroy(R);
    mkl_sparse_destroy(S);

    return py::make_tuple(Q_py, S_py);
}
