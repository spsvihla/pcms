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


// TODO: optimize linear algebra functions

// double precision element-wise vector-plus-vector
// b <- a + b
inline void
dewvpv(double* ptr, int size, double value)
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
inline SparseMatrix<double>
dspmm_transp(SparseMatrix<double>& A, SparseMatrix<double>& B, int n, int m, int nnz_max)
{
    std::vector<double> values;
    std::vector<int> rowidx;
    std::vector<int> colptr;

    values.reserve(nnz_max);
    rowidx.reserve(nnz_max);
    colptr.reserve(m);

    int idx = 0;
    int b_start = 0;
    for(int i = 1; i < m; ++i)
    {
        int b_end = B.colptr[i];
        int a_start = 0;
        for(int j = 1; j < n; ++j)
        {
            int a_end = A.colptr[j];

            if(b_end < a_start || a_end < b_start)
            {
                a_start = a_end;
                continue;
            }
            
            int start = std::max(a_start, b_start);
            int size = std::min(a_end, b_end) - start;

            double val = ddot(&A.values[start], &B.values[start], size);

            if(val != 0) 
            {
                values.push_back(val);
                rowidx.push_back(j);
                idx++;
            }
            a_start = a_end;
        }
        colptr.push_back(idx);
        b_start = b_end;
    }

    return SparseMatrix<double>{values, rowidx, colptr};
}

inline void 
add_column(const std::vector<double>& tbl,
           const std::vector<std::vector<double>>& hlw, 
           std::vector<double>& tbl_values, std::vector<double>& hlw_values,
           std::vector<int>& rowidx, std::vector<int>& colptr, int& col_counter, 
           long& val_idx, int& start_idx)
{
    for(int j = 0; j < hlw.size(); ++j)
    {
        for(int k = 0; k < static_cast<int>(hlw[j].size()); ++k)
        {
            int idx = k + start_idx;
            tbl_values[val_idx] = tbl[idx];
            hlw_values[val_idx] = hlw[j][idx];
            rowidx[val_idx] = idx;
            val_idx++;
        }
        colptr[++col_counter] = val_idx;
    }
}

inline void
populate_matrices(Tree* tree, std::vector<double> tbl_values, 
                  std::vector<double> hlw_values, std::vector<int> rowidx, 
                  std::vector<int> colptr)
{
    colptr[0] = 0;

    const int n_nodes = tree->get_size();
    const std::vector<int>& subtree_size = tree->get_subtree_size();

    std::vector<double> tbl(n_nodes, 0.0);
    long val_idx = 0;
    int col_counter = 0;
    int leaf_counter = 0;
    std::stack<int> stack;

    for(int i = 0; i < n_nodes; ++i)
    {
        if(tree->get_child(i) == -1)        // leaf node
        {
            stack.push(leaf_counter++);
        }
        else                                // interior node
        {
            add_column(tbl, tree->make_hlw(i), tbl_values, hlw_values, rowidx,
                       colptr, col_counter, val_idx, stack.top());
        }

        // accumulate trace length
        dewvpv(&tbl[stack.top()], subtree_size[i], tree->find_tbl(i));

        if(!tree->is_first(i))              // sibling node
        {
            stack.pop();
        }
    }

    dewvmv(hlw_values.data(), tbl_values.data(), hlw_values.size());
}

std::pair<SparseMatrix<double>, SparseMatrix<double>> 
sparsify(Tree* tree)
{
    // define sarse matrix data vectors
    int nnz = tree->find_nnz();
    std::vector<double> tbl_values(nnz, 0.0);
    std::vector<double> hlw_values(nnz, 0.0);
    std::vector<int> rowidx(nnz);
    std::vector<int> colptr(nnz);

    populate_matrices(tree, tbl_values, hlw_values, rowidx, colptr);

    // compute sparse matrix-matrix product
    int n_nodes = tree->get_size();
    SparseMatrix<double> Q = {hlw_values, rowidx, colptr};
    SparseMatrix<double> S = {tbl_values, rowidx, colptr};
    S = dspmm_transp(Q, S, n_nodes, n_nodes, 1); // TODO: solve for nnz_max

    return std::make_pair(Q, S);
}