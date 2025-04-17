/**
 * @file _sparse.hpp
 * @brief
 * @author Sean Svihla
 */
#ifndef _SPARSE_HPP
#define _SPARSE_HPP

// project-specific includes
#include "tree.hpp"

template <typename ValueType>
struct SparseMatrix {
    std::vector<ValueType> values;
    std::vector<int> rowidx;
    std::vector<int> colptr;
};

/**
 * @brief 
 * @param tree 
 * @return 
 */
std::pair<SparseMatrix<double>, SparseMatrix<double>> sparsify(Tree* tree);

#endif // _SPARSE_HPP