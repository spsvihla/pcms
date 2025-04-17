/**
 * @file haar-sparsify.hpp
 * @brief
 * @author Sean Svihla
 */
#ifndef HAAR_SPARSIFY_HPP
#define HAAR_SPARSIFY_HPP

// project-specific includes
#include "tree.hpp"

/**
 * @struct SparseMatrix
 * @brief 
 * @tparam ValueType 
 */
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

#endif // HAAR_SPARSIFY_HPP