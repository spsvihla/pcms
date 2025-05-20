/**
 * @file haar-sparsify.hpp
 * @brief
 * @author Sean Svihla
 */
#ifndef HAAR_SPARSIFY_HPP
#define HAAR_SPARSIFY_HPP

// project-specific includes
#include "tree.hpp"

// Pybind11 includes
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

/**
 * @struct SparseMatrix
 * @brief 
 * @tparam ValueType 
 */
template <typename ValueType>
struct SparseMatrix {
    py::array_t<ValueType> values;
    py::array_t<py::ssize_t> rowidx;
    py::array_t<py::ssize_t> colptr;
};

/**
 * @brief 
 * @param tree 
 * @return 
 */
std::pair<SparseMatrix<double>, SparseMatrix<double>> sparsify(Tree* tree);

#endif // HAAR_SPARSIFY_HPP