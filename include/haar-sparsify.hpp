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
 * @brief Computes the Haar-like wavelet basis and sparse covariance matrix of 
 *        a tree.
 * @param tree
 * @return A pair of tuples with the compressed sparse column representation
 *         of the basis matrix Q and the matrix SQ, where S is the dense
 *         covariance matrix.
 */
py::tuple sparsify(Tree* tree);

#endif // HAAR_SPARSIFY_HPP