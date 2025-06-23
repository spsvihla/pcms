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
 * @brief 
 * @param tree 
 * @return 
 */
py::tuple sparsify(Tree* tree);

#endif // HAAR_SPARSIFY_HPP