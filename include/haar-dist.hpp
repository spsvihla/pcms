/**
 * @file haar-dist.hpp
 * @brief Functions for computing the CDF of the Haar-like coordinate ⟨f,φ⟩.
 * 
 * This header defines the interface for estimating the cumulative distribution 
 * function (CDF) of the random variable ⟨f,φ⟩ where either 'f' or 'φ' is 
 * allowed to vary randomly.
 * 
 * @author Sean Svihla
 */

 #ifndef HAAR_DIST_H
 #define HAAR_DIST_H
 
 // Standard library includes
 #include <vector>
 
// project-specific includes
#include "tree.hpp"

 // Pybind11 includes for interoperability with NumPy and Python
 #include <pybind11/numpy.h>
 #include <pybind11/pybind11.h>
 
 namespace py = pybind11;
 
 
 /**
  * @brief Sample the Haar-like component λ_v(Δ_v^2) where T is drawn from
  *        the critical beta-splitting distribution with Exp(|L(v)|) edge
  *        lengths and f ~ Uniform(f(Π_n)), the set of permutations of f.
  * 
  * @param f A 1D NumPy array representing a function on the leaves of T(v).
  * @param n_samples The number of samples to draw.
  * @param seed Random seed for reproducibility. If zero, a random seed is used.
  * @return A 1D NumPy array of the same length as `ys`, containing estimated 
  *         CDF values.
  */
 py::array_t<double> rand_dh_component(const py::array_t<double>& f, 
                                       int n_samples, 
                                       std::vector<Tree*>& buffer,
                                       std::optional<unsigned int> seed);

 #endif // HAAR_DIST_H 