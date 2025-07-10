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
 
 // Pybind11 includes for interoperability with NumPy and Python
 #include <pybind11/numpy.h>
 #include <pybind11/pybind11.h>
 
 namespace py = pybind11;
 
 /**
  * @brief Estimates the CDF values of ⟨f,φ⟩ at points 'ys', where φ~CBST(n) is
  *        a Haar-like wavelet from a critical beta-splitting random tree.
  * 
  * This function estimates the value of the CDF using a stochastic estimator
  * in order to approximate summation over all permutations of 'f'. The value
  * of 'num_iter' should be chosen to ensure the desired convergence.
  * 
  * @param ys A 1D NumPy array of points at which to evaluate the estimated CDF.
  * @param func A 1D NumPy array representing the function f to project against 
  *             the random basis.
  * @param pmf A 1D NumPy array representing the CBST PMF.
  * @param num_iter Number of Monte Carlo samples to use for estimation.
  * @param seed Random seed for reproducibility. If zero, a random seed is used.
  * @return A 1D NumPy array of the same length as `ys`, containing estimated 
  *         CDF values.
  */
 py::array_t<double> cdf_cbst_topology(const py::array_t<double>& ys, 
                                       const py::array_t<double>& func, 
                                       const py::array_t<double>& pmf, 
                                       int num_iter, 
                                       std::optional<unsigned int> seed);

 #endif // HAAR_DIST_H 