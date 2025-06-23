// Project-specific includes
#include "tree.hpp"
#include "haar-dist.hpp"
#include "haar-sparsify.hpp"

// Pybind11 includes
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;


// pybind11 code to expose the module to Python
PYBIND11_MODULE(_haar, m) {
    m.def
    (
        "cdf_rand_basis",
        &cdf_rand_basis,
        py::arg("ys"),
        py::arg("func"),
        py::arg("pmf"),
        py::arg("num_iter"),
        py::arg("seed") = 0
    );
    m.def
    (
        "sparsify",
        &sparsify,
        py::arg("tree")
    );
}