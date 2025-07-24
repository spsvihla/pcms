// Project-specific includes
#include "tree.hpp"
#include "haar-dist.hpp"
#include "haar-sparsify.hpp"

// Pybind11 includes
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;


// Pybind11 code to expose the module to Python
PYBIND11_MODULE(_haar, m) {
    m.def
    (
        "sample_dh_component",
        &sample_dh_component,
        py::arg("func"),
        py::arg("n_samples"),
        py::arg("seed")
    );
    m.def
    (
        "sparsify",
        &sparsify,
        py::arg("tree")
    );
}