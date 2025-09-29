// Project-specific includes
#include "tree.hpp"
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
        "sparsify",
        &sparsify,
        py::arg("tree")
    );
}