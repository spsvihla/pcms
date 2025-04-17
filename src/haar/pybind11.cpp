// Project-specific includes
#include "tree.hpp"
#include "haar-dist.hpp"
#include "haar-sparsify.hpp"

// Pybind11 includes
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;


py::tuple
SparseMatrix2py_tuple(SparseMatrix<double>& S)
{
    const std::vector<double>& values = S.values;
    const std::vector<int>& rowidx = S.rowidx;
    const std::vector<int>& colptr = S.colptr;

    // cast int -> py::ssize_t
    std::vector<py::ssize_t> py_ssize_t_rowidx(rowidx.begin(), rowidx.end());
    std::vector<py::ssize_t> py_ssize_t_colptr(colptr.begin(), colptr.end());

    // cast std::vector -> py::array_t
    py::array_t<double> py_values(
        values.size(),
        values.data()
    );
    py::array_t<py::ssize_t> py_rowidx(
        py_ssize_t_rowidx.size(),
        py_ssize_t_rowidx.data()
    );
    py::array_t<py::ssize_t> py_colptr(
        py_ssize_t_colptr.size(),
        py_ssize_t_colptr.data()
    );

    return py::make_tuple(py_values, py_rowidx, py_colptr);
}

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
        [](Tree& tree) -> py::object {
            auto [Q, S] = sparsify(&tree);
            py::tuple Q_py_tuple = SparseMatrix2py_tuple(Q);
            py::tuple S_py_tuple = SparseMatrix2py_tuple(S);
            return py::make_tuple(Q_py_tuple, S_py_tuple);
        },
        py::arg("tree"),
        py::return_value_policy::take_ownership
    );
}