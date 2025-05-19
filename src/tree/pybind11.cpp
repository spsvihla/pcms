// Project-specific includes
#include "tree-dist.hpp"
#include "tree.hpp"

// Pybind11 includes
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;


// Pybind11 code to expose the Tree class to Python
PYBIND11_MODULE(tree, m) {
    py::class_<Tree>(m, "Tree")
        .def(
            py::init<int>(), 
            py::arg("n_nodes")
        )
        .def
        (
            "get_size", 
            &Tree::get_size
        )
        .def
        (
            "get_parent", 
            &Tree::get_parent, 
            py::arg("u")
        )
        .def
        (
            "get_child", 
            &Tree::get_child, 
            py::arg("u")
        )
        .def
        (
            "get_sibling", 
            &Tree::get_sibling, 
            py::arg("u")
        )
        .def
        (
            "get_subtree_size", 
            [](const Tree& tree, std::optional<int> u) -> py::object {
                if(u.has_value())
                {
                    return py::cast(tree.get_subtree_size(u.value()));
                }
                else
                {
                    const std::vector<int> arr = tree.get_subtree_size();
                    return py::array_t<int>(arr.size(), arr.data());
                }
            },
            py::arg("u") = std::nullopt
        )
        .def
        (
            "get_edge_length", 
            [](const Tree& tree, std::optional<int> u) -> py::object {
                if(u.has_value())
                {
                    return py::cast(tree.get_edge_length(u.value()));
                }
                else
                {
                    const std::vector<double> arr = tree.get_edge_length();
                    return py::array_t<double>(arr.size(), arr.data());
                }
            },
            py::arg("u") = std::nullopt
        )
        .def
        (
            "set_edge_length", 
            &Tree::set_edge_length, 
            py::arg("u"), 
            py::arg("value")
        )
        .def
        (
            "get_name", 
            &Tree::get_name, 
            py::arg("u")
        )
        .def
        (
            "set_name", 
            &Tree::set_name, 
            py::arg("u"), 
            py::arg("value")
        )
        .def
        (
            "link", 
            &Tree::link, 
            py::arg("u"), 
            py::arg("v")
        )
        .def
        (
            "cut", 
            &Tree::cut, 
            py::arg("u")
        )
        .def
        (
            "swap", 
            &Tree::swap, 
            py::arg("u"), 
            py::arg("v")
        )
        .def
        (
            "find_children",
            [](const Tree& tree, int u) {
                const std::vector<int>& children = tree.find_children(u);
                return py::array_t<int>(children.size(), children.data());
            }, 
            py::arg("u")
        )
        .def
        (
            "find_ancestors", 
            [](const Tree& tree, int u) {
                const std::vector<int>& ancestors = tree.find_ancestors(u);
                return py::array_t<int>(ancestors.size(), ancestors.data());
            }, 
            py::arg("u")
        )
        .def
        (
            "find_path", 
            [](const Tree& tree, int u, int v) {
                auto [a, b] = tree.find_path(u, v);
                py::array_t<int> arr_a(a.size(), a.data());
                py::array_t<int> arr_b(b.size(), b.data());
                return py::make_tuple(arr_a, arr_b);
            }, 
            py::arg("u"), 
            py::arg("v")
        )
        .def
        (
            "find_root",
            &Tree::find_root
        )
        .def(
            "find_leaves",
            [](const Tree& tree, std::optional<int> u) {
                std::vector<int> a, b;
                if (u.has_value()) {
                    std::tie(a, b) = tree.find_leaves(u.value());
                } else {
                    std::tie(a, b) = tree.find_leaves();
                }
                py::array_t<int> arr_a(a.size(), a.data());
                py::array_t<int> arr_b(b.size(), b.data());
                return py::make_tuple(arr_a, arr_b);
            },
            py::arg("u") = py::none()
        )
        .def
        (
            "find_epl", 
            &Tree::find_epl
        )
        .def(
            "find_tbl",
            [](const Tree& tree, std::optional<int> u, std::optional<int> v) -> py::object {
                if(!u.has_value() && !v.has_value()) 
                {
                    return py::cast(tree.find_tbl());
                } 
                else if(u.has_value() && !v.has_value()) 
                {
                    return py::cast(tree.find_tbl(u.value()));
                } 
                else if(u.has_value() && v.has_value()) 
                {
                    return py::cast(tree.find_tbl(u.value(), v.value()));
                } 
                else 
                {
                    throw std::invalid_argument("Invalid combination of arguments to find_tbl");
                }
            },
            py::arg("u") = std::nullopt,
            py::arg("v") = std::nullopt
        )
        .def
        (
            "print",
            &Tree::print,
            py::arg("label") = "none"
        );
    py::class_<critical_beta_splitting_distribution>(m, "CriticalBetaSplittingDistribution")
        .def(
            py::init<int>(), 
            py::arg("n")
        )
        .def(
            "__call__", 
            [](critical_beta_splitting_distribution& dist) {
                static thread_local std::mt19937 rng(std::random_device{}());
                return dist(rng);
            }
        )
        .def(
            "get_pmf",
            &critical_beta_splitting_distribution::get_pmf
        )
        .def(
            "get_cdf",
            &critical_beta_splitting_distribution::get_cdf
        );
    m.def
    (
        "nwk2tree",
        &nwk2tree,
        py::arg("filename"),
        py::return_value_policy::take_ownership    
    );
    m.def
    (
        "remy", 
        &remy,
        py::arg("n_leaves"), 
        py::arg("seed") = 0,
        py::return_value_policy::take_ownership
    );
    m.def
    (
        "cbst",
        &cbst,
        py::arg("n_leaves"),
        py::arg("seed") = 0,
        py::return_value_policy::take_ownership
    );
}