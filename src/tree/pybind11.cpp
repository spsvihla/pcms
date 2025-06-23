// Project-specific includes
#include "tree-dist.hpp"
#include "tree.hpp"

// Pybind11 includes
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;


// Pybind11 code to expose the Tree class to Python
PYBIND11_MODULE(_tree, m) {
    py::class_<Tree>(m, "Tree")
        .def(
            py::init<int>(), 
            py::arg("n_nodes")
        )
        .def
        (
            "get_n_nodes", 
            &Tree::get_n_nodes
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
            "get_is_first",
            &Tree::get_is_first,
            py::arg("u")
        )
        .def
        (
            "get_subtree_size", 
            [](const Tree& tree, std::optional<int> u) -> py::object {
                if(u.has_value())
                {
                    return py::int_(tree.get_subtree_size(u.value()));
                }
                else
                {
                    return tree.get_subtree_size();
                }
            },
            py::arg("u") = std::nullopt,
            py::return_value_policy::reference_internal
        )
        .def
        (
            "get_edge_length", 
            [](const Tree& tree, std::optional<int> u) -> py::object {
                if(u.has_value())
                {
                    return py::float_(tree.get_edge_length(u.value()));
                }
                else
                {
                    return tree.get_edge_length();
                }
            },
            py::arg("u") = std::nullopt,
            py::return_value_policy::reference_internal
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
            &Tree::find_children,
            py::arg("u")
        )
        .def
        (
            "find_ancestors", 
            &Tree::find_ancestors,
            py::arg("u")
        )
        .def
        (
            "find_path", 
            &Tree::find_path,
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
            &Tree::find_leaves,
            py::arg("u")
        )
        .def(
            "find_subtree_start_indices",
            &Tree::find_subtree_start_indices
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
                    return tree.find_tbl();
                } 
                else if(u.has_value() && !v.has_value()) 
                {
                    return py::float_(tree.find_tbl(u.value()));
                } 
                else if(u.has_value() && v.has_value()) 
                {
                    return py::float_(tree.find_tbl(u.value(), v.value()));
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
            "to_string",
            &Tree::to_string,
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
        py::arg("newick_string"),
        py::return_value_policy::take_ownership    
    );
    m.def
    (
        "remy", 
        &remy,
        py::arg("n_leaves"), 
        py::arg("seed") = std::nullopt,
        py::return_value_policy::take_ownership
    );
    m.def
    (
        "cbst",
        &cbst,
        py::arg("n_leaves"),
        py::arg("seed") = std::nullopt,
        py::return_value_policy::take_ownership
    );
}