// Project-specific includes
#include "tree-dist.hpp"
#include "tree.hpp"

// Pybind11 includes
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;


// Pybind11 code to expose the module to Python
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
        .def
        (
            "find_is_planted",
            &Tree::find_is_planted
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
        .def(
            "find_n_leaves",
            &Tree::find_n_leaves
        )
        .def
        (
            "find_epl", 
            &Tree::find_epl
        )
        .def(
            "find_tbl",
            [](const Tree& tree, py::object u, py::object v) -> py::object {
                if(u.is_none() && v.is_none()) 
                {
                    return tree.find_tbl();
                } 
                else if(!u.is_none() && v.is_none()) 
                {
                    return py::float_(tree.find_tbl(u.cast<int>()));
                } 
                else if(!u.is_none() && !v.is_none()) 
                {
                    return py::float_(tree.find_tbl(u.cast<int>(), v.cast<int>()));
                } 
                else 
                {
                    throw std::invalid_argument("Invalid combination of arguments to find_tbl");
                }
            },
            py::arg("u"),
            py::arg("v")
        )
        .def
        (
            "compute_wavelets", 
            &Tree::compute_wavelets, 
            py::arg("u")
        )
        .def
        (
            "compute_supports", 
            &Tree::compute_supports, 
            py::arg("u")
        )
        .def
        (
            "to_string",
            &Tree::to_string,
            py::arg("label")
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
        py::arg("ensure_planted"),
        py::return_value_policy::take_ownership    
    );
    m.def
    (
        "remy", 
        &remy,
        py::arg("n_leaves"), 
        py::arg("planted"),
        py::arg("seed"),
        py::return_value_policy::take_ownership
    );
    m.def
    (
        "cbst",
        &cbst,
        py::arg("n_leaves"),
        py::arg("planted"),
        py::arg("seed"),
        py::return_value_policy::take_ownership
    );
}