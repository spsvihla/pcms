// Project-specific includes
#include "dist.hpp"
#include "tree.hpp"

// Pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;


// Pybind11 code to expose the Tree class to Python
PYBIND11_MODULE(tree, m) {
    py::class_<Tree>(m, "Tree")
        .def(
            py::init<int>(), 
            py::arg("n_nodes") = 1
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
            static_cast<int (Tree::*)(int) const>(&Tree::get_subtree_size),
            py::arg("u")
        )
        .def
        (
            "get_subtree_sizes", 
            static_cast<const std::vector<int, AlignedAllocator<int, 32>>& (Tree::*)() const>(&Tree::get_subtree_size)
        )
        .def
        (
            "get_edge_length", 
            static_cast<double (Tree::*)(int) const>(&Tree::get_edge_length),
            py::arg("u")
        )
        .def
        (
            "get_edge_lengths", 
            static_cast<const std::vector<double, AlignedAllocator<double, 32>>& (Tree::*)() const>(&Tree::get_edge_length)
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
            "find_support", 
            &Tree::find_support, 
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
            "find_leaves", 
            &Tree::find_leaves
        )
        .def
        (
            "find_epl", 
            &Tree::find_epl
        )
        .def
        (
            "print",
            &Tree::print,
            py::arg("label") = "none"
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