/**
 * @file tree-dist.cpp
 * 
 * This file implements the `critical_beta_splitting_distribution` class and 
 * the functions `remy` and `cbst` for generating random trees.
 * 
 * @author Sean Svihla
 */

// Standard library includes
#include <cassert>
#include <chrono>
#include <cmath>
#include <memory>
#include <random>
#include <stack>
#include <tuple>
#include <unordered_map>

// Other system library includes
#include <gsl/gsl_integration.h>

// Project-specific includes
#include "tree-dist.hpp"
#include "tree.hpp"

// Pybind11 includes
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

// Constants
constexpr double GAMMA = 0.57721566490153286060; // Euler-Mascheroni constant   


/* Integrand to compute harmonic numbers by Euler's formula  */
double 
integrand(double x, void* params) {
    int n = *(int*)params;
    return (1 - pow(x, n)) / (1 - x);
}

/* Computes the harmonic number for a given n. */
double 
harmonic_number(int n) 
{
    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(1000);

    // define integrand
    gsl_function F;
    F.function = &integrand;
    F.params = &n;

    // numerical integration
    double result, error;
    gsl_integration_qags(
        &F,   // Function to integrate (gsl_function struct)
        0,    // Lower limit of integration
        1,    // Upper limit of integration
        1e-7, // Absolute error tolerance
        1e-7, // Relative error tolerance
        1000, // Maximum number of subintervals
        workspace, // Integration workspace
        &result,   // Pointer to store the integral result
        &error     // Pointer to store the estimated error
    );    

    gsl_integration_workspace_free(workspace);
    return result;
}

critical_beta_splitting_distribution::critical_beta_splitting_distribution(py::ssize_t n)
: n(n), pmf(n-1), cdf(n-1) 
{
    if(n < 2)
    {
        throw py::value_error("Number of leaves must be at least 2.");
    }

    double hn = harmonic_number(n - 1);

    auto pmf_ = pmf.mutable_unchecked<1>();
    auto cdf_ = cdf.mutable_unchecked<1>();

    double factor = n / (2.0 * hn);
    for(py::ssize_t i = 1; i < n; i++) 
    {
        pmf_[i - 1] = factor / static_cast<double>(i * (n - i));
    }

    cdf_(0) = pmf_[0];
    for(py::ssize_t i = 1; i < n - 2; i++) 
    {
        cdf_(i) = cdf_[i-1] + pmf_[i];
    }
    cdf_(n-2) = 1.0;
}

py::array_t<double>
critical_beta_splitting_distribution::get_pmf() const
{
    return pmf;
}

py::array_t<double>
critical_beta_splitting_distribution::get_cdf() const
{
    return cdf;
}

Tree*
remy(int n_leaves, std::optional<unsigned int> seed)
{
    unsigned int seed_ = seed.value_or(std::random_device{}());
    std::default_random_engine rng(seed_);
    std::uniform_int_distribution<int> bern(0,1); // Bernoulli(p=0.5)

    Tree* tree = new Tree(2 * n_leaves - 1);
    for(int node = 1; node < 2 * n_leaves - 1; node += 2)
    {
        std::uniform_int_distribution<int> unif(0, node-1);
        int rand_node = unif(rng);
        int b = bern(rng);

        int left_child = (node + 1) * b + rand_node * (1 - b);
        int right_child = rand_node * b + (node + 1) * (1 - b);
        int parent = tree->get_parent(rand_node);

        tree->cut(rand_node);
        tree->link(node, parent);
        tree->link(right_child, node);
        tree->link(left_child, node);
    }

    return tree;
}

void 
_cbst(Tree* tree, int start, int end, int n0, std::mt19937 &rng) 
{
    std::stack<std::tuple<int, int, int, int>> stack;
    stack.push({start, end, n0, 0});

    while (!stack.empty()) {
        auto [start, end, n0, depth] = stack.top();
        stack.pop();

        if (n0 == 1) {
            continue;
        }
        if (n0 == 2) {
            tree->link(start + 1, end);
            tree->link(start, end);
            continue;
        }

        critical_beta_splitting_distribution dist(n0);
        int x = dist(rng);
        int n1 = 2 * x - 1;

        tree->link(end - 1, end);
        tree->link(start + n1 - 1, end);

        stack.push({start, start + n1 - 1, x, depth + 1});
        stack.push({start + n1, end - 1, n0 - x, depth + 1});
    }
}

Tree*
cbst(int n_leaves, std::optional<unsigned int> seed)
{
    unsigned int seed_ = seed.value_or(std::random_device{}());
    std::mt19937 rng(seed_);

    Tree* tree = new Tree(2 * n_leaves - 1);
    _cbst(tree, 0, tree->get_n_nodes() - 1, n_leaves, rng);
    return tree;
}