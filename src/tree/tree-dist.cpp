/**
 * @file tree-dist.cpp
 * 
 * This file implements the `critical_beta_splitting_distribution` class and 
 * the functions `remy` and `cbst` for generating random trees.
 * 
 * @author Sean Svihla
 */

// Standard library includes
#include <algorithm>
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
#include "thread-pool.hpp"

// Pybind11 includes
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

// Constants
constexpr double GAMMA = 0.57721566490153286060; // Euler-Mascheroni constant   


// integrand to compute harmonic numbers by Euler's formula 
double 
integrand(double x, void* params) {
    int n = *(int*)params;
    return (1 - pow(x, n)) / (1 - x);
}

// compute the nth harmonic number
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

critical_beta_splitting_distribution::critical_beta_splitting_distribution(int n)
: n(n), pmf(n-1), cdf(n-1) 
{
    double hn = harmonic_number(n - 1);

    double factor = n / (2.0 * hn);
    for(py::ssize_t i = 1; i < n; i++) 
    {
        pmf[i - 1] = factor / static_cast<double>(i * (n - i));
    }

    cdf[0] = pmf[0];
    for(py::ssize_t i = 1; i < n - 2; i++) 
    {
        cdf[i] = cdf[i-1] + pmf[i];
    }
    cdf[n-2] = 1.0;
}

std::vector<double>
critical_beta_splitting_distribution::get_pmf_() const
{
    return pmf;
}

py::array_t<double>
critical_beta_splitting_distribution::get_pmf() const
{
    return py::array(pmf.size(), pmf.data());
}

std::vector<double>
critical_beta_splitting_distribution::get_cdf_() const
{
    return cdf;
}

py::array_t<double>
critical_beta_splitting_distribution::get_cdf() const
{
    return py::array(cdf.size(), cdf.data());
}

Tree*
remy(int n_leaves, bool planted, std::optional<unsigned int> seed)
{
    unsigned int seed_ = seed.value_or(std::random_device{}());
    std::default_random_engine rng(seed_);
    std::uniform_int_distribution<int> bern(0,1); // Bernoulli(p=0.5)

    int n_nodes = planted ? 2 * n_leaves : 2 * n_leaves - 1;
    int end = planted ? n_nodes - 2 : n_nodes - 1;

    Tree* tree = new Tree(n_nodes);
    if(planted)
    {
        tree->link(n_nodes - 2, n_nodes - 1);
    }
    
    for(int node = 1; node < end; node += 2)
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

// Exponential(rate=lam)
inline double 
rand_exponential(double lam, std::mt19937& rng) 
{
    double U = (rng() + 0.5) * (1.0 / (rng.max() + 1.0));   // Uniform(0, 1)
    return -std::log(U) / lam;
}

inline void
randomize_edge_lengths(Tree* tree, std::optional<unsigned int> seed)
{
    unsigned int seed_ = seed.value_or(std::random_device{}());
    std::mt19937 rng(seed_);
    for(int i = 0; i < tree->get_n_nodes(); ++i)
    {
        tree->set_edge_length(i, rand_exponential(tree->get_subtree_size(i), rng));
    }
}

inline void 
_cbst(Tree* tree, int start, int end, int n0, std::mt19937 &rng) 
{
    std::stack<std::tuple<int, int, int, int>> stack;
    stack.push({start, end, n0, 0});

    while(!stack.empty()) 
    {
        auto [start, end, n0, depth] = stack.top();
        stack.pop();

        if(n0 == 1) 
        {
            continue;
        }
        if(n0 == 2) 
        {
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
cbst(int n_leaves, bool planted, bool do_randomize_edge_lengths, std::optional<unsigned int> seed)
{
    unsigned int seed_ = seed.value_or(std::random_device{}());
    std::mt19937 rng(seed_);

    int n_nodes = planted ? 2 * n_leaves : 2 * n_leaves - 1;
    int end = planted ? n_nodes - 2 : n_nodes - 1;

    Tree* tree = new Tree(n_nodes);
    if(planted)
    {
        tree->link(n_nodes - 2, n_nodes - 1);
    }
    _cbst(tree, 0, end, n_leaves, rng);

    if(do_randomize_edge_lengths)
    {
        randomize_edge_lengths(tree, seed);
    }

    return tree;
}

py::tuple
cbst_batched(int n_leaves, bool planted, bool do_randomize_edge_lengths, int num_samples, std::optional<unsigned int> seed)
{
    // generate random seeds
    unsigned int seed_ = seed.value_or(std::random_device{}());
    std::mt19937 seed_rng(seed_);

    std::vector<unsigned int> seeds(num_samples);
    for(int i = 0; i < num_samples; ++i)
    {
        seeds[i] = seed_rng();
    }

    // run tasks
    std::size_t n_threads = std::thread::hardware_concurrency();
    n_threads = (n_threads == 0) ? 4 : n_threads;
    n_threads = std::min(n_threads, static_cast<std::size_t>(num_samples));

    ThreadPool pool(n_threads);

    std::vector<std::future<Tree*>> futures;
    futures.reserve(num_samples);

    for(int i = 0; i < num_samples; ++i) 
    {
        futures.push_back(pool.submit([=]() {
            return cbst(n_leaves, planted, do_randomize_edge_lengths, seeds[i]);
        }));
    }

    // collect results
    py::tuple result(num_samples);
    for(int i = 0; i < num_samples; ++i) 
    {
        Tree* tree = futures[i].get();
        result[i] = py::cast(tree);
    }

    return result;
}
