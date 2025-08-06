/**
 * @file tree-dist.cpp
 * 
 * This file implements the `critical_beta_splitting_distribution` class and 
 * the functions `remy` and `cbst` for generating random trees.
 * 
 * @author Sean Svihla
 */

// Standard library includes
#include <algorithm>    // std::min
#include <cmath>        // std::log, std::sqrt
#include <random>       // std::mt19937, RNG
#include <stack>        // std::stack
#include <tuple>        // std::tuple, structured bindings

// Project-specific includes
#include "tree-dist.hpp"
#include "tree.hpp"
#include "thread-pool.hpp"
#include "rand.hpp"

// Pybind11 includes
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

// Constants
constexpr double GAMMA = 0.57721566490153286060; // Euler-Mascheroni constant   


// compute the nth harmonic number
double 
harmonic_number(int n) 
{
    double hn;
    if(n <= 1e6)
    {
        hn = 0.0;
        for(int k = 1; k <= n; ++k)
        {
            hn += 1.0 / k;
        }
    }
    else
    {
        hn = std::log(n) + GAMMA + 1.0 / (2 * n) - 1.0 / (12 * n * n);
    }
    return hn;
}

// sample the critical beta split distribution 
double rand_critical_beta_split(int n, std::mt19937& rng)
{
    double hn = harmonic_number(n - 1);
    double factor = n / (2.0 * hn);

    double u = rand_uniform_double(0.0, 1.0, rng);
    double cdf_val = 0.0;
    for(int i = 1; i <= n-1; ++i)
    {
        cdf_val += factor / (i * (n - i));
        if(cdf_val >= u)
        {
            return i;
        }
    }
    return n-1; // this line shouldn't be reached
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

void
remy(Tree* tree, bool planted, bool do_randomize_edge_lengths, std::optional<unsigned int> seed)
{
    unsigned int seed_ = seed.value_or(std::random_device{}());
    std::mt19937 rng(seed_);

    int n_nodes = tree->get_n_nodes();
    int end = planted ? n_nodes - 2 : n_nodes - 1;

    if(planted)
    {
        tree->link_(n_nodes - 2, n_nodes - 1);
    }

    for(int node = 1; node < end; node += 2)
    {
        int rand_node = rand_uniform_int(0, node-1, rng);
        int b = rand_uniform_int(0, 1, rng);                // Bernoulli(p=0.5)

        int left_child = (node + 1) * b + rand_node * (1 - b);
        int right_child = rand_node * b + (node + 1) * (1 - b);
        int parent = tree->get_parent(rand_node);

        tree->cut_(rand_node);
        tree->link_(node, parent);
        tree->link_(right_child, node);
        tree->link_(left_child, node);
    }

    tree->update_subtree_size();

    if(do_randomize_edge_lengths)
    {
        randomize_edge_lengths(tree, seed);
    }
}

inline void 
cbst_(Tree* tree, int start, int end, int n0, std::mt19937 &rng) 
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
            tree->link_(start + 1, end);
            tree->link_(start, end);
            continue;
        }

        int x = rand_critical_beta_split(n0, rng);
        int n1 = 2 * x - 1;

        tree->link_(end - 1, end);
        tree->link_(start + n1 - 1, end);

        stack.push({start, start + n1 - 1, x, depth + 1});
        stack.push({start + n1, end - 1, n0 - x, depth + 1});
    }
}

void
cbst(Tree* tree, bool planted, bool do_randomize_edge_lengths, std::optional<unsigned int> seed)
{
    unsigned int seed_ = seed.value_or(std::random_device{}());
    std::mt19937 rng(seed_);

    int n_nodes = tree->get_n_nodes();
    int n_leaves = planted ? n_nodes / 2 : (n_nodes + 1) / 2;
    int end = planted ? n_nodes - 2 : n_nodes - 1;

    if(planted)
    {
        tree->link_(n_nodes - 2, n_nodes - 1);
    }
    cbst_(tree, 0, end, n_leaves, rng);

    tree->update_subtree_size();

    if(do_randomize_edge_lengths)
    {
        randomize_edge_lengths(tree, seed);
    }
}

void 
batched_tree_generator(
    std::function<void(Tree*, bool, bool, std::optional<unsigned int>)> tree_builder,
    std::vector<Tree*>& buffer,
    bool planted,
    bool do_randomize_edge_lengths,
    int n_samples,
    std::optional<unsigned int> seed)
{
    // generate random seeds
    unsigned int seed_ = seed.value_or(std::random_device{}());
    std::mt19937 seed_rng(seed_);

    std::vector<unsigned int> seeds(n_samples);
    for(int i = 0; i < n_samples; ++i) 
    {
        seeds[i] = seed_rng();
    }

    // run tasks
    std::size_t n_threads = std::thread::hardware_concurrency();
    n_threads = (n_threads == 0) ? 4 : n_threads;
    n_threads = std::min(n_threads, static_cast<std::size_t>(n_samples));

    ThreadPool pool(n_threads);

    std::vector<std::future<void>> futures;
    futures.reserve(n_samples);

    for(int i = 0; i < n_samples; ++i) 
    {
        futures.push_back(pool.submit([&, i]() {
            tree_builder(buffer[i], planted, do_randomize_edge_lengths, seeds[i]);
        }));
    }

    // wait for all
    for(auto& f : futures)
    {
        f.get();
    }
}

void
cbst_batched(std::vector<Tree*>& buffer, bool planted, bool do_randomize_edge_lengths,
             int n_samples, std::optional<unsigned int> seed)
{
    batched_tree_generator(cbst, buffer, planted, do_randomize_edge_lengths, 
                           n_samples, seed);
}

void 
remy_batched(std::vector<Tree*>& buffer, bool planted, bool do_randomize_edge_lengths,
             int n_samples, std::optional<unsigned int> seed)
{
    batched_tree_generator(remy, buffer, planted, do_randomize_edge_lengths, 
                           n_samples, seed);
}
