/**
 * @file tree-dist.hpp
 * @brief Defines the critical_beta_splitting_distribution class and utility 
 *        functions for sampling critical beta-splitting and uniformly random binary trees.
 * 
 * This file contains the class definition for `critical_beta_splitting_distribution`, 
 * which implements the critical beta-splitting distribution. The file also includes 
 * functions for sampling uniformly random binary trees (`remy`) and critical beta-
 * splitting trees (`cbst`).
 * 
 * The `critical_beta_splitting_distribution` class supports sampling random trees 
 * by using a probability distribution based on the harmonic number. The `remy` and 
 * `cbst` functions generate random trees based on uniform and critical beta-splitting 
 * distributions, respectively.
 * 
 * @section references References
 * - Aldous, David. 1996. “Probability Distributions on Cladograms.” In Random Discrete Structures, edited by David Aldous and Robin Pemantle, 76:1–18. The IMA Volumes in Mathematics and Its Applications. New York, NY: Springer New York. https://doi.org/10.1007/978-1-4612-0719-1_1.
 * - Aldous, David, and Boris Pittel. 2023. “The Critical Beta-Splitting Random Tree I: Heights and Related Results.” arXiv. https://doi.org/10.48550/arXiv.2302.05066.
 * - Aldous, David J., and Svante Janson. 2024. “The Critical Beta-Splitting Random Tree II: Overview and Open Problems.” arXiv. https://doi.org/10.48550/arXiv.2303.02529.
 * - Aldous, David J., and Svante Janson. 2024. “The Critical Beta-splitting Random Tree III: The exchangeable partition representation and the fringe tree.” arXiv. http://arxiv.org/abs/2412.09655.
 * - Aldous, David J., and Svante Janson. 2024. “The Critical Beta-splitting Random Tree IV: Mellin analysis of Leaf Height.” arXiv. http://arxiv.org/abs/2412.12319.
 * 
 * @author Sean Svihla
 */
#ifndef TREE_DIST_H
#define TREE_DIST_H

// Standard library includes
#include <random>
#include <unordered_map>

// Project-specific includes
#include "tree.hpp"

// Pybind11 includes
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;


/**
 * @brief Samples a uniformly random binary tree.
 * 
 * This function generates a random binary tree using a uniform distribution.
 * 
 * @param tree A Tree object to use as a buffer.
 * @param planted Whether the tree should be planted.
 * @param do_randomize_edge_lengths Whether to randomize edge lengths.
 * @param seed The seed for the random number generator (default is 0).
 */
void remy(Tree* tree, bool planted, bool do_randomize_edge_lengths, std::optional<unsigned int> seed);

/**
 * @brief Samples a uniformly random binary tree.
 * 
 * This function returns 'n_samples' independent samples from the uniform
 * distribution over the specified binary trees.
 * 
 * @param trees A vector of Tree pointers to use as a buffer.
 * @param planted Whether the tree should be planted.
 * @param do_randomize_edge_lengths Whether to randomize edge lengths.
 * @param n_samples The number of samples to generate.
 * @param seed The seed for the random number generator (default is 0).
 */
void remy_batched(std::vector<Tree*> trees, bool planted, bool do_randomize_edge_lengths, int n_samples, std::optional<unsigned int> seed);

/**
 * @brief Samples a critical beta-splitting tree.
 * 
 * This function generates a binary tree using the critical beta-splitting algorithm.
 * 
 * @param tree A Tree object to use as a buffer.
 * @param planted Whether the tree should be planted.
 * @param do_randomize_edge_lengths Whether to randomize edge lengths.
 * @param seed The seed for the random number generator (default is 0).
 */
void cbst(Tree* tree, bool planted, bool do_randomize_edge_lengths, std::optional<unsigned int> seed);

/**
 * @brief Batched version of cbst.
 * 
 * This function returns 'n_samples' independent samples from the specified
 * critical beta-splitting tree distribution.
 * 
 * @param trees A vector of Tree pointers to use as a buffer.
 * @param planted Whether the tree should be planted.
 * @param do_randomize_edge_lengths Whether to randomize edge lengths.
 * @param n_samples The number of samples to generate.
 * @param seed The seed for the random number generator (default is 0).
 */
void cbst_batched(std::vector<Tree*> trees, bool planted, bool do_randomize_edge_lengths, int n_samples, std::optional<unsigned int> seed);

#endif // TREE_DIST_H