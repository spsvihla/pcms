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

/**
 * @brief A distribution that implements the critical beta splitting distribution.
 * 
 * This class generates random numbers following the critical beta-splitting
 * distribution, which is used for sampling critical binary trees.
 */
class critical_beta_splitting_distribution
{
public:
    /**
     * @brief Constructs the critical beta splitting distribution for a given size.
     * 
     * @param n The number of leaf nodes in the tree.
     */
    critical_beta_splitting_distribution(int n);

    /**
     * @brief Samples from the critical beta splitting distribution.
     * 
     * @tparam Generator A random number generator type.
     * @param rng The random number generator to use for sampling.
     * @return A random index sampled according to the distribution.
     */
    template<typename Generator> int
    operator()(Generator &rng)
    {
        std::uniform_real_distribution<double> unif(0.0, 1.0);
        double u = unif(rng);
        for(int i = 0; i < n-1; i++)
        {
            if(cdf[i] >= u)
            {
                return i + 1;
            }
        }
        return n-1; // this line shouldn't be reached
    }

private:
    int n; ///< The number of leaves (size of the tree)
    std::vector<double> pmf; ///< Probability mass function (PMF) values
    std::vector<double> cdf; ///< Cumulative distribution function (CDF) values

    static std::unordered_map<int, std::pair<std::vector<double>, std::vector<double>>> cache;
};

/**
 * @brief Samples a uniformly random binary tree.
 * 
 * This function generates a random binary tree using a uniform distribution.
 * 
 * @param n_leaves The number of leaf nodes in the tree.
 * @param seed The seed for the random number generator (default is 0).
 * @return A pointer to the newly created tree.
 */
Tree *remy(int n_leaves, unsigned int seed = 0);

/**
 * @brief Samples a critical beta-splitting tree.
 * 
 * This function generates a binary tree using the critical beta-splitting algorithm.
 * 
 * @param n_leaves The number of leaf nodes in the tree.
 * @param seed The seed for the random number generator (default is 0).
 * @return A pointer to the newly created tree.
 */
Tree *cbst(int n_leaves, unsigned int seed = 0);

#endif // TREE_DIST_H