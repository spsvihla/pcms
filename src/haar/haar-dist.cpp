/**
 * @file haar-dist.cpp
 * @author Sean Svihla
 */

// Standard library includes
#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <vector>

// Other includes
#include <immintrin.h>

// Project-specific includes
#include "haar-dist.hpp"
#include "tree-dist.hpp"

// Pybind11 includes 
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;


// Exponential(1.0)
inline double 
rand_exponential(std::mt19937& rng) {
    double U = (rng() + 0.5) * (1.0 / (rng.max() + 1.0));   // Uniform(0, 1)
    return -std::log(U);
}

// Dirichlet(1,...,1)
inline std::vector<double>
rand_dirichlet_uniform(int N, int k, std::mt19937& rng) {
    std::vector<double> samples(N * k);
    for (int i = 0; i < N; ++i) {
        double sum = 0.0;
        for (int j = 0; j < k; ++j) {
            double e = rand_exponential(rng);
            samples[i * k + j] = e;
            sum += e;
        }
        for (int j = 0; j < k; ++j) {
            samples[i * k + j] /= sum;
        }
    }
    return samples;
}

// stack copies of an array
std::vector<double>
make_copies(double* arr, int size, int n_copies)
{
    std::vector<double> out(size * n_copies);
    #pragma omp parallel for
    for(std::size_t i = 0; i < static_cast<std::size_t>(n_copies); ++i) 
    {
        std::copy(arr, arr + size, out.begin() + i * size);
    }
    return out;
}

// in-place Fisher-Yates shuffle
inline void
fisher_yates(double* arr, int size, std::mt19937& rng)
{
    std::size_t max = rng.max();
    for(std::size_t j = size - 1; j + 1 >= 4; j -= 4) 
    {
        // iteration 1
        std::size_t k1;
        do 
        {
            k1 = rng();
        } while (k1 >= max - (max % (j + 1)));
        k1 %= (j + 1);
        std::swap(arr[j], arr[k1]);

        // iteration 2
        std::size_t k2;
        do 
        {
            k2 = rng();
        } while (k2 >= max - (max % (j)));
        k2 %= j;
        std::swap(arr[j - 1], arr[k2]);

        // iteration 3
        std::size_t k3;
        do 
        {
            k3 = rng();
        } while (k3 >= max - (max % (j - 1)));
        k3 %= (j - 1);
        std::swap(arr[j - 2], arr[k3]);

        // iteration 4
        std::size_t k4;
        do 
        {
            k4 = rng();
        } while (k4 >= max - (max % (j - 2)));
        k4 %= (j - 2);
        std::swap(arr[j - 3], arr[k4]);
    }

    // Remainder loop for leftovers
    for (std::size_t j = (size - 1) % 4; j > 0; --j) {
        std::size_t k;
        do {
            k = rng();
        } while (k >= max - (max % (j + 1)));
        k %= (j + 1);
        std::swap(arr[j], arr[k]);
    }

}

inline void
compute_wavelet(Tree* tree, double* out)
{
    std::vector<int> children = tree->find_children_(tree->get_n_nodes() - 2);

    int lsize = tree->get_subtree_size(children[0]);
    int rsize = tree->get_subtree_size(children[1]);
    int sum = lsize + rsize;

    double rsize_ = static_cast<double>(rsize);
    double lsize_ = static_cast<double>(lsize);
    double sum_ = static_cast<double>(sum);

    double lval =  sqrt(rsize_ / (lsize_ * sum_));
    double rval = -sqrt(lsize_ / (rsize_ * sum_));

    // left subtree
    #pragma omp simd
    for(int j = 0; j < lsize; ++j)
    {
        out[j] = lval;
    }

    // right subtree
    #pragma omp simd
    for(int j = lsize; j < sum; ++j)
    {
        out[j] = rval;
    }
}

inline void
compute_trace_length(Tree* tree, double* out)
{
    int n_leaves = tree->find_n_leaves();

    std::vector<int> subtree_starts = tree->find_subtree_start_indices_();

    #pragma omp simd
    for(int i = 0; i < n_leaves; ++i)
    {
        out[i] = 0.0;
    }

    for(int i = 0; i < tree->get_n_nodes() - 2; ++i)
    {
        int start = subtree_starts[i];
        double tbl = tree->find_tbl(i);

        #pragma omp simd
        for(int j = 0; j < tree->get_subtree_size(i); ++j)
        {
            out[start + j] += tbl;
        }
    }
}

py::array_t<double>
sample_dh_component(const py::array_t<double>& f, int n_samples, 
                    std::optional<unsigned int> seed)
{
    unsigned int seed_ = seed.value_or(std::random_device{}());
    std::mt19937 rng(seed_);

    int n_leaves = static_cast<int>(f.size());
    double* f_ = static_cast<double*>(f.request().ptr);

    // permute f
    std::vector<double> perms = make_copies(f_, n_leaves, n_samples);
    #pragma omp parallel for
    for(int i = 0; i < n_samples; ++i)
    {
        fisher_yates(perms.data() + i * n_leaves, n_leaves, rng);
    }

    // sample trees
    std::vector<Tree*> trees = cbst_batched(n_leaves, true, true, n_samples, seed_);

    std::vector<double> wavelets(n_leaves * n_samples);
    #pragma omp parallel for
    for(int i = 0; i < n_samples; ++i)
    {
        compute_wavelet(trees[i], wavelets.data() + i * n_leaves);
    }

    std::vector<double> trace_lengths(n_leaves * n_samples);
    #pragma omp parallel for
    for(int i = 0; i < n_samples; ++i)
    {
        compute_trace_length(trees[i], trace_lengths.data() + i * n_leaves);
    }

    // clean up trees
    for(Tree* t : trees) 
    {
        delete t;
    }

    // perform element-wise operations
    #pragma omp simd
    for(int i = 0; i < n_leaves * n_samples; ++i)
    {
        perms[i] *= wavelets[i];
        wavelets[i] *= wavelets[i];
        trace_lengths[i] *= wavelets[i]; 
    }

    // perform row-wise operations
    std::vector<double> samples(n_samples);
    #pragma omp parallel for
    for(int i = 0; i < n_samples; ++i)
    {
        double sum1 = 0.0;
        double sum2 = 0.0;
        #pragma omp simd reduction(+:sum1,sum2)
        for(int j = 0; j < n_leaves; ++j)
        {
            sum1 += perms[i * n_leaves + j];
            sum2 += trace_lengths[i * n_leaves + j];
        }
        samples[i] = sum1 * sum1 * sum2;
    }

    return py::array(samples.size(), samples.data());
}