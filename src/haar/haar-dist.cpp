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

// Fisher-Yates shuffle
inline void
fisher_yates(double* arr_, std::size_t n, std::mt19937& rng)
{
    std::size_t max = rng.max();
    for (std::size_t j = n - 1; j + 1 >= 4; j -= 4) 
    {
        // iteration 1
        std::size_t k1;
        do 
        {
            k1 = rng();
        } while (k1 >= max - (max % (j + 1)));
        k1 %= (j + 1);
        std::swap(arr_[j], arr_[k1]);

        // iteration 2
        std::size_t k2;
        do 
        {
            k2 = rng();
        } while (k2 >= max - (max % (j)));
        k2 %= j;
        std::swap(arr_[j - 1], arr_[k2]);

        // iteration 3
        std::size_t k3;
        do 
        {
            k3 = rng();
        } while (k3 >= max - (max % (j - 1)));
        k3 %= (j - 1);
        std::swap(arr_[j - 2], arr_[k3]);

        // iteration 4
        std::size_t k4;
        do 
        {
            k4 = rng();
        } while (k4 >= max - (max % (j - 2)));
        k4 %= (j - 2);
        std::swap(arr_[j - 3], arr_[k4]);
    }

    // Remainder loop for leftovers
    for (std::size_t j = (n - 1) % 4; j > 0; --j) {
        std::size_t k;
        do {
            k = rng();
        } while (k >= max - (max % (j + 1)));
        k %= (j + 1);
        std::swap(arr_[j], arr_[k]);
    }

}

// project f_ onto the Haar-like basis and populate out_ with its coordinates
inline void
project_haar(const double* f_, std::vector<double>& out_, std::size_t n)
{
    double n_ = static_cast<double>(n);

    // helper variables for SIMD vectorized pre-compute loops
    __m256d n_vec = _mm256_set1_pd(n_);
    __m256d one_vec = _mm256_set1_pd(1.0);
    std::size_t j;

    // pre-compute reciprocals to avoid expensive divisions
    std::vector<double> recip_jp1(n);       // reciprocals of (j+1)
    std::vector<double> recip_nmj_1(n);     // reciprocals of (n - j - 1)
    for(j = 1; j + 3 < n - 1; j += 4) 
    {
        __m256d j_vec = _mm256_set_pd(
            static_cast<double>(j),
            static_cast<double>(j + 1),
            static_cast<double>(j + 2),
            static_cast<double>(j + 3)
        );

        // 1 / (j+1)
        __m256d denom_jp1 = _mm256_add_pd(j_vec, one_vec);
        __m256d recip_jp1_vec = _mm256_div_pd(one_vec, denom_jp1);
        _mm256_storeu_pd(&recip_jp1[j], recip_jp1_vec);

        // 1 / (n - j - 1)
        __m256d denom_nmj_1 = _mm256_sub_pd(_mm256_sub_pd(n_vec, j_vec), one_vec);
        __m256d recip_nmj_1_vec = _mm256_div_pd(one_vec, denom_nmj_1);
        _mm256_storeu_pd(&recip_nmj_1[j], recip_nmj_1_vec);
    }
    // scalar fallback for remaining elements
    for (; j < n - 1; ++j) 
    {
        recip_jp1[j] = 1.0 / (j + 1);
        recip_nmj_1[j] = 1.0 / (n_ - j - 1);
    }

    // pre-compute coefficients
    std::vector<double> coefs(n);
    for(j = 1; j + 3 < n - 1; j += 4) 
    {
        __m256d j_vec = _mm256_set_pd(
            static_cast<double>(j),
            static_cast<double>(j + 1),
            static_cast<double>(j + 2),
            static_cast<double>(j + 3)
        );

        __m256d n_minus_j = _mm256_sub_pd(n_vec, j_vec);    // n - j
        __m256d prod = _mm256_mul_pd(j_vec, n_minus_j);     // j * (n - j)
        __m256d frac = _mm256_div_pd(prod, n_vec);          // j * (n - j) / n
        __m256d sqrt_vec = _mm256_sqrt_pd(frac);            // sqrt(j * (n - j) / n)

        _mm256_storeu_pd(&coefs[j], sqrt_vec);
    }
    // scalar fallback for remaining elements
    for (; j < n - 1; ++j) 
    {
        coefs[j] = std::sqrt(j * (n_ - j) / n_);
    }

    // pre-compute prefix and suffix sums
    std::vector<double> prefix_sum(n + 1, 0.0);
    for(std::size_t j = 0; j < n; ++j) 
    {
        prefix_sum[j + 1] = prefix_sum[j] + f_[j];
    }

    std::vector<double> suffix_sum(n + 1, 0.0);
    for(std::size_t j = n; j-- > 0;) 
    {
        suffix_sum[j] = suffix_sum[j + 1] + f_[j];
    }

    // compute Haar-like projections
    for(j = 1; j + 3 < n - 1; j += 4) {
        // left averages
        __m256d prefix_vals = _mm256_set_pd(
            prefix_sum[j + 4], 
            prefix_sum[j + 3], 
            prefix_sum[j + 2], 
            prefix_sum[j + 1]
        );
        __m256d idxs = _mm256_set_pd(
            static_cast<double>(j + 3 + 1),
            static_cast<double>(j + 2 + 1),
            static_cast<double>(j + 1 + 1),
            static_cast<double>(j + 0 + 1)
        );
        __m256d left_avg = _mm256_div_pd(prefix_vals, idxs);

        // right averages
        __m256d suffix_vals = _mm256_set_pd(
            suffix_sum[j + 4], 
            suffix_sum[j + 3], 
            suffix_sum[j + 2], 
            suffix_sum[j + 1]
        );
        __m256d right_denoms = _mm256_set_pd(
            static_cast<double>(n - (j + 3) - 1),
            static_cast<double>(n - (j + 2) - 1),
            static_cast<double>(n - (j + 1) - 1),
            static_cast<double>(n - (j + 0) - 1)
        );
        __m256d right_avg = _mm256_div_pd(suffix_vals, right_denoms);

        // coef * (left_avg - right_avg)
        __m256d coefs_vec = _mm256_loadu_pd(&coefs[j]);
        __m256d diff = _mm256_sub_pd(left_avg, right_avg);
        __m256d result = _mm256_mul_pd(coefs_vec, diff);
        _mm256_storeu_pd(&out_[j], result);
    }
    // scalar fallback for remaining elements
    for(; j < n - 1; ++j) 
    {
        double left_avg = prefix_sum[j + 1] / (j + 1);
        double right_avg = suffix_sum[j + 1] / (n - j - 1);
        out_[j] = coefs[j] * (left_avg - right_avg);
    }
}

// evaluate cdf of <f,φ> when φ ~ CBST(n)
py::array_t<double> 
cdf_proj_cbst(const py::array_t<double>& ys, const py::array_t<double>& f, 
              const py::array_t<double>& pmf, int num_iter, 
              std::optional<unsigned int> seed)
{   
    unsigned int seed_ = seed.value_or(std::random_device{}());
    std::mt19937 rng(seed_);

    std::size_t n = static_cast<std::size_t>(f.shape(0));
    std::size_t m = static_cast<std::size_t>(ys.shape(0));

    auto ys_  = static_cast<double*>(ys.request().ptr);
    auto f_   = static_cast<double*>(f.request().ptr);
    auto pmf_ = static_cast<double*>(pmf.request().ptr);

    // initialize output array 
    py::array_t<double> total(m);
    auto total_ = static_cast<double*>(total.request().ptr);
    for(std::size_t i = 0; i < m; ++i)
    {
        total_[i] = 0.0;
    }

    // estimate cdf value
    std::vector<double> tmp(n - 1); // temporary buffer for Haar-like coordinates
    for(int i = 0; i < num_iter; ++i)
    {
        fisher_yates(f_, n, rng);
        project_haar(f_, tmp, n);
        
        // update total
        for(std::size_t j = 0; j < m; ++j)
        {
            double val = 0.0;
            for(std::size_t k = 0; k < n - 1; ++k)
            {
                val += pmf_[k] * static_cast<double>(tmp[k] <= ys_[j]);
            }
            total_[j] += val;
        }
    }

    // normalize in place
    for(std::size_t i = 0; i < m; ++i)
    {
        total_[i] /= static_cast<double>(num_iter);
    }

    return total;
}