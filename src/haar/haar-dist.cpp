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

// Project-specific includes
#include "haar-dist.hpp"
#include "tree-dist.hpp"

// Pybind11 includes 
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;


// Fisher-Yates shuffle
inline void
permute_indices(py::detail::unchecked_mutable_reference<int, 1L>& perm_, 
                py::ssize_t n, std::mt19937& rng)
{
    for(py::ssize_t i = n - 1; i > 0; --i)
    {
        std::uniform_int_distribution<int> dist(0, i);
        int j = dist(rng);
        int tmp = perm_[i];
        perm_(i) = perm_[j];
        perm_(j) = tmp;
    }
}

// populate out_ with Haar-like wavelet coordinates
inline void
populate_coordinates(const py::detail::unchecked_reference<double, 1L>& func_, 
                     const py::detail::unchecked_mutable_reference<int, 1L>& perm_, 
                     py::detail::unchecked_mutable_reference<double, 1L>& out_,
                     py::ssize_t n)
{
    double n_ = static_cast<double>(n);

    // compute coordinate when split n -> (1, n-1)
    double l = func_[perm_[0]];
    double r = 0.0;
    for(py::ssize_t i = 1; i < n; ++i)
    {
        r += func_[perm_[i]];
    }
    double coef = std::sqrt((n_ - 1) / n_);
    out_(0) = coef * std::abs(l - r);

    // compute coordinates when split n -> (i+1, n-(i+1))
    for(py::ssize_t i = 1; i < n - 1; ++i)
    {
        double f = func_[perm_[i]];
        l = (i * l + f) / (i + 1);
        r = ((n_ - i) * r - f) / (n_ - i - 1);
        coef = std::sqrt(i * (n_ - i) / n_);
        out_(i) = coef * std::abs(l - r);
    }
}

inline void
update_total(const py::detail::unchecked_reference<double, 1L>& ys_,
             const py::detail::unchecked_mutable_reference<double, 1L>& coords_,
             const py::detail::unchecked_reference<double, 1L>& pmf_,
             py::detail::unchecked_mutable_reference<double, 1L>& total_,
             py::ssize_t n, py::ssize_t m)
{
    for(py::ssize_t j = 0; j < m; ++j)
    {
        double val = 0.0;
        // TODO: SIMD inner loop
        for(py::ssize_t k = 0; k < n; ++k)
        {
            val += pmf_[k] * static_cast<double>(coords_[k] <= ys_(j));
        }
        total_(j) += val;
    }
}

// evaluate cdf of <f,φ> when φ ~ CBST(n)
py::array_t<double> 
cdf_rand_basis(const py::array_t<double>& ys, const py::array_t<double>& func, 
               const py::array_t<double>& pmf, int num_iter, 
               std::optional<unsigned int> seed)
{   
    // random number generator
    unsigned int seed_ = seed.value_or(std::random_device{}());
    std::mt19937 rng(seed_);

    // get sizes
    py::ssize_t n = func.shape(0);
    py::ssize_t m = ys.shape(0);

    // unchecked accessors
    auto ys_ = ys.unchecked<1>();
    auto func_ = func.unchecked<1>();
    auto pmf_ = pmf.unchecked<1>();

    // allocate output array 
    py::array_t<double> total(m);
    auto total_ = total.mutable_unchecked<1>();

    // allocate arrays for coordinates and permuation
    py::array_t<double> coords(n - 1);
    py::array_t<int> perm(n);
    auto coords_ = coords.mutable_unchecked<1>();
    auto perm_ = perm.mutable_unchecked<1>();

    // initialize output array
    for(py::ssize_t i = 0; i < m; ++i)
    {
        total_(i) = 0.0;
    }

    // initialize permuation array 
    for(py::ssize_t i = 0; i < n; ++i)
    {
        perm_(i) = static_cast<int>(i);
    }

    // estimate cdf value
    // TODO: OpenMP loop
    for(int i = 0; i < num_iter; ++i)
    {
        permute_indices(perm_, n, rng);
        populate_coordinates(func_, perm_, coords_, n);
        update_total(ys_, coords_, pmf_, total_, n, m);
    }

    // normalize
    double num_iter_ = static_cast<double>(num_iter);
    for(py::ssize_t i = 0; i < m; ++i)
    {
        total_(i) /= num_iter_;
    }

    return total;
}