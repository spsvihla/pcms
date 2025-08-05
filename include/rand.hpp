#ifndef RAND_HPP
#define RAND_HPP

// standard library includes
#include <random>

// Exponential(rate=lam)
inline double 
rand_exponential(double lam, std::mt19937& rng) 
{
    double U = (rng() + 0.5) * (1.0 / (rng.max() + 1.0));   // Uniform(0, 1)
    return -std::log(U) / lam;
}

// Dirichlet(1,...,1)
inline std::vector<double>
rand_dirichlet_uniform(int N, int k, std::mt19937& rng) 
{
    std::vector<double> samples(N * k);
    for(int i = 0; i < N; ++i) 
    {
        double sum = 0.0;
        for(int j = 0; j < k; ++j) 
        {
            double e = rand_exponential(1.0, rng);
            samples[i * k + j] = e;
            sum += e;
        }
        for(int j = 0; j < k; ++j) 
        {
            samples[i * k + j] /= sum;
        }
    }
    return samples;
}

#endif // RAND_HPP