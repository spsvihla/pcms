#ifndef RAND_HPP
#define RAND_HPP

// standard library includes
#include <random>


// Uniform(a, b)
inline double
rand_uniform_double(double a, double b, std::mt19937& rng)
{
    return a + ((rng() + 0.5) / (rng.max() + 1.0)) * (b - a);
}

// Uniform(a, b)
inline int 
rand_uniform_int(int a, int b, std::mt19937& rng)
{
    int range = b - a + 1;
    using result_t = std::mt19937::result_type;

    result_t max = rng.max();
    result_t bucket_size = (max + 1) / static_cast<result_t>(range);

    result_t r;
    do 
    {
        r = rng();
    }
    while(r >= bucket_size * static_cast<result_t>(range));

    return a + static_cast<int>(r / bucket_size);
}

// Exponential(rate=lam)
inline double 
rand_exponential(double lam, std::mt19937& rng) 
{
    double u = rand_uniform_double(0.0, 1.0, rng);
    return -std::log(u) / lam;
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