#ifndef RAND_HPP
#define RAND_HPP

// standard library includes
#include <random>

// Constants
constexpr double GAMMA = 0.57721566490153286060; // Euler-Mascheroni constant   


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

// compute the nth harmonic number
inline double 
harmonic_number(int n) 
{
    double hn;
    if(n <= 30)
    {
        hn = 0.0;
        for(int k = 1; k <= n; ++k)
        {
            hn += 1.0 / k;
        }
    }
    else
    {
        // error < 10e-8 for n > 30
        hn = std::log(n) + GAMMA + 1.0 / (2 * n) - 1.0 / (12 * n * n);
    }
    return hn;
}

// sample the critical beta split distribution 
inline double 
rand_critical_beta_split(int n, std::mt19937& rng)
{
    double hn = harmonic_number(n - 1);
    double factor = n / (2.0 * hn);

    double u = rand_uniform_double(0.0, 1.0, rng);
    double thresh = u / factor;
    double cdf_val = 0.0;
    for(int i = 1; i <= n - 1; ++i)
    {
        int denom = i * (n - i);
        cdf_val += 1.0 / denom;
        if(cdf_val >= thresh)
        {
            return i;
        }
    }
    return n-1; // this line shouldn't be reached
}

#endif // RAND_HPP