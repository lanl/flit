//
// © 2024. Triad National Security, LLC. All rights reserved.
//
// This program was produced under U.S. Government contract 89233218CNA000001 
// for Los Alamos National Laboratory (LANL), which is operated by 
// Triad National Security, LLC for the U.S. Department of Energy/National Nuclear 
// Security Administration. All rights in the program are reserved by 
// Triad National Security, LLC, and the U.S. Department of Energy/National 
// Nuclear Security Administration. The Government is granted for itself and 
// others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide 
// license in this material to reproduce, prepare. derivative works, 
// distribute copies to the public, perform publicly and display publicly, 
// and to permit others to do so.
//
// Author:
//    Kai Gao, kaigao@lanl.gov
//

#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <cmath>
#include <fstream>

#include "pcg_random.hpp"

using namespace std;

extern "C" {
    void uniform_rand_int(int *w, int nw, int lb, int ub);
    void uniform_rand_seed_int(int *w, int nw, int lb, int ub, int seed);
    void normal_rand_int(int *w, int nw, float mu, float sigma);
    void normal_rand_seed_int(int *w, int nw, float mu, float sigma, int seed);

    void uniform_rand(float *w, int nw, float lb, float ub);
    void uniform_rand_seed(float *w, int nw, float lb, float ub, int seed);
    void normal_rand(float *w, int nw, float mu, float sigma);
    void normal_rand_seed(float *w, int nw, float mu, float sigma, int seed);
    void cauchy_rand(float *w, int nw, float a, float b);
    void cauchy_rand_seed(float *w, int nw, float a, float b, int seed);
    void poisson_rand(float *w, int nw, float mean);
    void poisson_rand_seed(float *w, int nw, float mean, int seed);
    void exponential_rand(float *w, int nw, float lambda);
    void exponential_rand_seed(float *w, int nw, float lambda, int seed);

    void uniform_rand_double(double *w, int nw, double lb, double ub);
    void uniform_rand_seed_double(double *w, int nw, double lb, double ub, int seed);
    void normal_rand_double(double *w, int nw, double mu, double sigma);
    void normal_rand_seed_double(double *w, int nw, double mu, double sigma, int seed);
    void cauchy_rand_double(double *w, int nw, double a, double b);
    void cauchy_rand_seed_double(double *w, int nw, double a, double b, int seed);
    void poisson_rand_double(double *w, int nw, double mean);
    void poisson_rand_seed_double(double *w, int nw, double mean, int seed);
    void exponential_rand_double(double *w, int nw, double lambda);
    void exponential_rand_seed_double(double *w, int nw, double lambda, int seed);
}

// uniform distribution (device random)
void uniform_rand_int(int *w, int nw, int lb, int ub)
{

    // Seed with a real random value, if available
    pcg_extras::seed_seq_from<std::random_device> seed_source;

    // Start a random number generator
    pcg32 rng(seed_source);

    // Initialize the distribution
    std::uniform_int_distribution<int> dist(lb, ub);

    // Generate random numbers
    for (size_t n = 0; n < nw; ++n)
    {
        w[n] = dist(rng);
    }
}

// uniform distribution (device random)
void uniform_rand(float *w, int nw, float lb, float ub)
{

    // Seed with a real random value, if available
    pcg_extras::seed_seq_from<std::random_device> seed_source;

    // Start a random number generator
    pcg32 rng(seed_source);

    // Initialize the distribution
    std::uniform_real_distribution<float> dist(lb, ub);

    // Generate random numbers
    for (size_t n = 0; n < nw; ++n)
    {
        w[n] = dist(rng);
    }
}

// uniform distribution (with given random seed)
void uniform_rand_seed_int(int *w, int nw, int lb, int ub, int seed)
{

    // Start a random number generator
    pcg32 rng(seed);

    // Initialize the distribution
    std::uniform_int_distribution<int> dist(lb, ub);

    // Generate random numbers
    for (size_t n = 0; n < nw; ++n)
    {
        w[n] = dist(rng);
    }
}

// uniform distribution (with given random seed)
void uniform_rand_seed(float *w, int nw, float lb, float ub, int seed)
{

    // Start a random number generator
    pcg32 rng(seed);

    // Initialize the distribution
    std::uniform_real_distribution<float> dist(lb, ub);

    // Generate random numbers
    for (size_t n = 0; n < nw; ++n)
    {
        w[n] = dist(rng);
    }
}

// normal distribution (device random)
void normal_rand_int(int *w, int nw, float mu, float sigma)
{

    // Seed with a real random value, if available
    pcg_extras::seed_seq_from<std::random_device> seed_source;

    // Start a random number generator
    pcg32 rng(seed_source);

    // Initialize the distribution
    std::normal_distribution<float> dist(mu, sigma);

    // Generate random numbers
    for (size_t n = 0; n < nw; ++n)
    {
        w[n] = int(dist(rng));
    }
}

// normal distribution (device random)
void normal_rand(float *w, int nw, float mu, float sigma)
{

    // Seed with a real random value, if available
    pcg_extras::seed_seq_from<std::random_device> seed_source;

    // Start a random number generator
    pcg32 rng(seed_source);

    // Initialize the distribution
    std::normal_distribution<float> dist(mu, sigma);

    // Generate random numbers
    for (size_t n = 0; n < nw; ++n)
    {
        w[n] = dist(rng);
    }
}

// normal distribution (device random)
void normal_rand_seed_int(int *w, int nw, float mu, float sigma, int seed)
{

    // Start a random number generator
    pcg32 rng(seed);

    // Initialize the distribution
    std::normal_distribution<float> dist(mu, sigma);

    // Generate random numbers
    for (size_t n = 0; n < nw; ++n)
    {
        w[n] = int(dist(rng));
    }
}

// normal distribution (device random)
void normal_rand_seed(float *w, int nw, float mu, float sigma, int seed)
{

    // Start a random number generator
    pcg32 rng(seed);

    // Initialize the distribution
    std::normal_distribution<float> dist(mu, sigma);

    // Generate random numbers
    for (size_t n = 0; n < nw; ++n)
    {
        w[n] = dist(rng);
    }
}

// Cauchy distribution (device random)
void cauchy_rand(float *w, int nw, float a, float b)
{

    // Seed with a real random value, if available
    pcg_extras::seed_seq_from<std::random_device> seed_source;

    // Start a random number generator
    pcg32 rng(seed_source);

    // Initialize the distribution
    std::cauchy_distribution<float> dist(a, b);

    // Generate random numbers
    for (size_t n = 0; n < nw; ++n)
    {
        w[n] = dist(rng);
    }
}

// Cauchy distribution (device random)
void cauchy_rand_seed(float *w, int nw, float a, float b, int seed)
{

    // Start a random number generator
    pcg32 rng(seed);

    // Initialize the distribution
    std::cauchy_distribution<float> dist(a, b);

    // Generate random numbers
    for (size_t n = 0; n < nw; ++n)
    {
        w[n] = dist(rng);
    }
}

// Poisson distribution (device random)
void poisson_rand(float *w, int nw, float mean)
{

    // Seed with a real random value, if available
    pcg_extras::seed_seq_from<std::random_device> seed_source;

    // Start a random number generator
    pcg32 rng(seed_source);

    // Initialize the distribution
    std::poisson_distribution<int> dist(mean);

    // Generate random numbers
    for (size_t n = 0; n < nw; ++n)
    {
        w[n] = float(dist(rng));
    }
}

// Poisson distribution (device random)
void poisson_rand_seed(float *w, int nw, float mean, int seed)
{

    // Start a random number generator
    pcg32 rng(seed);

    // Initialize the distribution
    std::poisson_distribution<int> dist(mean);

    // Generate random numbers
    for (size_t n = 0; n < nw; ++n)
    {
        w[n] = float(dist(rng));
    }
}

// Exponential distribution (device random)
void exponential_rand(float *w, int nw, float lambda)
{

    // Seed with a real random value, if available
    pcg_extras::seed_seq_from<std::random_device> seed_source;

    // Start a random number generator
    pcg32 rng(seed_source);

    // Initialize the distribution
    std::exponential_distribution<float> dist(lambda);

    // Generate random numbers
    for (size_t n = 0; n < nw; ++n)
    {
        w[n] = float(dist(rng));
    }
}

// Exponential distribution (device random)
void exponential_rand_seed(float *w, int nw, float lambda, int seed)
{

    // Start a random number generator
    pcg32 rng(seed);

    // Initialize the distribution
    std::exponential_distribution<float> dist(lambda);

    // Generate random numbers
    for (size_t n = 0; n < nw; ++n)
    {
        w[n] = float(dist(rng));
    }
}


// uniform distribution (device random)
void uniform_rand_double(double *w, int nw, double lb, double ub)
{

    // Seed with a real random value, if available
    pcg_extras::seed_seq_from<std::random_device> seed_source;

    // Start a random number generator
    pcg32 rng(seed_source);

    // Initialize the distribution
    std::uniform_real_distribution<double> dist(lb, ub);

    // Generate random numbers
    for (size_t n = 0; n < nw; ++n)
    {
        w[n] = dist(rng);
    }
}

// uniform distribution (with given random seed)
void uniform_rand_seed_double(double *w, int nw, double lb, double ub, int seed)
{

    // Start a random number generator
    pcg32 rng(seed);

    // Initialize the distribution
    std::uniform_real_distribution<double> dist(lb, ub);

    // Generate random numbers
    for (size_t n = 0; n < nw; ++n)
    {
        w[n] = dist(rng);
    }
}

// normal distribution (device random)
void normal_rand_double(double *w, int nw, double mu, double sigma)
{

    // Seed with a real random value, if available
    pcg_extras::seed_seq_from<std::random_device> seed_source;

    // Start a random number generator
    pcg32 rng(seed_source);

    // Initialize the distribution
    std::normal_distribution<double> dist(mu, sigma);

    // Generate random numbers
    for (size_t n = 0; n < nw; ++n)
    {
        w[n] = dist(rng);
    }
}

// normal distribution (device random)
void normal_rand_seed_double(double *w, int nw, double mu, double sigma, int seed)
{

    // Start a random number generator
    pcg32 rng(seed);

    // Initialize the distribution
    std::normal_distribution<double> dist(mu, sigma);

    // Generate random numbers
    for (size_t n = 0; n < nw; ++n)
    {
        w[n] = dist(rng);
    }
}

// Cauchy distribution (device random)
void cauchy_rand_double(double *w, int nw, double a, double b)
{

    // Seed with a real random value, if available
    pcg_extras::seed_seq_from<std::random_device> seed_source;

    // Start a random number generator
    pcg32 rng(seed_source);

    // Initialize the distribution
    std::cauchy_distribution<double> dist(a, b);

    // Generate random numbers
    for (size_t n = 0; n < nw; ++n)
    {
        w[n] = dist(rng);
    }
}

// Cauchy distribution (device random)
void cauchy_rand_seed_double(double *w, int nw, double a, double b, int seed)
{

    // Start a random number generator
    pcg32 rng(seed);

    // Initialize the distribution
    std::cauchy_distribution<double> dist(a, b);

    // Generate random numbers
    for (size_t n = 0; n < nw; ++n)
    {
        w[n] = dist(rng);
    }
}

// Poisson distribution (device random)
void poisson_rand_double(double *w, int nw, double mean)
{

    // Seed with a real random value, if available
    pcg_extras::seed_seq_from<std::random_device> seed_source;

    // Start a random number generator
    pcg32 rng(seed_source);

    // Initialize the distribution
    std::poisson_distribution<int> dist(mean);

    // Generate random numbers
    for (size_t n = 0; n < nw; ++n)
    {
        w[n] = double(dist(rng));
    }
}

// Poisson distribution (device random)
void poisson_rand_seed_double(double *w, int nw, double mean, int seed)
{

    // Start a random number generator
    pcg32 rng(seed);

    // Initialize the distribution
    std::poisson_distribution<int> dist(mean);

    // Generate random numbers
    for (size_t n = 0; n < nw; ++n)
    {
        w[n] = double(dist(rng));
    }
}

// Exponential distribution (device random)
void exponential_rand_double(double *w, int nw, double lambda)
{

    // Seed with a real random value, if available
    pcg_extras::seed_seq_from<std::random_device> seed_source;

    // Start a random number generator
    pcg32 rng(seed_source);

    // Initialize the distribution
    std::exponential_distribution<double> dist(lambda);

    // Generate random numbers
    for (size_t n = 0; n < nw; ++n)
    {
        w[n] = double(dist(rng));
    }
}

// Exponential distribution (device random)
void exponential_rand_seed_double(double *w, int nw, double lambda, int seed)
{

    // Start a random number generator
    pcg32 rng(seed);

    // Initialize the distribution
    std::exponential_distribution<double> dist(lambda);

    // Generate random numbers
    for (size_t n = 0; n < nw; ++n)
    {
        w[n] = double(dist(rng));
    }
}
