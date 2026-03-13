#include "GaussianMarginal.h"

GaussianMarginal::GaussianMarginal(double mu, double sigma, unsigned int seed)
    : mu(mu), sigma(sigma)
{
    gsl_rng_set(eng_.get(), seed);
}

Vector GaussianMarginal::rvs(std::size_t n) {
    Vector z(n);
    for (std::size_t i = 0; i < n; ++i) 
        z[i] = mu + gsl_ran_gaussian(eng_.get(), sigma);
    return z;
}

double GaussianMarginal::logpdf(double x) {
    return -0.5 * (std::log(2 * PI * sigma * sigma) + std::pow((x - mu) / sigma, 2));
    // return -0.5 * std::pow((x - mu) / sigma, 2);
}

double GaussianMarginal::cdf(double x) {
    return gsl_cdf_gaussian_P(x - mu, sigma);
}

double GaussianMarginal::ppf(double p) {
    return mu + gsl_cdf_gaussian_Pinv(p, sigma);
}

double GaussianMarginal::mean() {
    return mu;
}

double GaussianMarginal::std() {
    return sigma;
}