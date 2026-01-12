#include "GaussianMarginal.h"

GaussianMarginal::GaussianMarginal(unsigned int seed, double mu, double sigma)
    : mu(mu), sigma(sigma)
{
    eng_ = gsl_rng_alloc(rng_tp);
    gsl_rng_set(eng_, seed);
}

Vector GaussianMarginal::rvs(std::size_t n) {
    Vector z(n);
    for (std::size_t i = 0; i < n; ++i) 
        z[i] = mu + gsl_ran_gaussian(eng_, sigma);
    return z;
}

double GaussianMarginal::logpdf(double x) {
    return std::log(gsl_ran_gaussian_pdf(x - mu, sigma));
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