#include "SplitGaussianMarginal.h"

SplitGaussianMarginal::SplitGaussianMarginal(double mu, double sigma_p, double sigma_m, unsigned int seed)
    : mu(mu), sigma_p(sigma_p), sigma_m(sigma_m)
{
    N = 2.0 / (sigma_p + sigma_m);
    w = sigma_m / (sigma_p + sigma_m);

    eng_ = gsl_rng_alloc(rng_tp);
    gsl_rng_set(eng_, seed);
}

Vector SplitGaussianMarginal::rvs(std::size_t n) {
    Vector z(n);
    for (std::size_t i = 0; i < n; ++i) 
        z[i] = ppf(gsl_ran_flat(eng_, 0, 1));
    return z;
}

double SplitGaussianMarginal::logpdf(double x) {
    return std::log(N * gsl_ran_gaussian_pdf((x - mu) / (x > mu ? sigma_p : sigma_m), 1));
}

double SplitGaussianMarginal::cdf(double x) {
    if (x > mu) {
        return 2 * w * gsl_cdf_gaussian_P((x - mu) / sigma_m, 1);
    } else {
        return w + 2 * (1 - w) * (gsl_cdf_gaussian_P((x - mu) / sigma_p, 1) - 0.5);
    }
}

double SplitGaussianMarginal::ppf(double p) {
    double u, s;
    if (p > w) {
        u = (1 - 2 * w + p) / (2 * (1 - w));
        s = sigma_p;
    } else {
        u = p / (2 * w);
        s = sigma_m;
    }
    return mu + s * gsl_cdf_gaussian_Pinv(u, 1);
}

double SplitGaussianMarginal::mean() {
    return mu + std::sqrt(2 / PI) * (sigma_p - sigma_m);
}

double SplitGaussianMarginal::std() {
    return std::sqrt((1 - 2 / PI) * std::pow(sigma_p - sigma_m, 2) + sigma_p * sigma_m);
}