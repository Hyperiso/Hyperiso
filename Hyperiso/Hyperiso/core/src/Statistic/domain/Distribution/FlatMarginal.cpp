#include "FlatMarginal.h"

FlatMarginal::FlatMarginal(double a, double b, unsigned int seed)
    : a(a), b(b)
{
    eng_ = gsl_rng_alloc(rng_tp);
    gsl_rng_set(eng_, seed);
}

Vector FlatMarginal::rvs(std::size_t n) {
    Vector z(n);
    for (std::size_t i = 0; i < n; ++i) 
        z[i] = gsl_ran_flat(eng_, a, b);
    return z;
}

double FlatMarginal::logpdf(double x) {
    return std::log(gsl_ran_flat_pdf(x, a, b));
}

double FlatMarginal::cdf(double x) {
    return gsl_cdf_flat_P(x, a, b);
}

double FlatMarginal::ppf(double p) {
    return gsl_cdf_flat_Pinv(p, a, b);
}

double FlatMarginal::mean() {
    return 0.5 * (a + b);
}

double FlatMarginal::std() {
    return (b - a) / std::sqrt(12.); 
}