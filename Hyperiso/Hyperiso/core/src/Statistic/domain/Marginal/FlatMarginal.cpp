#include "FlatMarginal.h"

FlatMarginal::FlatMarginal(double a, double b, unsigned int seed)
    : a(a), b(b)
{
    gsl_rng_set(eng_.get(), seed);
}

std::vector<double> FlatMarginal::rvs(std::size_t n) {
    std::vector<double> z(n);
    for (std::size_t i = 0; i < n; ++i) 
        z[i] = gsl_ran_flat(eng_.get(), a, b);
    return z;
}

double FlatMarginal::logpdf(double x) {
    return (x > a && x < b) ? std::log(gsl_ran_flat_pdf(x, a, b)) : -1e100;
}

PDFDiff FlatMarginal::f_df_ddf(double x) {
    double f = logpdf(x);
    return {f, 0.0, 0.0};
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