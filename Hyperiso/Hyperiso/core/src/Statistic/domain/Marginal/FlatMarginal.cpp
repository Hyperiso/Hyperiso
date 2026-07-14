#include "FlatMarginal.h"

FlatMarginal::FlatMarginal(double a, double b, unsigned int seed)
    : a(a), b(b)
{
    if (!std::isfinite(a) || !std::isfinite(b) || !(a < b)) {
        throw std::invalid_argument("FlatMarginal requires finite bounds with a < b");
    }
    if (!eng_) {
        throw std::runtime_error("FlatMarginal could not allocate its GSL RNG");
    }
    gsl_rng_set(eng_.get(), seed);
}

std::vector<double> FlatMarginal::rvs(std::size_t n) {
    std::vector<double> z(n);
    for (std::size_t i = 0; i < n; ++i) 
        z[i] = gsl_ran_flat(eng_.get(), a, b);
    return z;
}

double FlatMarginal::logpdf(double x) {
    return (x >= a && x <= b) ? std::log(gsl_ran_flat_pdf(x, a, b)) : -std::numeric_limits<double>::infinity();
}

PDFDiff FlatMarginal::f_df_ddf(double x) {
    const double pdf = (x >= a && x <= b) ? gsl_ran_flat_pdf(x, a, b) : 0.0;
    return {pdf, 0.0, 0.0};
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