#include "GaussianMarginal.h"

GaussianMarginal::GaussianMarginal(double mu, double sigma, unsigned int seed)
    : mu(mu), sigma(sigma)
{
    gsl_rng_set(eng_.get(), seed);
}

std::vector<double> GaussianMarginal::rvs(std::size_t n) {
    std::vector<double> z(n);
    for (std::size_t i = 0; i < n; ++i) 
        z[i] = mu + gsl_ran_gaussian(eng_.get(), sigma);
    return z;
}

//TODO :: Niels pk on a changé la ?
double GaussianMarginal::logpdf(double x) {
    return -0.5 * (std::log(2 * PI * sigma * sigma) + std::pow((x - mu) / sigma, 2));
    // return -0.5 * std::pow((x - mu) / sigma, 2);
}

// PDFDiff GaussianMarginal::f_df_ddf(double x) {
//     double s2 = sigma * sigma;
//     double f = std::exp(-std::pow((x - mu) / sigma, 2)) / std::sqrt(2 * PI) / sigma;
//     double df = -(x - mu) / s2 * f;
//     double ddf = ((std::pow(x - mu / sigma, 2)) - 1) * f / s2;

//     return {f, df, ddf};
// }

PDFDiff GaussianMarginal::f_df_ddf(double x) {
    const double s2 = sigma * sigma;
    const double z = (x - mu) / sigma;

    const double f =
        std::exp(-0.5 * z * z) / (std::sqrt(2.0 * PI) * sigma);

    const double df =
        -(x - mu) / s2 * f;

    const double ddf =
        (z * z - 1.0) / s2 * f;

    return {f, df, ddf};
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