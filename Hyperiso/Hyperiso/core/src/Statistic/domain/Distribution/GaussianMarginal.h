#ifndef STANDARD_NORMAL_H
#define STANDARD_NORMAL_H

#include <random>
#include "IMarginalDistribution.h"
#include "AbstractConfig.h"
#include "Math.h"

struct GaussianMarginalCfg : public AbstractConfig {
    double mu;
    double sigma;

    GaussianMarginalCfg(double mu, double sigma) : mu(mu), sigma(sigma) {}
};

class GaussianMarginal final : public IMarginalDistribution {
public:
    explicit GaussianMarginal(double mu, double sigma, unsigned int seed = std::random_device{}());

    Vector rvs(std::size_t n) override;
    double logpdf(double x) override;
    double cdf(double x) override;
    double ppf(double p) override;
    double mean() override;
    double std() override;

private:
    double mu, sigma;
    const gsl_rng_type* rng_tp {gsl_rng_mt19937};
    gsl_rng* eng_;
};

#endif