#ifndef STANDARD_FLAT_H
#define STANDARD_FLAT_H

#include <random>
#include <cmath>
#include "IMarginalDistribution.h"
#include "Include.h"
#include "AbstractConfig.h"

struct FlatMarginalCfg : public AbstractConfig {
    double a;
    double b;

    FlatMarginalCfg(double a, double b) : a(a), b(b) {}
};

class FlatMarginal final : public IMarginalDistribution {
public:
    explicit FlatMarginal(double a, double b, unsigned int seed = std::random_device{}());

    Vector rvs(std::size_t n) override;
    double logpdf(double x) override;
    double cdf(double x) override;
    double ppf(double p) override;
    double mean() override;
    double std() override;

private:
    double a, b;
    const gsl_rng_type* rng_tp {gsl_rng_mt19937};
    gsl_rng* eng_;
};

#endif