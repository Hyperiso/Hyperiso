#ifndef LIKELIHOOD_DISCRETE_H
#define LIKELIHOOD_DISCRETE_H

#include <random>
#include <vector>
#include <numeric>
#include <stdexcept>
#include <cmath>
#include "IMarginalDistribution.h"
#include "AbstractConfig.h"

struct LikelihoodMarginalCfg : public AbstractConfig {
    Vector values;
    Vector weights;
};

// Discrete likelihood sampler using Vose's alias method.
// values[i] drawn with probability proportional to weights[i].
class LikelihoodMarginal final : public IMarginalDistribution {
public:
    LikelihoodMarginal(Vector values,
                       Vector weights,
                       unsigned int seed = std::random_device{}(),
                       bool standardize = false);

    Vector rvs(std::size_t n) override;
    double logpdf(double x) override { return 0.0; }
    double cdf(double x) override { return 0.0; }
    double ppf(double p) override { return 0.0; }
    double mean() override { return 0.0; }
    double std() override { return 0.0; }

private:
    void build_alias_tables(std::vector<double> weights);

    std::mt19937 eng_;
    std::uniform_real_distribution<double> u01_;
    std::vector<double> values_;

    // Alias tables:
    std::vector<double> prob_;   // in [0,1]
    std::vector<std::size_t> alias_;

    bool standardize_;
    double mean_ = 0.0;
    double std_  = 1.0;
};

#endif
