#ifndef LIKELIHOOD_DISCRETE_H
#define LIKELIHOOD_DISCRETE_H

#include <random>
#include <vector>
#include <numeric>
#include <stdexcept>
#include <cmath>
#include "IDistribution.h"

// Discrete likelihood sampler using Vose's alias method.
// values[i] drawn with probability proportional to weights[i].
class LikelihoodDiscrete final : public IDistribution {
public:
    LikelihoodDiscrete(std::vector<double> values,
                       std::vector<double> weights,
                       unsigned int seed = std::random_device{}(),
                       bool standardize = false);

    Vector sample(std::size_t n) override;

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
