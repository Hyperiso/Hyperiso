#ifndef DISTRIBUTION_FACTORY_H
#define DISTRIBUTION_FACTORY_H

#include <memory>
#include <random>
#include <algorithm>
#include <variant>

#include "MarginalType.h"
#include "IMarginalDistribution.h"
#include "GaussianMarginal.h"
#include "FlatMarginal.h"
#include "LikelihoodMarginal.h"
#include "SplitGaussianMarginal.h"

using MarginalConfig = std::variant<FlatMarginalCfg, GaussianMarginalCfg, SplitGaussianMarginalCfg, LikelihoodMarginalCfg>;

class DistributionFactory {
public:
    static std::unique_ptr<IMarginalDistribution> create(MarginalType name, MarginalConfig cfg, unsigned int seed = std::random_device{}());
};

#endif