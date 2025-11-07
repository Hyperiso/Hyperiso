#ifndef DISTRIBUTION_FACTORY_H
#define DISTRIBUTION_FACTORY_H

#include <memory>
#include <random>
#include <algorithm>

#include "IDistribution.h"

class DistributionFactory {
public:
    static std::unique_ptr<IDistribution> create(const std::string& name,
                                                 unsigned int seed = std::random_device{}());
};

#endif