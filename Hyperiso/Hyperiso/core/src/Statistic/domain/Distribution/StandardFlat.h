#ifndef STANDARD_FLAT_H
#define STANDARD_FLAT_H

#include <random>
#include <cmath>
#include "IDistribution.h"

class StandardFlat final : public IDistribution {
public:
    explicit StandardFlat(unsigned int seed = std::random_device{}());

    Vector sample(std::size_t n) override;

private:
    std::mt19937 eng_;
    std::uniform_real_distribution<double> dist_;
};

#endif