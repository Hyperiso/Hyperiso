#ifndef STANDARD_NORMAL_H
#define STANDARD_NORMAL_H

#include <random>
#include "IDistribution.h"

class StandardNormal final : public IDistribution {
public:
    explicit StandardNormal(unsigned int seed = std::random_device{}());

    Vector sample(std::size_t n) override;

private:
    std::mt19937 eng_;
    std::normal_distribution<double> dist_;
};

#endif