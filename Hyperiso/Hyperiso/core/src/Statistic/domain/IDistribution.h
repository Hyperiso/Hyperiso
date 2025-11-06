#ifndef IDISTRIBUTION_H
#define IDISTRIBUTION_H

#include "RNGHelper.h"

struct IDistribution {
    virtual ~IDistribution() = default;

    virtual Vector sample(std::size_t n) = 0;
};

#endif