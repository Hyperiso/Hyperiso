#pragma once
#include <vector>
#include <random>
#include "Matrix.h"
#include "Include.h"


using Vec = std::vector<double>;

// Port: any sampler that can draw η ~ P(η; η̄, Σ)
class INuisanceSampler {
public:
    virtual ~INuisanceSampler() = default;
    virtual std::map<ParamId, double> sample() const = 0;
    virtual std::vector<std::map<ParamId, double>> sample(std::size_t n) const = 0;
};