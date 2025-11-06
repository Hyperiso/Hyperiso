#pragma once
#include <vector>
#include <random>


using Vec = std::vector<double>;
using Matrix = std::vector<std::vector<double>>;


// Port: any sampler that can draw η ~ P(η; η̄, Σ)
class INuisanceSampler {
public:
virtual ~INuisanceSampler() = default;
virtual Vec sample(const Vec& mean, const Matrix& cov, std::mt19937& rng) const = 0;
};