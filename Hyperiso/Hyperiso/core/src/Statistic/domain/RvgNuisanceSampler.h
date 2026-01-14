#pragma once
#include <vector>
#include <random>
#include <stdexcept>
#include <cmath>

#include "INuisanceSampler.h"
#include "RNGHelper.h"
#include "JointDistribution.h"
#include "MarginalType.h"
#include "MarginalFactory.h"
#include "CopulaFactory.h"
#include "Indexing.h"

class RvgNuisanceSampler final : public INuisanceSampler {
public:
    explicit RvgNuisanceSampler(const std::vector<ParamId>& ids, std::unique_ptr<JointDistribution> rvg);
    std::map<ParamId, double> sample() const override;
    std::vector<std::map<ParamId, double>> sample(std::size_t n) const override;

private:
    std::unique_ptr<JointDistribution> rvg_;
    std::vector<ParamId> ids_;
};