#include "RvgNuisanceSampler.h"

RvgNuisanceSampler::RvgNuisanceSampler(const std::vector<ParamId> &ids, std::unique_ptr<JointDistribution> rvg)
   : ids_(ids), rvg_(std::move(rvg))    
{}

std::map<ParamId, double> RvgNuisanceSampler::sample() const {
    std::vector<double> s = rvg_->sample();
    return zip(ids_, s);
}

std::vector<std::map<ParamId, double>> RvgNuisanceSampler::sample(std::size_t n) const {
    std::vector<std::vector<double>> samples = rvg_->sample(n);
    std::vector<std::map<ParamId, double>> out;
    for (const auto& s : samples) {
        out.emplace_back(zip(ids_, s));
    }
    return out;
}
