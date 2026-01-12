#include "RvgNuisanceSampler.h"

RvgNuisanceSampler::RvgNuisanceSampler(const std::vector<ParamId> &ids, std::unique_ptr<JointDistribution> rvg)
   : ids_(ids), rvg_(std::move(rvg))    
{}

std::map<ParamId, double> RvgNuisanceSampler::sample() const {
    Vector s = rvg_->sample();
    return zip(ids_, s);
}
