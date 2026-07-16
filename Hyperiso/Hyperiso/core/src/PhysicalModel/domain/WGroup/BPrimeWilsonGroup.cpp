#include "BPrimeWilsonGroup.h"

#include <stdexcept>

BPrimeCoefficientGroup::BPrimeCoefficientGroup(WilsonGroupAdapterConfig adapters) : CoefficientGroup(adapters) {
    LOG_TRACE("In BPrimeCoefficientGroup constructor");
    this->id = GroupMapper::to_id(WGroup::BPrime);
}

std::shared_ptr<CoefficientGroup> BPrimeCoefficientGroup::clone() const {
    return std::make_shared<BPrimeCoefficientGroup>(*this);
}

std::unordered_map<WCoefId, scalar_t> BPrimeCoefficientGroup::base_1_LO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching,
    const BlockSrc& src
)
{
    auto ids = WCoefMapper::get_group(WGroup::BPrime);
    if (ids.size() != BRP::exp_prime_running.size()) {
        throw std::runtime_error("BPrime running exponent table size does not match the BPrime coefficient group size");
    }

    double eta = src.get_val("WPARAM_RUN_SM",2);
    
    std::unordered_map<WCoefId, scalar_t> Ci_run_map {};
    for (size_t k = 0; k < ids.size(); k++) {
        const auto id = WCoefMapper::to_id(ids[k]);
        Ci_run_map[id] = pow(eta, BRP::exp_prime_running[k]) * coef_matching.at(QCDOrder::LO).at(id);
    }

    return Ci_run_map;
}