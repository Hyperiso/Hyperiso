#include "BPrimeWilsonGroup.h"

BPrimeCoefficientGroup::BPrimeCoefficientGroup(WilsonGroupAdapterConfig adapters, bool force_sm) : CoefficientGroup(adapters) {
    LOG_TRACE("In BPrimeCoefficientGroup constructor");
    this->id = GroupMapper::to_id(WGroup::BPrime);
}

std::shared_ptr<CoefficientGroup> BPrimeCoefficientGroup::clone() const {
    return std::make_shared<BPrimeCoefficientGroup>(*this);
}

std::unordered_map<WCoef, scalar_t> BPrimeCoefficientGroup::base_1_LO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching,
    const BlockSrc& src
)
{
    std::array<complex_t, 12> CPi_match = {};
    auto ids = WCoefMapper::get_group(WGroup::BPrime);
    for (size_t k = 0; k < ids.size(); k++) {
        CPi_match[k] = coef_matching.at(QCDOrder::LO).at(ids[k]);
    }

    double eta = src.get_val("WPARAM_RUN_SM",2);
    
    std::unordered_map<WCoef, scalar_t> Ci_run_map {};
    for (size_t k = 0; k < ids.size(); k++) {
        Ci_run_map[ids[k]] = pow(eta, BRP::exp_prime_running[k]) * CPi_match[k];
    }

    return Ci_run_map;
}