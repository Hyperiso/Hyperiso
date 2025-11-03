#include "ChargedCurrentsWilsonGroup.h"

std::unordered_map<WCoef, scalar_t> BclnuCoefficientGroup::base_1_LO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching,
    const BlockSrc& src
) {
    auto ids = WCoefMapper::get_group(WGroup::CC_bc);

    std::unordered_map<WCoef, scalar_t> Ci_run_map {};
    for (size_t k = 0; k < ids.size(); k++) {
        Ci_run_map[ids[k]] = coef_matching.at(QCDOrder::LO).at(ids[k]);
    }
    return Ci_run_map;
}

BclnuCoefficientGroup::BclnuCoefficientGroup(WilsonGroupAdapterConfig adapters, bool force_sm) : CoefficientGroup(adapters) {
    this->id = GroupMapper::to_id(WGroup::CC_bc);
}


std::shared_ptr<CoefficientGroup> BclnuCoefficientGroup::clone() const {
    return std::make_shared<BclnuCoefficientGroup>(*this);
}


std::unordered_map<WCoef, scalar_t> BulnuCoefficientGroup::base_1_LO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching,
    const BlockSrc& src
) {
    auto ids = WCoefMapper::get_group(WGroup::CC_bu);

    std::unordered_map<WCoef, scalar_t> Ci_run_map {};
    for (size_t k = 0; k < ids.size(); k++) {
        Ci_run_map[ids[k]] = coef_matching.at(QCDOrder::LO).at(ids[k]);
    }

    return Ci_run_map;
}

BulnuCoefficientGroup::BulnuCoefficientGroup(WilsonGroupAdapterConfig adapters, bool force_sm) : CoefficientGroup(adapters) {
    this->id = GroupMapper::to_id(WGroup::CC_bu);
}


std::shared_ptr<CoefficientGroup> BulnuCoefficientGroup::clone() const {
    return std::make_shared<BulnuCoefficientGroup>(*this);
}