#include "DChargedCurrentsWilsonGroup.h"

std::unordered_map<WCoefId, scalar_t> DslnuCoefficientGroup::base_1_LO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching,
    const BlockSrc&
) {
    auto ids = WCoefMapper::get_group(WGroup::CC_cs);

    std::unordered_map<WCoefId, scalar_t> Ci_run_map {};
    for (size_t k = 0; k < ids.size(); k++) {
        Ci_run_map[WCoefMapper::to_id(ids[k])] = coef_matching.at(QCDOrder::LO).at(WCoefMapper::to_id(ids[k]));
    }

    return Ci_run_map;
}

DslnuCoefficientGroup::DslnuCoefficientGroup(WilsonGroupAdapterConfig adapters) : CoefficientGroup(adapters) {
    this->id = GroupMapper::to_id(WGroup::CC_cs);
}


std::shared_ptr<CoefficientGroup> DslnuCoefficientGroup::clone() const {
    return std::make_shared<DslnuCoefficientGroup>(*this);
}

std::unordered_map<WCoefId, scalar_t> DdlnuCoefficientGroup::base_1_LO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching,
    const BlockSrc&
) {
    auto ids = WCoefMapper::get_group(WGroup::CC_cd);

    std::unordered_map<WCoefId, scalar_t> Ci_run_map {};
    for (size_t k = 0; k < ids.size(); k++) {
        Ci_run_map[WCoefMapper::to_id(ids[k])] = coef_matching.at(QCDOrder::LO).at(WCoefMapper::to_id(ids[k]));
    }

    return Ci_run_map;
}

DdlnuCoefficientGroup::DdlnuCoefficientGroup(WilsonGroupAdapterConfig adapters) : CoefficientGroup(adapters) {
    this->id = GroupMapper::to_id(WGroup::CC_cd);
}


std::shared_ptr<CoefficientGroup> DdlnuCoefficientGroup::clone() const {
    return std::make_shared<DdlnuCoefficientGroup>(*this);
}