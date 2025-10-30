#include "DChargedCurrentsWilsonGroup.h"

std::unordered_map<WCoef, scalar_t> DslnuCoefficientGroup::base_1_LO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching,
    const std::unordered_map<std::string, std::shared_ptr<Block>>& src
) {
    auto ids = WCoefMapper::get_group(WGroup::BCC_cs);

    std::unordered_map<WCoef, scalar_t> Ci_run_map {};
    for (size_t k = 0; k < ids.size(); k++) {
        Ci_run_map[ids[k]] = coef_matching.at(QCDOrder::LO).at(ids[k]);
    }

    return Ci_run_map;
}

DslnuCoefficientGroup::DslnuCoefficientGroup(WilsonGroupAdapterConfig adapters, bool force_sm) : CoefficientGroup(adapters) {
    this->id = GroupMapper::to_id(WGroup::BCC_cs);
}


std::shared_ptr<CoefficientGroup> DslnuCoefficientGroup::clone() const {
    return std::make_shared<DslnuCoefficientGroup>(*this);
}

std::unordered_map<WCoef, scalar_t> DdlnuCoefficientGroup::base_1_LO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching,
    const std::unordered_map<std::string, std::shared_ptr<Block>>& src
) {
    auto ids = WCoefMapper::get_group(WGroup::BCC_cd);

    std::unordered_map<WCoef, scalar_t> Ci_run_map {};
    for (size_t k = 0; k < ids.size(); k++) {
        Ci_run_map[ids[k]] = coef_matching.at(QCDOrder::LO).at(ids[k]);
    }

    return Ci_run_map;
}

DdlnuCoefficientGroup::DdlnuCoefficientGroup(WilsonGroupAdapterConfig adapters, bool force_sm) : CoefficientGroup(adapters) {
    this->id = GroupMapper::to_id(WGroup::BCC_cd);
}


std::shared_ptr<CoefficientGroup> DdlnuCoefficientGroup::clone() const {
    return std::make_shared<DdlnuCoefficientGroup>(*this);
}