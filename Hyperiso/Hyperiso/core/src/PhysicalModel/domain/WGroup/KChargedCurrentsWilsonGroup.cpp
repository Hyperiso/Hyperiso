#include "KChargedCurrentsWilsonGroup.h"

std::unordered_map<WCoef, scalar_t> KulnuCoefficientGroup::base_1_LO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching,
    const BlockSrc& src
) {
    auto ids = WCoefMapper::get_group(WGroup::CC_su);

    std::unordered_map<WCoef, scalar_t> Ci_run_map {};
    for (size_t k = 0; k < ids.size(); k++) {
        Ci_run_map[ids[k]] = coef_matching.at(QCDOrder::LO).at(ids[k]);
    }

    return Ci_run_map;
}

KulnuCoefficientGroup::KulnuCoefficientGroup(WilsonGroupAdapterConfig adapters, bool force_sm) : CoefficientGroup(adapters) {
    this->id = GroupMapper::to_id(WGroup::CC_su);
}


std::shared_ptr<CoefficientGroup> KulnuCoefficientGroup::clone() const {
    return std::make_shared<KulnuCoefficientGroup>(*this);
}