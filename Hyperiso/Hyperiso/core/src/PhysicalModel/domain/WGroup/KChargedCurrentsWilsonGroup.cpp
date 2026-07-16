#include "KChargedCurrentsWilsonGroup.h"

std::unordered_map<WCoefId, scalar_t> KulnuCoefficientGroup::base_1_LO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching,
    const BlockSrc&
) {
    auto ids = WCoefMapper::get_group(WGroup::CC_su);

    std::unordered_map<WCoefId, scalar_t> Ci_run_map {};
    for (size_t k = 0; k < ids.size(); k++) {
        Ci_run_map[WCoefMapper::to_id(ids[k])] = coef_matching.at(QCDOrder::LO).at(WCoefMapper::to_id(ids[k]));
    }

    return Ci_run_map;
}

KulnuCoefficientGroup::KulnuCoefficientGroup(WilsonGroupAdapterConfig adapters) : CoefficientGroup(adapters) {
    this->id = GroupMapper::to_id(WGroup::CC_su);
}


std::shared_ptr<CoefficientGroup> KulnuCoefficientGroup::clone() const {
    return std::make_shared<KulnuCoefficientGroup>(*this);
}