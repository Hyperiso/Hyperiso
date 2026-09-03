#include "KNuNuWilsonGroup.h"

namespace {

scalar_t perturbative_identity_rescaling(QCDOrder order, const BlockSrc& src) {
    if (order == QCDOrder::LO) return 1.0;
    const double alpha_match = src.get_val("WPARAM_MATCH_SM", 1);
    const double alpha_run = src.get_val("WPARAM_RUN_SM", 1);
    if (alpha_run == 0.0) {
        LOG_ERROR("DivisionByZero", "KNuNu identity running received alpha_s(mu_h)=0");
    }
    const double ratio = alpha_match / alpha_run;
    if (order == QCDOrder::NLO) return ratio;
    if (order == QCDOrder::NNLO) return ratio * ratio;
    return 1.0;
}

std::unordered_map<WCoefId, scalar_t> identity_for_order(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& matching,
    QCDOrder order,
    const BlockSrc& src)
{
    const auto members = WCoefMapper::get_group(WGroup::KNuNu);
    const scalar_t rescaling = perturbative_identity_rescaling(order, src);
    std::unordered_map<WCoefId, scalar_t> out;
    out.reserve(members.size());
    const auto order_it = matching.find(order);

    for (WCoef coefficient : members) {
        const auto id = WCoefMapper::to_id(coefficient);
        scalar_t value = 0.0;
        if (order_it != matching.end()) {
            const auto value_it = order_it->second.find(id);
            if (value_it != order_it->second.end()) value = value_it->second;
        }
        out.emplace(id, value * rescaling);
    }
    return out;
}

} // namespace

std::unordered_map<WCoefId, scalar_t> KNuNuCoefficientGroup::base_1_LO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& matching,
    const BlockSrc& src)
{ return identity_for_order(matching, QCDOrder::LO, src); }

std::unordered_map<WCoefId, scalar_t> KNuNuCoefficientGroup::base_1_NLO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& matching,
    const BlockSrc& src)
{ return identity_for_order(matching, QCDOrder::NLO, src); }

std::unordered_map<WCoefId, scalar_t> KNuNuCoefficientGroup::base_1_NNLO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& matching,
    const BlockSrc& src)
{ return identity_for_order(matching, QCDOrder::NNLO, src); }

KNuNuCoefficientGroup::KNuNuCoefficientGroup(WilsonGroupAdapterConfig adapters)
    : CoefficientGroup(adapters)
{
    this->id = GroupMapper::to_id(WGroup::KNuNu);
}

std::shared_ptr<CoefficientGroup> KNuNuCoefficientGroup::clone() const {
    return std::make_shared<KNuNuCoefficientGroup>(*this);
}
