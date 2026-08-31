#include "BNuNuWilsonGroup.h"

namespace {
std::unordered_map<WCoefId, scalar_t> identity_for_order(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& matching,
    QCDOrder order)
{
    const auto members = WCoefMapper::get_group(WGroup::BNuNu);
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
        out.emplace(id, value);
    }
    return out;
}
}

std::unordered_map<WCoefId, scalar_t> BNuNuCoefficientGroup::base_1_LO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& matching,
    const BlockSrc&)
{ return identity_for_order(matching, QCDOrder::LO); }

std::unordered_map<WCoefId, scalar_t> BNuNuCoefficientGroup::base_1_NLO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& matching,
    const BlockSrc&)
{ return identity_for_order(matching, QCDOrder::NLO); }

std::unordered_map<WCoefId, scalar_t> BNuNuCoefficientGroup::base_1_NNLO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& matching,
    const BlockSrc&)
{ return identity_for_order(matching, QCDOrder::NNLO); }

BNuNuCoefficientGroup::BNuNuCoefficientGroup(WilsonGroupAdapterConfig adapters)
    : CoefficientGroup(adapters)
{
    this->id = GroupMapper::to_id(WGroup::BNuNu);
}

std::shared_ptr<CoefficientGroup> BNuNuCoefficientGroup::clone() const {
    return std::make_shared<BNuNuCoefficientGroup>(*this);
}
