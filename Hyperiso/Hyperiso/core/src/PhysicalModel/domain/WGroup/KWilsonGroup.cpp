#include "KWilsonGroup.h"

namespace {

std::unordered_map<WCoefId, scalar_t> k_identity_running_for_order(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching,
    QCDOrder order
) {
    // The K group currently has no non-trivial QCD evolution implemented in this
    // basis. The safest running map is therefore the identity at each
    // perturbative order:
    //
    //   C_i(mu_h)|_order = C_i(mu_W)|_order
    //
    // Important: do not use a fixed C1..C10 array here. WGroup::K contains its
    // own ids: CK9, CK10, CKQ1, CKQ2, CK_L, CPK9, CPK10, CPKQ1, CPKQ2.
    // It does not contain C1..C10.
    //
    // Also be tolerant with missing coefficients/orders: some K coefficients are
    // only defined at LO, while CK10 / CK_L can be requested at NLO depending on
    // the decay. Missing entries are interpreted as zero rather than throwing
    // unordered_map::at.
    const auto members = WCoefMapper::get_group(WGroup::K);

    std::unordered_map<WCoefId, scalar_t> out;
    out.reserve(members.size());

    const auto order_it = coef_matching.find(order);

    for (const auto coef : members) {
        const auto id = WCoefMapper::to_id(coef);
        scalar_t value = 0.0;

        if (order_it != coef_matching.end()) {
            const auto coef_it = order_it->second.find(id);
            if (coef_it != order_it->second.end()) {
                value = coef_it->second;
            }
        }

        out.emplace(id, value);
    }

    return out;
}

} // namespace

std::unordered_map<WCoefId, scalar_t> KCoefficientGroup::base_1_LO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching,
    const BlockSrc&
)
{
    return k_identity_running_for_order(coef_matching, QCDOrder::LO);
}

std::unordered_map<WCoefId, scalar_t> KCoefficientGroup::base_1_NLO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching,
    const BlockSrc&
)
{
    return k_identity_running_for_order(coef_matching, QCDOrder::NLO);
}

std::unordered_map<WCoefId, scalar_t> KCoefficientGroup::base_1_NNLO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching,
    const BlockSrc&
)
{
    return k_identity_running_for_order(coef_matching, QCDOrder::NNLO);
}