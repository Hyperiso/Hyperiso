#include "BNuNuWilson.h"
#include "wcoef_ids.hpp"

namespace {
constexpr double CL_SM_BNUNU = -6.32;
}

bool BNuNuWilson::is_left(WCoef coefficient) {
    switch (coefficient) {
        case WCoef::CNU_L_EE:
        case WCoef::CNU_L_EMU:
        case WCoef::CNU_L_ETAU:
        case WCoef::CNU_L_MUE:
        case WCoef::CNU_L_MUMU:
        case WCoef::CNU_L_MUTAU:
        case WCoef::CNU_L_TAUE:
        case WCoef::CNU_L_TAUMU:
        case WCoef::CNU_L_TAUTAU:
            return true;
        default:
            return false;
    }
}

bool BNuNuWilson::is_sm_diagonal_left(WCoef coefficient) {
    return coefficient == WCoef::CNU_L_EE
        || coefficient == WCoef::CNU_L_MUMU
        || coefficient == WCoef::CNU_L_TAUTAU;
}

BNuNuWilson::BNuNuWilson(WCoef coefficient, ContributionType contribution)
    : WilsonCoefficient(
          WCoefMapper::str(coefficient),
          GroupMapper::str(WGroup::BNuNu, ScaleType::MATCHING)),
      coefficient_(coefficient)
{
    set_contribution_type(contribution);

    const scalar_t base_value =
        contribution == ContributionType::SM && is_sm_diagonal_left(coefficient)
            ? CL_SM_BNUNU
            : 0.0;

    // HyperIso's getFull*() sums perturbative increments as
    // C_LO + alpha_s/(4pi) C_NLO + ... .  The quoted -6.32 is already the
    // phenomenological corrected SM coefficient, hence it belongs entirely in
    // the base/LO slot and higher increments are deliberately zero.
    matching_info[QCDOrder::LO] = {
        {},
        [base_value](const ParamSrc&) -> scalar_t { return base_value; },
        WCoefMapper::flha_full(coefficient, QCDOrder::LO, contribution)
    };
    matching_info[QCDOrder::NLO] = {
        {},
        [](const ParamSrc&) -> scalar_t { return 0.0; },
        WCoefMapper::flha_full(coefficient, QCDOrder::NLO, contribution)
    };
    matching_info[QCDOrder::NNLO] = {
        {},
        [](const ParamSrc&) -> scalar_t { return 0.0; },
        WCoefMapper::flha_full(coefficient, QCDOrder::NNLO, contribution)
    };
}
