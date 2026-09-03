#include "KNuNuWilson.h"
#include "KWilson.h"
#include "wcoef_ids.hpp"

bool KNuNuWilson::is_left(WCoef coefficient) {
    return coefficient == WCoef::CKNU_L_EE
        || coefficient == WCoef::CKNU_L_MUMU
        || coefficient == WCoef::CKNU_L_TAUTAU;
}

bool KNuNuWilson::is_sm_diagonal_left(WCoef coefficient) {
    return is_left(coefficient);
}

KNuNuWilson::KNuNuWilson(WCoef coefficient, ContributionType contribution)
    : WilsonCoefficient(
          WCoefMapper::str(coefficient),
          GroupMapper::str(WGroup::KNuNu, ScaleType::MATCHING)),
      coefficient_(coefficient)
{
    set_contribution_type(contribution);

    const bool native_sm_left =
        contribution == ContributionType::SM && is_sm_diagonal_left(coefficient);

    if (native_sm_left) {
        // Use the already validated native SM CK_L matching formula as the
        // implementation of the SM short-distance function, but expose it in
        // the flavour-resolved KNuNu basis.  KPinunuDecay itself no longer
        // reads the legacy universal CK_L coefficient.
        matching_info[QCDOrder::LO] = {
            {
                {"WPARAM_MATCH_SM", LhaID(2, 1)},
                {"WPARAM_SI_SM", 4}
            },
            CK_L::compute_LO,
            WCoefMapper::flha_full(coefficient, QCDOrder::LO, contribution)
        };
        matching_info[QCDOrder::NLO] = {
            {
                {"WPARAM_MATCH_SM", LhaID(2, 1)},
                {"WPARAM_SI_SM", 4},
                {ParameterType::SM, "MASS", 24},
                {ParameterType::SM, "QCD", 6}
            },
            CK_L::compute_NLO,
            WCoefMapper::flha_full(coefficient, QCDOrder::NLO, contribution)
        };
    } else {
        matching_info[QCDOrder::LO] = {
            {},
            [](const ParamSrc&) -> scalar_t { return 0.0; },
            WCoefMapper::flha_full(coefficient, QCDOrder::LO, contribution)
        };
        matching_info[QCDOrder::NLO] = {
            {},
            [](const ParamSrc&) -> scalar_t { return 0.0; },
            WCoefMapper::flha_full(coefficient, QCDOrder::NLO, contribution)
        };
    }

    matching_info[QCDOrder::NNLO] = {
        {},
        [](const ParamSrc&) -> scalar_t { return 0.0; },
        WCoefMapper::flha_full(coefficient, QCDOrder::NNLO, contribution)
    };
}
