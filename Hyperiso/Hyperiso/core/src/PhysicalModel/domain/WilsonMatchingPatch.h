#ifndef HYPERISO_WILSON_MATCHING_PATCH_H
#define HYPERISO_WILSON_MATCHING_PATCH_H

#include <functional>
#include <string>
#include <unordered_set>

#include "Include.h"
#include "Math.h"
#include "SourcesView.h"

/**
 * @file WilsonMatchingPatch.h
 * @brief Small extension point for additive Wilson matching terms.
 *
 * A WilsonMatchingPatch decorates one coefficient/order/contribution with an
 * extra matching-scale term.  It is intentionally local to the matching layer:
 * the usual SM/BSM/TOTAL composition and hadronic running then see the patched
 * matching coefficient exactly like any other Wilson coefficient.
 */
using WilsonMatchingPatchFunction = std::function<scalar_t(const ParamSrc&)>;

struct WilsonMatchingPatch {
    WGroupId group{};
    WCoefId coefficient{};
    QCDOrder order{QCDOrder::LO};
    ContributionType contribution{ContributionType::SM};
    std::unordered_set<ParamId> sources{};
    WilsonMatchingPatchFunction compute{};
    std::string label{};
    bool enabled{true};

    /**
     * If true, apply the patch only to MARTY-backed coefficients.  This avoids
     * double counting when the same SM contribution is supplied by Hyperiso's
     * built-in analytic backend.
     */
    bool marty_only{false};
};

/**
 * @brief Optional HyperIso/SuperIso SM photon piece for C9 at LO.
 *
 * The production MARTY path now uses the builtin HyperIso/SuperIso SM C9,
 * so this helper is not registered automatically.  It is kept as an explicit
 * building block for workflows that deliberately compute the non-photon SM
 * C9 pieces with MARTY and want to add the finite photon matching term by hand:
 *
 *     -D0t(x_t) + 38/27 - 4/9 L,    L = log(mu_W^2 / M_W^2).
 *
 * For BSM photon penguins there is no universal replacement; users should add
 * a model-specific WilsonMatchingPatch.
 */
inline WilsonMatchingPatch make_hyperiso_c9_sm_photon_patch() {
    WilsonMatchingPatch patch;
    patch.group = GroupMapper::to_id(WGroup::B);
    patch.coefficient = WCoefMapper::to_id(WCoef::C9);
    patch.order = QCDOrder::LO;
    patch.contribution = ContributionType::SM;
    patch.sources = {
        ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2, 1)},
        ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(3)}
    };
    patch.compute = [](const ParamSrc& src) -> scalar_t {
        const double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2, 1));
        const scalar_t L = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(3));
        return scalar_t(-D0t(xt) + 38. / 27.) - scalar_t(4. / 9.) * L;
    };
    patch.label = "Hyperiso:C9:SM-photon";
    patch.marty_only = true;
    return patch;
}

// Backward-compatible alias for older user code.  The old name only described
// the light u/c part; the helper now returns the full SM photon matching piece.
inline WilsonMatchingPatch make_hyperiso_c9_light_photon_sm_patch() {
    return make_hyperiso_c9_sm_photon_patch();
}

/**
 * @brief THDM charged-Higgs photon-penguin piece for C9 at LO.
 *
 * The MARTY C9 BSM template vetoes the massless photon linker to avoid the
 * reg_prop-dependent four-fermion projection.  In the THDM this finite photon
 * matching contribution is known analytically and is the -D9H0(y_t, lambda_u)
 * term in the HyperIso/SuperIso C9_THDM coefficient.
 */
inline WilsonMatchingPatch make_hyperiso_c9_thdm_photon_patch() {
    WilsonMatchingPatch patch;
    patch.group = GroupMapper::to_id(WGroup::B);
    patch.coefficient = WCoefMapper::to_id(WCoef::C9);
    patch.order = QCDOrder::LO;
    patch.contribution = ContributionType::BSM;
    patch.sources = {
        ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", LhaID(7)},
        ParamId{ParameterType::WILSON, "WPARAM_MATCH_BSM", LhaID(1)}
    };
    patch.compute = [](const ParamSrc& src) -> scalar_t {
        const double lu = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", LhaID(7));
        const double yt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", LhaID(1));
        return scalar_t(-D9H0(yt, lu));
    };
    patch.label = "Hyperiso:C9:THDM-photon";
    patch.marty_only = true;
    return patch;
}

/**
 * @brief THDM charged-Higgs photon-penguin piece for CP9 at LO.
 *
 * This is only the finite photon term removed from the MARTY CP9 template,
 * not the full CP9_THDM coefficient.
 */
inline WilsonMatchingPatch make_hyperiso_cp9_thdm_photon_patch() {
    WilsonMatchingPatch patch;
    patch.group = GroupMapper::to_id(WGroup::BPrime);
    patch.coefficient = WCoefMapper::to_id(WCoef::CP9);
    patch.order = QCDOrder::LO;
    patch.contribution = ContributionType::BSM;
    patch.sources = {
        ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", LhaID(8)},
        ParamId{ParameterType::WILSON, "WPARAM_MATCH_BSM", LhaID(1)},
        ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(6)},
        ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)},
        ParamId{ParameterType::SM, "MASS", LhaID(3)}
    };
    patch.compute = [](const ParamSrc& src) -> scalar_t {
        const double ld     = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", LhaID(8));
        const double yt     = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_BSM", LhaID(1));
        const double mt_muW = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(6));
        const double mb_muW = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1));
        const double ms     = src.get_val(ParameterType::SM, "MASS", LhaID(3));
        return scalar_t(-ms * mb_muW / (mt_muW * mt_muW) * D9H0(yt, ld));
    };
    patch.label = "Hyperiso:CP9:THDM-photon";
    patch.marty_only = true;
    return patch;
}

#endif // HYPERISO_WILSON_MATCHING_PATCH_H
