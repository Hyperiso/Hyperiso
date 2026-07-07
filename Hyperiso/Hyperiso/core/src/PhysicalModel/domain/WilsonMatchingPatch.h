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
 * @brief Hyperiso/SuperIso SM light-photon piece for C9 at LO.
 *
 * The MARTY C9 template removes photon-penguin diagrams with light up-type
 * quarks (u,c) in the loop.  This patch adds back the finite EFT matching term
 * that replaces those diagrams in the standard NDR convention:
 *
 *     + 38/27 - 4/9 L,    L = log(mu_W^2 / M_W^2)
 *
 * This is only the light photon-penguin finite piece.  It does not include the
 * SM charm/up box contribution (+1/(4 s_W^2)), which remains in the MARTY
 * calculation because only photon-mediated light diagrams are filtered.
 */
inline WilsonMatchingPatch make_hyperiso_c9_light_photon_sm_patch() {
    WilsonMatchingPatch patch;
    patch.group = GroupMapper::to_id(WGroup::B);
    patch.coefficient = WCoefMapper::to_id(WCoef::C9);
    patch.order = QCDOrder::LO;
    patch.contribution = ContributionType::SM;
    patch.sources = {
        ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(3)}
    };
    patch.compute = [](const ParamSrc& src) -> scalar_t {
        const scalar_t L = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 3);
        return scalar_t(38. / 27.) - scalar_t(4. / 9.) * L;
    };
    patch.label = "Hyperiso:C9:SM-light-photon";
    patch.marty_only = true;
    return patch;
}

#endif // HYPERISO_WILSON_MATCHING_PATCH_H
