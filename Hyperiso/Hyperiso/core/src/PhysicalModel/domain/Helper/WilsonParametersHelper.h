#ifndef WILSON_PARAMETER_HELPER_H
#define WILSON_PARAMETER_HELPER_H

#include <array>
#include <memory>
#include <unordered_set>

#include "IBlockComposer.h"
#include "Include.h"
#include "IWilsonParameters.h"
#include "BWilsonRunningParameters.h"
#include "MesonMixingRunningParameters.h"
#include "QCDHelper.h"

using BRP = BWilsonRunningParameters;
using MMRP = MesonMixingRunningParameters;

/**
 * @file WilsonParameterHelper.h
 * @brief Default (SM) implementation of @ref IWilsonParameterHelper.
 *
 * This helper composes the standard Wilson auxiliary blocks used by built-in
 * computations:
 *
 * - "WPARAM_SI_SM"    : scale-independent constants (xh, generation, mf, sw2, beta0, ...)
 * - "WPARAM_MATCH_SM" : matching-scale dependent quantities at mu_W (alphas, xt, logs, masses, ...)
 * - "WPARAM_RUN_SM"   : running-scale dependent factors (alphas(mu), eta_5, eta_4, eta_3, ...)
 *
 * And depending on the group, it also composes running matrices blocks:
 * - For B/BPrime/BScalar groups:  "ETA_POWS", "U_MATRIX", "V_MATRIX"
 * - For MesonMixing group:        "ETA_POWS_MIXING", "UM_MATRIX_5", "UM_MATRIX_4"
 *
 * The blocks are registered via @ref IBlockComposer::compose_block.
 *
 * @see IWilsonParameterHelper
 * @see IBlockComposer
 * @see QCDHelper
 */
class WilsonParameterHelper : public IWilsonParameterHelper {
public:
    /**
     * @brief Constructs the helper.
     * @param ibc Block composer used to register dependent blocks.
     */
    WilsonParameterHelper(std::shared_ptr<IBlockComposer> ibc) : IWilsonParameterHelper(ibc) {}

    /**
     * @brief Initializes all required SM Wilson helper blocks.
     *
     * If called multiple times, the helper should not register duplicates.
     *
     * @param gen Generation index (used by some helper parameters).
     * @param grp Wilson group identifier (steers which running matrices are built).
     */
    void init(int gen, WGroupId grp) override;

    /**
     * @brief Resets the helper internal state.
     *
     * Resets the core/running initialization guards.
     */
    void cleanup() override;

    /**
     * @brief Composes the meson-mixing evolution matrices.
     *
     * Kept here so the physics formulae stay in one place, but called from the
     * MesonMixing group-definition setup hook rather than implicitly from every
     * WilsonParameterHelper::init() call.
     */
    static void compose_meson_mixing_running_blocks(const std::shared_ptr<IBlockComposer>& iblock_c);
    
private:
    /// Builds generation-dependent, scale-independent helper block ("WPARAM_SI_SM").
    void init_scale_independent_block(int gen) override;

    /// Builds matching-scale helper block ("WPARAM_MATCH_SM").
    void init_matching_block() override;

    /// Builds running-scale helper block ("WPARAM_RUN_SM") and relevant matrices blocks.
    void init_running_block(WGroupId grp) override;

    /// Builds B-sector running matrices ("ETA_POWS", "U_MATRIX", "V_MATRIX").
    void init_running_parameter_blocks_B();

    /// Builds meson-mixing running matrices ("ETA_POWS_MIXING", "UM_MATRIX_5", "UM_MATRIX_4").
    void init_running_parameter_blocks_MM();

    bool running_block_initialized{false};
    bool b_running_matrices_initialized{false};
    std::unordered_set<WGroupId> initialized_running_groups;

};

#endif // WILSON_PARAMETER_HELPER_H
