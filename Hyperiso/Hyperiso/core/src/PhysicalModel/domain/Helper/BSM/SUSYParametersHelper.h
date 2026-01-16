#ifndef SUSY_PARAMETERS_HELPER_H
#define SUSY_PARAMETERS_HELPER_H

#include <algorithm>
#include <array>
#include <functional>

#include "Math_SUSY.h"
#include "Logger.h"
#include "ParameterProxy.h"
#include "WilsonParamComposer.h"
#include "IWilsonParameters.h"

/**
 * @file SUSYParametersHelper.h
 * @brief SUSY-specific implementation of @ref IWilsonParameterHelper.
 *
 * This helper is responsible for building auxiliary Wilson parameter blocks
 * required when the active model is SUSY.
 *
 * It composes:
 *  - a SUSY scale-independent block (implementation-specific),
 *  - a SUSY matching block (implementation-specific),
 *  - an "EPSILON_SUSY" block used in SUSY loop corrections / effective couplings.
 *
 * The "EPSILON_SUSY" block is built from a mix of:
 *  - SM blocks (MASS, SMINPUTS),
 *  - SUSY/BSM blocks (MASS, GAUGE, HMIX, MSOFT, AD, AU, YD, YU, mixings, ...),
 *  - and Wilson helper blocks (e.g. WPARAM_SI_SM).
 *
 * Composition is performed via the injected @ref IBlockComposer.
 *
 * @see IWilsonParameterHelper
 * @see IBlockComposer
 * @see DependentBlock
 */

constexpr double Pi = 3.14159265358979323846;

constexpr int N_UL_UR = 7, M_UL_UR = 4, N_NL_NR = 4, M_NL_NR = 4, N_Gamma_U = 7, M_Gamma_U = 7, N_X = 3, N_Mch = 3, N_MsqU = 6, N_MsqD = 6, N_Msn = 4;

using Array2D_7x4 = std::array<std::array<double, M_UL_UR>, N_UL_UR>;
using Array2D_4x4 = std::array<std::array<double, M_NL_NR>, N_NL_NR>;
using Array2D_7x7 = std::array<std::array<double, M_Gamma_U>, N_Gamma_U>;
using Array3D_3x7x4 = std::array<Array2D_7x4, N_X>;
using Array1D_4 = std::array<double, 4>;
using Array1D_3 = std::array<double, N_Mch>;
using Array1D_7 = std::array<double, N_MsqU>;
using Array2D_4x4_I = std::array<std::array<complex_t, M_NL_NR>, N_NL_NR>;

/**
 * @class susy_parameters
 * @brief Wilson helper blocks builder for the SUSY model.
 *
 * The helper is idempotent: init() returns immediately if already initialized.
 *
 * Note:
 * - init_running_block() is currently empty for this helper (no running blocks),
 *   but the hook is kept to respect the interface.
 */
class SUSYParameterHelper : public IWilsonParameterHelper {
public:
	/**
     * @brief Constructs the helper.
     * @param ibc Block composer used to register dependent blocks/params.
     */
	SUSYParameterHelper(std::shared_ptr<IBlockComposer> iblock_c) : IWilsonParameterHelper(iblock_c) {}

	/**
     * @brief Initializes SUSY-specific helper blocks.
     *
     * Should be idempotent: if already initialized, does nothing.
     *
     * @param gen Generation index (may be forwarded to some helper blocks).
     * @param grp Wilson group identifier (may control which running blocks are needed).
     */
	void init(int gen, WGroupId grp) override;

	/**
     * @brief Cleans up internal state and/or unregisters blocks if applicable.
     */
	void cleanup() override {}

protected:
	void init_epsilon_block();

	/**
     * @brief Builds SUSY scale-independent helper blocks.
     */
    void init_scale_independent_block(int gen) override;

	/**
     * @brief Builds SUSY matching-scale helper blocks.
     */
    void init_matching_block() override;

	/**
     * @brief Builds SUSY running-scale helper blocks.
     */
	void init_running_block(WGroupId) override {};
};

#endif