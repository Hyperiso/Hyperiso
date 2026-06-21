#ifndef THDM_PARAMETERS_HELPER_H
#define THDM_PARAMETERS_HELPER_H

#include <algorithm>
#include <array>
#include <functional>
#include <iostream>

#include "WilsonParamComposer.h"
#include "IWilsonParameters.h"
#include "Logger.h"
#include "SourcesView.h"

/**
 * @file THDMParametersHelper.h
 * @brief Wilson helper blocks specialized for the THDM model.
 *
 * This helper adds THDM-specific auxiliary parameter blocks used by
 * matching/running computations and some Wilson coefficient groups.
 *
 * Responsibilities:
 *  - build THDM scale-independent helper blocks (if any),
 *  - build THDM matching-scale helper blocks (if any),
 *  - build THDM running-scale helper blocks (if any),
 *  - guarantee idempotent initialization.
 *
 * It composes dependent blocks through the injected @ref IBlockComposer.
 *
 * @see IWilsonParameterHelper
 * @see IBlockComposer
 * @see WilsonParameterHelper
 */
class THDMParameterHelper : public IWilsonParameterHelper {
public:
    /**
     * @brief Constructs the helper.
     * @param ibc Block composer used to register dependent blocks/params.
     */
    THDMParameterHelper(std::shared_ptr<IBlockComposer> iblock_c) : IWilsonParameterHelper(iblock_c) {}

    /**
     * @brief Initializes THDM-specific helper blocks.
     *
     * Should be idempotent: if already initialized, does nothing.
     *
     * @param gen Generation index.
     * @param grp Wilson group identifier.
     */
    void init(int gen, WGroupId grp) override;

    /**
     * @brief Cleans up internal state and/or unregisters blocks if applicable.
     */
    void cleanup() override {}
protected:
    /**
     * @brief Builds THDM scale-independent helper blocks.
     */
    void init_scale_independent_block(int gen) override;

    /**
     * @brief Builds THDM matching-scale helper blocks.
     */
    void init_matching_block() override;

    /**
     * @brief Builds THDM running-scale helper blocks.
     */
    void init_running_block(WGroupId) override {}

};

#endif // THDM_PARAMETERS_HELPER_H