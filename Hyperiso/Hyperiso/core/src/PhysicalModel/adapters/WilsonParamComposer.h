#ifndef WILSONPARAMCOMPOSER_H
#define WILSONPARAMCOMPOSER_H

#include "IBlockComposer.h"
#include "CompositeParamCreator.h"

/**
 * @file WilsonParamComposer.h
 * @brief Block composer specialized for the Wilson parameter scope.
 *
 * This header defines @ref WilsonParamComposer, a concrete implementation
 * of @ref IBlockComposer that composes dependent blocks and parameters into
 * the WILSON Parameters instance (ParameterType::WILSON).
 *
 * Internally, this class delegates to @ref CompositeParamAdapter, which itself
 * adapts the high-level dependency operations to @ref DependentBlockManager.
 *
 * @see IBlockComposer
 * @see CompositeParamAdapter
 * @see DependentBlockManager
 */

/**
 * @class WilsonParamComposer
 * @ingroup DependencyManagementModule
 * @brief Composer that registers dependent entities into ParameterType::WILSON.
 *
 * Responsibilities:
 * - Compose dependent blocks into the Wilson scope (dest = ParameterType::WILSON).
 * - Compose dependent parameters, ensuring the target ParamId is typed as WILSON.
 * - Track created block names through IBlockComposer::composed_blocks so they can be removed later.
 *
 * Notes on typing behavior:
 * - compose_parameter() ensures the destination pid.type is ParameterType::WILSON.
 * - for source ParamIds:
 *   - if the source type is missing, it is assumed to be ParameterType::WILSON,
 *   - otherwise, the provided type is preserved.
 */
class WilsonParamComposer : public IBlockComposer {
public:
    /**
     * @brief Creates a dependent block in the Wilson scope.
     *
     * Delegates to:
     * @code
     * CompositeParamAdapter().add_block_dependency(..., ParameterType::WILSON, ...);
     * @endcode
     *
     * Also records @p block_name into IBlockComposer::composed_blocks.
     */
    void compose_block(const std::string& block_name, const std::unordered_map<ParameterType, std::vector<std::string>>& source_names, const DepUpdateFunc& update_func) override;

    /**
     * @brief Creates a dependent parameter in the Wilson scope.
     *
     * - Forces destination pid.type to ParameterType::WILSON.
     * - Normalizes sources:
     *   - missing type => assumed WILSON
     *   - specified type => preserved
     *
     * Also records the destination block name into IBlockComposer::composed_blocks.
     */
    void compose_parameter(const ParamId&, const std::unordered_set<ParamId>&, const DepParamUpdateFunc&) override;

    /**
     * @brief Removes a dependent block from the Wilson scope and from the registry.
     */
    void remove_block(const std::string&) override;

    /**
     * @brief Triggers an update on a dependent block in the Wilson scope.
     */
    void update(const std::string& block_name) override;

    /**
     * @brief Removes all dependent blocks composed via this composer.
     *
     * Iterates over IBlockComposer::composed_blocks and removes each block from
     * the Wilson scope, then clears the registry.
     */
    void remove_all_composed_blocks() override;

};

#endif // WILSONPARAMCOMPOSER_H
