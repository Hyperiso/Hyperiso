#ifndef IBLOCKCOMPOSER_H
#define IBLOCKCOMPOSER_H

#include "Include.h"
#include "Parameter.h"
#include "DependentParameter.h"
#include "Block.h"

/**
 * @file IBlockComposer.h
 * @brief Interface for composing (creating/removing/updating) dependent blocks and parameters.
 *
 * This header defines @ref IBlockComposer, an abstraction used to build
 * dependent computation objects:
 *  - Dependent blocks (via @ref DependentBlock / @ref DependentBlockManager),
 *  - Dependent parameters (via @ref DependentParameter).
 *
 * The goal is to provide a *high-level orchestration API* for "composition":
 * - define what depends on what (sources),
 * - register the dependent entities into a destination ParameterType,
 * - allow later updates/removals in a uniform way.
 *
 * Implementations typically target a specific destination scope, e.g.:
 * - Wilson scope composer (@ref WilsonParamComposer) composing into ParameterType::WILSON.
 *
 * @see CompositeParamAdapter
 * @see DependentBlockManager
 * @see DependentParameter
 */

/**
 * @class IBlockComposer
 * @ingroup DependencyManagementModule
 * @brief Abstract interface for composing dependent blocks and parameters.
 *
 * A composer is responsible for:
 * - creating dependent blocks given a name and sources,
 * - creating dependent parameters given a ParamId and source ParamIds,
 * - triggering updates,
 * - removing composed entities,
 * - optionally removing all composed entities at once.
 *
 * The interface also provides a shared registry @ref composed_blocks used
 * by implementations to track which blocks were created during composition.
 */
class IBlockComposer {
public:
    virtual ~IBlockComposer() = default;

    /**
     * @brief Compose (create/register) a dependent block.
     *
     * @param block_name   Name of the dependent block to create.
     * @param source_names Map associating each source ParameterType to a list of source block names.
     * @param update_func  Function used to recalculate the dependent block.
     */
    virtual void compose_block(const std::string&, const std::unordered_map<ParameterType, std::vector<std::string>>&, const DepUpdateFunc&) = 0;

    /**
     * @brief Compose (create/register) a dependent parameter.
     *
     * @param pid         Target parameter identifier (block + code; type may be completed by implementation).
     * @param source_pids Set of parameter identifiers this parameter depends on.
     * @param update_func Function used to recalculate the dependent parameter value.
     */
    virtual void compose_parameter(const ParamId&, const std::unordered_set<ParamId>&, const DepParamUpdateFunc&) = 0;

    /**
     * @brief Remove a previously composed dependent block.
     *
     * @param block_name Name of the block to remove.
     */
    virtual void remove_block(const std::string&) = 0;

    /**
     * @brief Manually trigger an update for a previously composed block.
     *
     * @param block_name Name of the block to update.
     */
    virtual void update(const std::string&) = 0;

    /**
     * @brief Remove all blocks recorded as composed by this composer family.
     *
     * Implementations use @ref composed_blocks to track created blocks.
     */
    virtual void remove_all_composed_blocks() = 0;

    /**
     * @brief Returns whether a block has already been composed by any composer instance.
     *
     * This is intentionally lightweight and is used by setup hooks to stay idempotent
     * when a group is built more than once, e.g. active group + SM intermediate group.
     */
    bool has_composed_block(const std::string& block_name) const {
        return composed_blocks.find(block_name) != composed_blocks.end();
    }

protected:
    /**
     * @brief Registry of blocks composed through implementations of IBlockComposer.
     *
     * This shared set allows "bulk cleanup" operations (remove_all_composed_blocks).
     *
     * Note: being `static`, this is shared across *all instances* of composer classes.
     */
    inline static std::unordered_set<std::string> composed_blocks;
};

#endif // IBLOCKCOMPOSER_H
