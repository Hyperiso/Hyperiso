#ifndef DEPENDANT_BLOCK_INFO_PROVIDER_H
#define DEPENDANT_BLOCK_INFO_PROVIDER_H

#include <string>
#include <vector>

#include "Include.h"
#include "Parameters.h"

/**
 * @file DependantBlockInfoProvider.h
 * @brief High-level adapter for dependent-block graph introspection.
 *
 * This provider exposes a compact API to inspect block-level dependencies from
 * Python and C++ without requiring callers to manipulate @ref Parameters,
 * @ref BlockAccessor, @ref Block, or @ref DependentBlock objects directly.
 *
 * The inspected graph is the runtime graph currently held by the selected
 * @ref Parameters repository:
 * - upstream/source links are read from @ref DependentBlock::get_source_blocks,
 * - downstream/dependent links are read from @ref Block::getObservers.
 *
 * The API is intentionally name-based: callers provide a @ref ParameterType and
 * a block name, and receive either a boolean or lists of block names.
 */

/**
 * @class DependantBlockInfoProvider
 * @ingroup DataProvidersModule
 * @brief Adapter exposing dependent-block metadata for one parameter namespace.
 *
 * Typical usage:
 * @code
 * DependantBlockInfoProvider provider;
 * bool is_dep = provider.is_dependent_block(ParameterType::SM, "EW");
 * auto sources = provider.get_source_blocks(ParameterType::SM, "EW");
 * auto dependents = provider.get_dependent_blocks(ParameterType::SM, "SMINPUTS");
 * @endcode
 *
 * The direct methods return only immediate graph neighbours. The `get_all_*`
 * methods recursively walk the graph and deduplicate block names.
 */
class DependantBlockInfoProvider {
public:
    /**
     * @brief Checks whether a block is a @ref DependentBlock.
     *
     * @param type Parameter namespace containing the block to inspect.
     * @param block_name Block name or alias to inspect.
     * @return True if the resolved block is a @ref DependentBlock.
     *
     * @throws std::invalid_argument if the block cannot be resolved.
     */
    bool is_dependent_block(ParameterType type, const std::string& block_name) const;

    /**
     * @brief Returns the direct source blocks of a block.
     *
     * This answers: "which blocks does this block depend on directly?"
     * Plain blocks return an empty list.
     *
     * @param type Parameter namespace containing the block to inspect.
     * @param block_name Block name or alias to inspect.
     * @return Sorted list of direct upstream/source block names.
     *
     * @throws std::invalid_argument if the block cannot be resolved.
     */
    std::vector<std::string> get_source_blocks(ParameterType type, const std::string& block_name) const;

    /**
     * @brief Returns the direct blocks depending on a block.
     *
     * This answers: "which blocks directly depend on this block?"
     * It follows direct block observers registered on the inspected block.
     *
     * @param type Parameter namespace containing the block to inspect.
     * @param block_name Block name or alias to inspect.
     * @return Sorted list of direct downstream/dependent block names.
     *
     * @throws std::invalid_argument if the block cannot be resolved.
     */
    std::vector<std::string> get_dependent_blocks(ParameterType type, const std::string& block_name) const;

    /**
     * @brief Returns all transitive source blocks of a block.
     *
     * This recursively follows upstream dependency links and returns each block
     * name at most once.
     *
     * @param type Parameter namespace containing the block to inspect.
     * @param block_name Block name or alias to inspect.
     * @return Sorted list of all upstream/source block names.
     *
     * @throws std::invalid_argument if the block cannot be resolved.
     */
    std::vector<std::string> get_all_source_blocks(ParameterType type, const std::string& block_name) const;

    /**
     * @brief Returns all transitive blocks depending on a block.
     *
     * This recursively follows downstream observer links and returns each block
     * name at most once.
     *
     * @param type Parameter namespace containing the block to inspect.
     * @param block_name Block name or alias to inspect.
     * @return Sorted list of all downstream/dependent block names.
     *
     * @throws std::invalid_argument if the block cannot be resolved.
     */
    std::vector<std::string> get_all_dependent_blocks(ParameterType type, const std::string& block_name) const;
};

#endif // DEPENDANT_BLOCK_INFO_PROVIDER_H
