#ifndef BLOCK_PROXY_H
#define BLOCK_PROXY_H

#include <unordered_set>
#include <string>

#include "IBlockProxy.h"
#include "BlockProvider.h"
#include "APIAdapter.h"

/**
 * @file BlockProxy.h
 * @brief Facade combining block existence/logging (BlockProvider) and block listing (APIAdapter).
 *
 * This header defines ::BlockProxy, a lightweight proxy that unifies:
 *  - block existence checks and logging via ::BlockProvider,
 *  - block list retrieval via ::APIAdapter.
 *
 * Rationale:
 * - ::BlockProvider is focused on *logging* and existence checks bound to ::Parameters.
 * - ::APIAdapter already implements robust block listing logic (DataMonitoringModule),
 *   so ::BlockProxy delegates "get block list" to it.
 *
 * Note about logging:
 * - Block logging ultimately relies on ::Parameters::print_block(blockname), which streams
 *   the underlying ::BlockAccessor content (via `operator<<` on `std::shared_ptr<BlockAccessor>`).
 *   This prints:
 *   - each block name,
 *   - each parameter entry (id/value) inside the block.
 *
 * @see BlockProvider
 * @see APIAdapter
 * @see Parameters
 * @see BlockAccessor
 */

/**
 * @class BlockProxy
 * @ingroup DataMonitoringModule
 * @brief Unified proxy for block queries and diagnostics.
 *
 * ::BlockProxy exposes a single entry point to:
 *  - test whether a block exists for a given ::ParameterType,
 *  - print either all blocks of a model or a specific block,
 *  - retrieve the list of block names for a model.
 *
 * Internally:
 *  - existence + logs are delegated to ::BlockProvider,
 *  - listing is delegated to ::APIAdapter (which uses the monitoring layer).
 */
class BlockProxy : public IBlockProxy {
public:
    /**
     * @brief Default constructor.
     *
     * Initializes the internal ::BlockProvider and ::APIAdapter members.
     */
    BlockProxy() : bp(), api() {}

    /**
     * @brief Checks whether a block exists for the given ParameterType.
     *
     * Delegates to ::BlockProvider::exists().
     *
     * @param blockname Name of the block to check.
     * @param pt ParameterType for which to check the block.
     * @return true if the block exists, false otherwise.
     */
    bool exists(const std::string& blockname, ParameterType pt) override;

    /**
     * @brief Logs all blocks for a given ParameterType.
     *
     * Delegates to ::BlockProvider::log_all_blocks().
     *
     * Depending on the logging backend, this typically prints the content of the
     * associated ::Parameters instance (which includes its ::BlockAccessor).
     *
     * @param pt ParameterType whose blocks should be logged.
     */
    void log_all_blocks(ParameterType pt) override;

    /**
     * @brief Logs a single block for a given ParameterType.
     *
     * Delegates to ::BlockProvider::log_block(), which itself calls
     * ::Parameters::print_block(blockname).
     *
     * Since `print_block()` streams `this->blockAccessor->at(blockname)`,
     * the output is produced using the `operator<<` overload for
     * `std::shared_ptr<BlockAccessor>`.
     *
     * @param pt ParameterType owning the block.
     * @param blockname Name of the block to log.
     */
    void log_block(ParameterType pt, const std::string& blockname) override;

    /**
     * @brief Retrieves the set of block names for a given ParameterType.
     *
     * Delegates to ::APIAdapter::get_blocks_list(pt), which provides the
     * monitoring-layer implementation for retrieving block names.
     *
     * @param pt ParameterType to query.
     * @return Unordered set of block names owned by @p pt.
     */
    std::unordered_set<BlockName> get_block_list(ParameterType pt) override;

    /**
     * @brief Retrieves the content of a block for a given ParameterType.
     *
     * Delegates to ::BlockProvider::get_block(pt), which provides the
     * monitoring-layer implementation for retrieving block contents.
     *
     * @param pt ParameterType to query.
     * @param blockname Name of the block.
     * @return map of theblock content owned by @p pt.
     */
    std::map<LhaID, scalar_t> get_block(ParameterType pt, const std::string& blockname) override;

private:
    /// Provider used for existence checks and logging.
    BlockProvider bp;

    /// Adapter used for retrieving block lists from the monitoring layer.
    APIAdapter api;
};

#endif