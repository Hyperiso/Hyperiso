#ifndef IBLOCK_PROXY_H
#define IBLOCK_PROXY_H

#include <unordered_set>
#include <string>

#include "Include.h"
#include "scalar.h"

/**
 * @file IBlockProxy.h
 * @brief Interface (port) for unified access to parameter block diagnostics.
 *
 * This header defines @ref IBlockProxy, an abstract interface exposing a
 * unified API to:
 *  - test the existence of parameter blocks,
 *  - log blocks or entire parameter sets,
 *  - retrieve the list of available blocks for a given model.
 *
 * The goal of this interface is to decouple *clients* from concrete
 * implementations such as:
 *  - @ref BlockProvider (low-level access via Parameters),
 *  - @ref APIAdapter (monitoring-layer access).
 *
 * Typical implementations (e.g. @ref BlockProxy) act as façades combining
 * multiple internal services.
 *
 * @see BlockProxy
 * @see BlockProvider
 * @see APIAdapter
 * @see Parameters
 */

/**
 * @class IBlockProxy
 * @ingroup DataMonitoringModule
 * @brief Abstract interface for block inspection and logging.
 *
 * @ref IBlockProxy defines the contract for classes that expose block-level
 * diagnostics for Hyperiso models. It provides:
 *  - existence checks for blocks,
 *  - logging of blocks or full parameter sets,
 *  - retrieval of block name lists.
 *
 * This interface is intended to be used by higher-level components that need
 * read-only introspection and diagnostics without direct access to
 * ::Parameters or ::MemoryManager.
 */
class IBlockProxy {
public:
    virtual ~IBlockProxy() = default;

    /**
     * @brief Checks whether a block exists for a given ParameterType.
     *
     * @param blockname Name of the block to check.
     * @param pt ParameterType to query.
     * @return true if the block exists, false otherwise.
     */
    virtual bool exists(const std::string& blockname, ParameterType pt) = 0;

    /**
     * @brief Logs all blocks for a given ParameterType.
     *
     * Implementations typically delegate to a lower-level logging facility
     * (e.g. ::BlockProvider or ::Parameters), which prints the full content
     * of the associated block accessor.
     *
     * @param pt ParameterType whose blocks should be logged.
     */
    virtual void log_all_blocks(ParameterType pt) = 0;

    /**
     * @brief Logs the content of a single block for a given ParameterType.
     *
     * Implementations are expected to rely on the underlying block storage
     * (e.g. ::Parameters::print_block), which streams the block content
     * using the `operator<<` overload for `std::shared_ptr<BlockAccessor>`.
     *
     * @param pt ParameterType owning the block.
     * @param blockname Name of the block to log.
     */
    virtual void log_block(ParameterType pt, const std::string& blockname) = 0;

    /**
     * @brief Retrieves the list of block names for a given ParameterType.
     *
     * @param pt ParameterType to query.
     * @return Unordered set of block names.
     */
    virtual std::unordered_set<BlockName> get_block_list(ParameterType pt) = 0;

    /**
     * @brief Retrieves the content of the block names for a given ParameterType.
     *
     * @param pt ParameterType to query.
     * @param blockname Name of the block.
     * @return map of block content.
     */
    virtual std::map<LhaID, scalar_t> get_block(ParameterType pt, const std::string& blockname) = 0;
};

#endif // IBLOCK_PROXY_H
