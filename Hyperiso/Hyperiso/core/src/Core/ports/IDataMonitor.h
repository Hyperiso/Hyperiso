#ifndef IDATAMONITOR_H
#define IDATAMONITOR_H

#include <unordered_set>
#include <map>
#include <vector>

#include "Include.h"
#include "scalar.h"

/**
 * @file IDataMonitor.h
 * @brief Interface for monitoring and retrieving information about parameter blocks.
 *
 * This header defines IDataMonitor, which provides a read-only, structured API
 * to inspect:
 *  - the set of available blocks,
 *  - the content of specific blocks,
 *  - the types owning a given block.
 */

/**
 * @example data_monitor_example.cpp
 * @brief Example usage of APIAdapter to inspect all parameter blocks.
 *
 * @defgroup DataMonitoringModule Data Monitoring System
 * @brief Provides interfaces and classes to monitor, retrieve, and inspect parameter data.
 *
 * This module defines:
 *  - Interfaces for monitoring parameter blocks and flags.
 *  - Concrete implementations to provide unified access to block lists, parameter types, and paths.
 *
 * ## Related Classes
 * - @ref IDataMonitor
 * - @ref APIAdapter
 */

/**
 * @class IDataMonitor
 * @ingroup DataMonitoringModule
 * @brief Interface for monitoring and retrieving information about parameter blocks.
 *
 * Implementations provide read-only access to:
 *  - the list of all blocks across models,
 *  - per-type block lists,
 *  - the map of parameter values inside a block,
 *  - the mapping between blocks and owning ParameterTypes.
 */
class IDataMonitor {
public:
    /// Virtual destructor for polymorphic use.
    virtual ~IDataMonitor() = default;

    /**
     * @brief Retrieves the names of all parameter blocks across all models.
     *
     * Implementations may aggregate block lists from multiple ParameterType
     * instances as well as additional sources (e.g. input cache).
     *
     * @return A set of block names.
     */
    virtual std::unordered_set<BlockName> get_all_blocks() = 0;

    /**
     * @brief Retrieves the list of blocks for a specific model type.
     *
     * Typically delegates to Parameters::GetInstance(param_type)->get_blocks_list().
     *
     * @param param_type The type of model (default: SM).
     * @return A set of block names.
     */
    virtual std::unordered_set<BlockName> get_blocks_list(ParameterType param_type = ParameterType::SM) = 0;

    /**
     * @brief Retrieves all parameter values within a block.
     *
     * @param block Name of the block.
     * @param param_type Model type associated with the block (default: SM).
     * @return A map of LHA IDs to their corresponding values.
     */
    virtual std::map<LhaID, scalar_t> get_block_infos(const BlockName& block, ParameterType param_type = ParameterType::SM) = 0;


    /**
     * @brief Retrieves the model types owning a specific block.
     *
     * For a given block, a model can appear multiple times in different
     * contexts (SM, BSM, WILSON, etc.). Implementations generally rely on
     * ParamRouter::GetType(block) to determine the corresponding types.
     *
     * @param block Name of the block.
     * @return A vector of model types associated with the block.
     */
    virtual std::vector<ParameterType> get_type_of_block(const BlockName& block) = 0;
};


#endif // IDATAMONITOR_H
