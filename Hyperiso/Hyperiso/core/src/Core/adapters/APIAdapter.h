#ifndef APIADAPTER_H
#define APIADAPTER_H

#include <unordered_set>
#include <map>

#include "Include.h"
#include "IMonitor.h"
#include "IPathProvider.h"
#include "IDataMonitor.h"
#include "MemoryManager.h"
#include "Parameters.h"

/**
 * @file APIAdapter.h
 * @brief High-level monitoring and access adapter for parameter and path information.
 *
 * This header declares APIAdapter, which aggregates:
 *  - flag monitoring via IMonitor<ExternalFlag>,
 *  - path retrieval via IPathProvider<APIPath>,
 *  - parameter inspection via IDataMonitor.
 */

/**
 * @enum APIPath
 * @brief Enumerates API-exposed paths.
 *
 * - LHA_PATH: path to the currently used LHA input (or spectrum) file.
 */
enum class APIPath {
    LHA_PATH,    ///< Path to the active LHA file.
};


/**
 * @class APIAdapter
 * @ingroup DataMonitoringModule
 * @brief Concrete adapter combining monitoring and path providing functionalities.
 *
 * This class provides a unified, high-level API to:
 *  - query external flags (via IMonitor<ExternalFlag>),
 *  - retrieve important paths like the LHA file path (via IPathProvider<APIPath>),
 *  - inspect parameter blocks and values (via IDataMonitor).
 *
 * It is intended for external tools (CLI, GUI, REST API) that need a
 * read-only view over the internal parameter state managed by Hyperiso.
 *
 * Typical usage:
 * @code
 *   APIAdapter api;
 *   auto all_blocks = api.get_all_blocks();
 *   auto sm_mass_blocks = api.get_blocks_list(ParameterType::SM);
 *   auto mass_values = api.get_block_infos("MASS", ParameterType::SM);
 *   auto lha_path = api.get_path(APIPath::LHA_PATH);
 * @endcode
 *
 * @see IDataMonitor
 * @see IMonitor
 * @see IPathProvider
 */
class APIAdapter : public IMonitor<ExternalFlag>, IPathProvider<APIPath>, IDataMonitor {
public:
    /**
     * @brief Checks the status of a specific external flag.
     *
     * Delegates to MemoryManager::getMemoryCache().config.flags.
     *
     * @param flag The flag to check.
     * @return True if the flag is active, false otherwise.
     */
    bool check_flag(ExternalFlag flag);

    /**
     * @brief Retrieves a specific filesystem path (e.g., LHA file path).
     *
     * For APIPath::LHA_PATH, this returns the LHA or spectrum file path stored in
     * MemoryManager::getMemoryCache().lha_path.
     *
     * @param path_name The requested path.
     * @return The filesystem path corresponding to the enum.
     */
    fs::path get_path(APIPath path_name);

    /**
     * @brief Retrieves all block names across all available model types.
     *
     * This may iterate over all ParameterType instances allowed by MemoryManager
     * and collect their block names.
     *
     * @return A set of all block names.
     */
    std::unordered_set<BlockName> get_all_blocks();

    /**
     * @brief Retrieves block names for a specific parameter type.
     *
     * Delegates to Parameters::GetInstance(param_type)->get_blocks_list().
     *
     * @param param_type The parameter type (default: SM).
     * @return A set of block names for the given type.
     */
    std::unordered_set<BlockName> get_blocks_list(ParameterType param_type = ParameterType::SM);

    /**
     * @brief Retrieves all parameter values inside a given block.
     *
     * Delegates to Parameters::GetInstance(param_type)->get_block_infos(block).
     *
     * @param block Name of the block.
     * @param param_type Type of the model (default: SM).
     * @return A map from LHA IDs to parameter values.
     */
    std::map<LhaID, scalar_t> get_block_infos(const BlockName& block, ParameterType param_type = ParameterType::SM);

    /**
     * @brief Retrieves the list of parameter types that own a specific block.
     *
     * Internally relies on ParamRouter::GetType(block) to list all ParameterType
     * values for which the block is registered.
     *
     * @param block The name of the block.
     * @return A vector of parameter types that own this block.
     */
    std::vector<ParameterType> get_type_of_block(const BlockName& block);
};


#endif // APIADAPTER_H
