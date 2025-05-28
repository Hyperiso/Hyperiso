#ifndef APIADAPTER_H
#define APIADAPTER_H

#include <unordered_set>
#include "Include.h"
#include "IMonitor.h"
#include "IPathProvider.h"
#include "IDataMonitor.h"
#include "MemoryManager.h"
#include "Parameters.h"

enum class APIPath {
    LHA_PATH,
};


/**
 * @class APIAdapter
 * @ingroup DataMonitoringModule
 * @brief Concrete adapter class combining monitoring and path providing functionalities.
 *
 * Provides unified access to parameter block lists, parameter values, flags, and important paths like the LHA file path.
 */
class APIAdapter : public IMonitor<ExternalFlag>, IPathProvider<APIPath>, IDataMonitor {
public:

    /**
     * @brief Checks the status of a specific external flag.
     * @param flag The flag to check.
     * @return True if the flag is active, false otherwise.
     */
    bool check_flag(ExternalFlag flag);

    /**
     * @brief Retrieves a specific filesystem path (e.g., LHA file path).
     * @param path_name The requested path.
     * @return The filesystem path.
     */
    fs::path get_path(APIPath path_name);

    /**
     * @brief Retrieves all block names across all available model types.
     * @return A set of all block names.
     */
    std::unordered_set<BlockName> get_all_blocks();

    /**
     * @brief Retrieves block names for a specific parameter type.
     * @param param_type The parameter type (default: SM).
     * @return A set of block names.
     */
    std::unordered_set<BlockName> get_blocks_list(ParameterType param_type = ParameterType::SM);

    /**
     * @brief Retrieves all parameter values inside a given block.
     * @param block Name of the block.
     * @param param_type Type of the model (default: SM).
     * @return A map from LHA IDs to parameter values.
     */
    std::map<LhaID, scalar_t> get_block_infos(const BlockName& block, ParameterType param_type = ParameterType::SM);

    /**
     * @brief Retrieves the list of parameter types that own a specific block.
     * @param block The name of the block.
     * @return A vector of parameter types.
     */
    std::vector<ParameterType> get_type_of_block(const BlockName& block);
};


#endif // APIADAPTER_H
