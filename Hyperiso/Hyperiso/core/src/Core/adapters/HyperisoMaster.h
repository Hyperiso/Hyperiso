#ifndef HYPERISO_MASTER_H
#define HYPERISO_MASTER_H

#include "IMonitor.h"
#include "MemoryManager.h"
#include "ParamBlockLoader.h"
#include "CorrelationAdapter.h"
#include "SpectrumCalculator.h"
#include "DefaultPathsProvider.h"

/**
 * @file ConfigProvider.h
 * @brief High-level helpers for initializing and monitoring the Hyperiso framework.
 *
 * This header declares:
 *  - @ref HyperisoMaster: a thin façade over MemoryManager to initialize Hyperiso
 *    and query high-level configuration flags and model information.
 *
 * It relies on:
 *  - MemoryManager and its loaders (ParamBlockLoader, CorrelationLoader)
 *  - SpectrumCalculator for spectrum generation
 *  - DefaultPathsProvider for default asset and configuration paths
 */

/**
 * @class HyperisoMaster
 * @ingroup MonitoringSystemModule
 * @brief High-level interface to initialize and monitor the main framework configuration.
 *
 * This class is a convenience façade over the low-level infrastructure:
 *  - It sets up a MemoryManager with concrete loaders and path providers.
 *  - It initializes Hyperiso from an LHA file and configuration.
 *  - It exposes monitoring/inspection utilities for external flags and the active model.
 *
 * Typical usage:
 * @code
 *   HyperisoMaster master;
 *   Config config;
 *   config.model = Model::SM;
 *   master.init("input.lha", config);
 *
 *   if (master.check_flag(ExternalFlag::IS_LHA_SPECTRUM)) {
 *       // ...
 *   }
 *
 *   Model active_model = master.get_model();
 * @endcode
 *
 * @see MemoryManager
 * @see ParamBlockLoader
 * @see CorrelationLoader
 * @see SpectrumCalculator
 * @see DefaultPathsProvider
 */
class HyperisoMaster : public IMonitor<ExternalFlag> {
public:

    /**
     * @brief Initializes Hyperiso using a LHA file and a full Config object.
     *
     * Internally:
     *  - Creates and wires the default IDataLoader and ISpectrumCalculator instances
     *    (ParamBlockLoader, CorrelationLoader, SpectrumCalculator).
     *  - Creates a DefaultPathsProvider.
     *  - Builds the MemoryManager singleton using these components.
     *  - Calls MemoryManager::init(lhaFile, config).
     *
     * @param lhaFile Path to the LHA input file.
     * @param config Configuration settings to drive initialization.
     */
    void init(const std::string &lhaFile, Config config);

    /**
     * @brief Initializes Hyperiso with only the LHA file, using default config.
     *
     * Equivalent to calling:
     * @code
     *   Config cfg;               // default constructed (SM, no special flags)
     *   master.init(lhaFile, cfg);
     * @endcode
     *
     * @param lhaFile Path to the LHA input file.
     */
    void init(const std::string &lhaFile);

    /**
     * @brief Checks if a specific external flag is active.
     *
     * Delegates to the MemoryManager::getMemoryCache().config.flags map.
     *
     * @param flag The flag to check.
     * @return True if active, false otherwise.
     */
    bool check_flag(ExternalFlag flag);

    /**
     * @brief Retrieves the current model type used by Hyperiso.
     *
     * Delegates to MemoryManager::getMemoryCache().config.model.
     *
     * @return The active model.
     */
    Model get_model();

    /**
     * @brief Switches the LHA file and applies a new configuration.
     *
     * Internally calls MemoryManager::switch_lha(lhaFile, config), which:
     *  - Restores the default cache snapshot,
     *  - Re-reads spectrum / LHA data,
     *  - Updates the internal Config and flags.
     *
     * @param lhaFile Path to the LHA input file.
     * @param config Configuration to use for Hyperiso.
     */
    void switch_lha(const std::string &lhaFile, Config config);
};

#endif // CONFIGPROVIDER_H
