#ifndef HYPERISO_MASTER_H
#define HYPERISO_MASTER_H

#include <vector>
#include <string>
#include <map>

#include "IMonitor.h"
#include "MemoryManager.h"
#include "ParamBlockLoader.h"
#include "CorrelationAdapter.h"
#include "SpectrumCalculator.h"
#include "DefaultPathsProvider.h"
#include "MartyRuntimeConfig.h"
#include "SoftSusy.h"

/**
 * @file HyperisoMaster.h
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
private:
    std::map<APIPath, fs::path> path_overrides; ///< Validated pre-init path overrides.

    void ensure_memory_manager_created();
    bool should_validate_marty_runtime(const HyperisoConfig& config) const;
    bool validate_marty_runtime_if_needed(const HyperisoConfig& config, const std::string& context) const;

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
    void init(const std::string &lhaFile, HyperisoConfig config);

    /**
     * @brief Registers one additional LHA block prototype before init().
     *
     * Must be called before init() to affect the first LHA parsing.
     */
    void pre_init_add_block(BlockName blockName,
                            size_t itemCount=2,
                            size_t valueIdx=1,
                            int scaleIdx=-1,
                            int rgIdx=-1,
                            int binIdx=-1,
                            bool globalScale=false);

    /**
     * @brief Registers several additional LHA block prototypes before init().
     */
    void pre_init_add_blocks(const std::vector<LhaPrototypeSpec>& prototypes);

    /**
     * @brief Registers an existing MARTY installation before init().
     *
     * The path may point to the MARTY install prefix itself, or to a nearby
     * include/, lib/, MARTY_INSTALL/, install/, marty.h, or libmarty file. The
     * normalized installation is validated immediately and reused by the MARTY
     * code-generation layer instead of the bundled Third_party/MARTY path.
     */
    void pre_init_set_marty_path(const std::string& martyInstallPath);

    /**
     * @brief Registers a SOFTSUSY executable or installation directory before init().
     *
     * The path may point directly to softpoint.x, or to a directory containing
     * softpoint.x, bin/softpoint.x, or src/SOFTSUSY/softpoint.x. This runtime
     * path is used even when Hyperiso was not built with BUILD_WITH_SOFTSUSY.
     */
    void pre_init_set_softsusy_path(const std::string& softsusyPath);

    /**
     * @brief Overrides selected HyperISO filesystem paths before init().
     *
     * Each entry is validated before being installed: default input files must
     * exist and use the .json extension, user input files must exist and use
     * .yaml or .yml, and directory entries must exist as directories.
     * LHA_PATH is intentionally excluded because the active LHA file is provided
     * directly to init() or switch_lha().
     *
     * Calling this after init() is supported only for future reload/switch
     * operations and emits a warning.
     *
     * @param pathOverrides Map from APIPath values to replacement filesystem paths.
     */
    void pre_init_set_paths(const std::map<APIPath, std::string>& pathOverrides);

    /**
     * @brief Sets the writable MARTY generated-code/cache directory before init().
     *
     * The directory is created if it does not exist.
     */
    void pre_init_set_marty_cache_dir(const std::string& cacheDir);

    /**
     * @brief Sets the writable spectrum cache directory before init().
     *
     * The directory is created if it does not exist.
     */
    void pre_init_set_spectrum_cache_dir(const std::string& cacheDir);

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
    void switch_lha(const std::string &lhaFile, HyperisoConfig config);
};

#endif // HYPERISO_MASTER_H
