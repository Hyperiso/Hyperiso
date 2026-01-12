#ifndef MARTYADAPTER_H
#define MARTYADAPTER_H

#include "IMonitor.h"
#include "IPathProvider.h"
#include "MemoryManager.h"
#include "DefaultPathsProvider.h"

/**
 * @file MartyAdapter.h
 * @brief Adapter giving access to MARTY-related configuration and internal flags.
 *
 * This header declares MartyAdapter, a small helper that:
 *  - monitors internal flags via IMonitor<InternalFlag>,
 *  - exposes MARTY-related filesystem paths via IPathProvider<MartyPath>,
 *  - provides access to the configured MARTY model name.
 */

 /**
 * @enum MartyPath
 * @brief Enumerates MARTY-related filesystem resources.
 *
 * - MODEL_FILE: path to the MARTY model header file (e.g. `<Model>.h`).
 * - TEMPLATE_DIR: path to a directory of MARTY templates.
 * - PARAM_MAPPING_DIR: path to parameter-mapping files used by MARTY.
 */
enum class MartyPath {
    MODEL_FILE,         ///< Path to the MARTY model header file.
    TEMPLATE_DIR,       ///< Path to the directory containing MARTY templates.
    PARAM_MAPPING_DIR,  ///< Path to the directory containing parameter mappings.
};

/**
 * @class MartyAdapter
 * @ingroup MonitoringSystemModule
 * @brief Adapter providing access to MARTY-specific configuration paths and flags.
 *
 * Responsibilities:
 *  - Provide filesystem paths to various MARTY resources using IPathProvider<MartyPath>.
 *  - Expose the state of internal flags via IMonitor<InternalFlag>.
 *  - Give access to the configured MARTY model name, if any.
 *
 * Typical use:
 * @code
 *   MartyAdapter m;
 *   auto model_path = m.get_path(MartyPath::MODEL_FILE);
 *   bool params_changed = m.check_flag(InternalFlag::PARAMS_CHANGED);
 *   std::string model_name = m.get_marty_model_name();
 * @endcode
 *
 * Internally, it uses MemoryManager::getMemoryCache() and the configured IPathsProvider.
 */
class MartyAdapter : public IMonitor<InternalFlag>, IPathProvider<MartyPath> {
public:

    /**
     * @brief Retrieves the path to a specific MARTY-related resource.
     *
     * The mapping between MartyPath and actual filesystem paths is typically
     * provided by an implementation of IPathsProvider (e.g., DefaultPathsProvider),
     * queried through MemoryManager.
     *
     * @param path_name The MARTY resource to retrieve.
     * @return Filesystem path to the requested resource.
     */
    fs::path get_path(MartyPath path_name) override;

    /**
     * @brief Checks the status of a specific internal flag.
     *
     * Delegates to MemoryManager::getMemoryCache().flags.
     *
     * @param flag The flag to check.
     * @return True if active, false otherwise.
     */
    bool check_flag(InternalFlag flag) override;

    /**
     * @brief Retrieves the name of the MARTY model currently configured.
     *
     * Typically comes from HyperisoConfig::mty_model_name stored in MemoryManager's cache.
     *
     * @return The name of the MARTY model, or an empty string if not set.
     */
    std::string get_marty_model_name() const;
};


#endif // MARTYADAPTER_H
