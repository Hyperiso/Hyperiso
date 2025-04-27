#ifndef MARTYADAPTER_H
#define MARTYADAPTER_H

#include "IMonitor.h"
#include "IPathProvider.h"
#include "MemoryManager.h"


enum class MartyPath {
    MODEL_FILE,
    TEMPLATE_DIR,
    PARAM_MAPPING_DIR,
};

/**
 * @class MartyAdapter
 * @ingroup MonitoringSystemModule
 * @brief Adapter providing access to MARTY-specific configuration paths and flags.
 */
class MartyAdapter : public IMonitor<InternalFlag>, IPathProvider<MartyPath> {
private:
    /* data */
public:

    /**
     * @brief Retrieves the path to a specific MARTY-related resource.
     * @param path_name The resource to retrieve.
     * @return Filesystem path to the requested resource.
     */
    fs::path get_path(MartyPath path_name) override;

    /**
     * @brief Checks the status of a specific internal flag.
     * @param flag The flag to check.
     * @return True if active, false otherwise.
     */
    bool check_flag(InternalFlag flag) override;

    /**
     * @brief Retrieves the name of the MARTY model currently configured.
     * @return The name of the MARTY model.
     */
    std::string get_marty_model_name();
};


#endif // MARTYADAPTER_H
