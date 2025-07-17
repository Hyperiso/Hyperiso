#ifndef CONFIGPROVIDER_H
#define CONFIGPROVIDER_H

#include "IMonitor.h"
#include "MemoryManager.h"
#include "ParamBlockLoader.h"

/**
 * @class HyperisoMaster
 * @ingroup MonitoringSystemModule
 * @brief High-level interface to initialize and monitor the main framework configuration.
 *
 * Provides initialization routines and allows querying important configuration flags and model type.
 */
class HyperisoMaster : public IMonitor<ExternalFlag> {
public:

    /**
     * @brief Initializes Hyperiso using a LHA file and a full Config object.
     * @param lhaFile Path to the LHA input file.
     * @param config Configuration settings.
     */
    void init(const std::string &lhaFile, Config config);

    /**
     * @brief Initializes Hyperiso with only the LHA file, using default config.
     * @param lhaFile Path to the LHA input file.
     */
    void init(const std::string &lhaFile);

    /**
     * @brief Checks if a specific external flag is active.
     * @param flag The flag to check.
     * @return True if active, false otherwise.
     */
    bool check_flag(ExternalFlag flag);

    /**
     * @brief Retrieves the current model type used.
     * @return The active model.
     */
    Model get_model();

    /**
     * @brief Switch lha in Hyperiso, and use a new configuration.
     * @param lhaFile Path to the LHA input file.
     * @param config Configuration used for hyperiso.
     */
    void switch_lha(const std::string &lhaFile, Config config);
};

#endif // CONFIGPROVIDER_H
