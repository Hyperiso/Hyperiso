#ifndef IMONITOR_H
#define IMONITOR_H

/**
 * @example examples/monitoring_system_example.cpp
 * @brief Example of using monitoring and access classes like HyperisoMaster, MartyAdapter, APIAdapter.
 * @defgroup MonitoringSystemModule Monitoring and Access System
 * @brief Provides unified access to internal framework data such as flags, paths, parameters and models.
 */

 /**
 * @class IMonitor
 * @ingroup MonitoringSystemModule
 * @brief Interface for monitoring flags inside the framework (internal or external).
 */
template<typename FlagType>
class IMonitor {
public:
    virtual ~IMonitor() = default;

    /**
     * @brief Checks if a specific flag is active.
     * @param flag The flag to check.
     * @return True if active, false otherwise.
     */
    virtual bool check_flag(FlagType flag) = 0;
};

#endif // __IMONITOR_H__
