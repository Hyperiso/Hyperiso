#ifndef IMONITOR_H
#define IMONITOR_H

/**
 * @file IMonitor.h
 * @brief Generic interface for monitoring flags within the framework.
 *
 * This header defines the templated IMonitor interface used by:
 *  - @ref HyperisoMaster (monitoring external flags, see ExternalFlag)
 *  - @ref MartyAdapter   (monitoring internal flags, see InternalFlag)
 */

/**
 * @example examples/monitoring_system_example.cpp
 * @brief Example of using monitoring and access classes like HyperisoMaster, MartyAdapter, APIAdapter.
 *
 * @defgroup MonitoringSystemModule Monitoring and Access System
 * @brief Provides unified access to internal framework data such as flags, paths, parameters and models.
 */

/**
 * @class IMonitor
 * @ingroup MonitoringSystemModule
 * @brief Interface for monitoring flags inside the framework (internal or external).
 *
 * @tparam FlagType Enum-like type describing the set of flags (e.g. ExternalFlag, InternalFlag).
 *
 * Implementations provide a type-safe way to query whether a particular flag is set
 * in the underlying system (e.g. MemoryManager cache, runtime configuration).
 */
template<typename FlagType>
class IMonitor {
public:
    /// Virtual destructor for proper polymorphic use.
    virtual ~IMonitor() = default;

    /**
     * @brief Checks if a specific flag is active.
     *
     * Concrete implementations should route this query to the actual storage of
     * flags (for example MemoryManager::getMemoryCache().flags or config.flags).
     *
     * @param flag The flag to check.
     * @return True if active, false otherwise.
     */
    virtual bool check_flag(FlagType flag) = 0;
};

#endif // IMONITOR_H
