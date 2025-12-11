#ifndef ICORE_API_H
#define ICORE_API_H

/**
 * @file ICoreAPI.h
 * @brief Declares a minimal interface for core application state access.
 *
 * This interface is used to expose lightweight read access to core
 * framework state (e.g. current ::Model) to higher-level components
 * without coupling them to concrete implementations.
 */

/**
 * @class ICoreAPI
 * @brief Generic interface for querying core state values.
 *
 * @tparam T Type of the core value to access (e.g. ::Model).
 *
 * Implementations provide a single method ::get() that returns the
 * current value. It is intentionally simple to avoid tight coupling
 * between layers.
 */
template <typename T>
class ICoreAPI{
public:
    virtual ~ICoreAPI() = default;
    
    /**
     * @brief Returns the current state value.
     *
     * @return The current value of type @p T.
     */
    inline virtual T get() = 0;
};

#endif