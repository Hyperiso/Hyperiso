#ifndef ICORE_API_H
#define ICORE_API_H

/**
 * @file ICoreAPI.h
 * @brief Minimal interface for querying core framework state in a type-safe way.
 *
 * This header defines @ref ICoreAPI, a small generic interface used to expose
 * a single "core" value from the framework (e.g. current Model, a boolean flag, etc.).
 *
 * Implementations typically wrap higher-level façade classes such as @ref HyperisoMaster
 * to avoid exposing low-level singletons (e.g. MemoryManager) directly to clients.
 *
 * @tparam T Type returned by the API (e.g., Model, bool).
 */

/**
 * @class ICoreAPI
 * @brief Generic interface returning a single value of type T from the core.
 *
 * Implementations provide a unified way to query a core state value:
 * @code
 *   std::shared_ptr<ICoreAPI<Model>> api = std::make_shared<ModelAPI>();
 *   Model m = api->get();
 * @endcode
 *
 * @tparam T Type of the returned value.
 */
template <typename T>
class ICoreAPI{
public:
    virtual ~ICoreAPI() = default;

    /**
     * @brief Returns the current value exposed by this API.
     *
     * The semantics depend on the implementation: it may query a singleton-backed
     * cache (via HyperisoMaster / MemoryManager), or return a derived property.
     *
     * @return Current value of type @p T.
     */
    virtual T get() = 0;
};

#endif // ICORE_API_H