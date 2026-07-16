#ifndef STAT_PARAM_SOURCES_PROXY_H
#define STAT_PARAM_SOURCES_PROXY_H

#include "IStatSourcesProxy.h"
#include "ParamSourcesProvider.h"

/**
 * @file StatParamSourcesProxy.h
 * @brief Statistics-layer adapter for retrieving leaf parameter sources.
 *
 * This class is a thin proxy over @ref ParamSourcesProvider.
 * It adapts the dependency-source traversal service to the statistics module
 * through the @ref IStatSourcesProxy interface.
 *
 * The proxy is intentionally lightweight:
 * - it owns a provider instance,
 * - it forwards leaf-source queries,
 * - it exposes a stable interface to higher-level statistics code.
 *
 * @see IStatSourcesProxy
 * @see ParamSourcesProvider
 */

/**
 * @class StatParamSourcesProxy
 * @brief Adapter exposing parameter leaf-source resolution to the statistics layer.
 *
 * This class delegates all work to an internal @ref ParamSourcesProvider.
 * It is used when the statistics module needs to identify the ultimate input
 * parameters behind one or more derived quantities.
 */
class StatParamSourcesProxy : public IStatSourcesProxy {
public:
    /**
     * @brief Default constructor.
     *
     * Constructs the proxy with its internal @ref ParamSourcesProvider.
     */
    StatParamSourcesProxy() = default;
    
    /**
     * @brief Returns all leaf/root sources contributing to a set of parameters.
     *
     * Delegates to:
     * @code
     * psp.get_all_leaf_sources(param_ids);
     * @endcode
     *
     * @param param_ids Set of starting parameter identifiers.
     * @return Set of leaf/root parameter identifiers.
     */
    std::unordered_set<ParamId> get_all_leaf_sources(const std::unordered_set<ParamId>& param_ids) 
                                      const override;
private:
    ParamSourcesProvider psp;   /// Underlying provider performing the actual dependency traversal.
};

#endif