#ifndef ISTAT_SOURCES_PROXY_H
#define ISTAT_SOURCES_PROXY_H

#include <unordered_set>
#include "Include.h"

/**
 * @file IStatSourcesProxy.h
 * @brief Statistical-layer port for retrieving leaf parameter sources.
 *
 * This interface exposes a statistics-oriented view over the parameter
 * dependency graph. Its main responsibility is to resolve a set of input
 * parameters down to their ultimate “leaf” sources, i.e. the parameters
 * that do not themselves depend on other parameters or dependent blocks.
 *
 * It is intended for the statistics / sensitivity / dependency-analysis layer,
 * where one often wants to answer questions such as:
 * - “Which fundamental inputs does this observable depend on?”
 * - “What are the root nuisance or model parameters behind this derived quantity?”
 *
 * Concrete implementations typically delegate this work to a lower-level
 * provider that walks:
 * - parameter-level dependencies,
 * - block-level dependencies,
 * until only leaf parameters remain.
 *
 * Typical usage:
 * @code
 *   std::shared_ptr<IStatSourcesProxy> proxy = std::make_shared<StatParamSourcesProxy>();
 *
 *   std::unordered_set<ParamId> seeds = {
 *       ParamId{ParameterType::OBSERVABLE, "FOBS", some_id}
 *   };
 *
 *   auto leaves = proxy->get_all_leaf_sources(seeds);
 * @endcode
 *
 * @see StatParamSourcesProxy
 * @see ParamSourcesProvider
 * @see BlockAccessor::get_all_source_parameters
 */

/**
 * @class IStatSourcesProxy
 * @brief Abstract statistics-layer interface for dependency source resolution.
 *
 * This interface provides a single high-level operation:
 * given a set of parameter identifiers, return the complete set of leaf
 * source parameters that ultimately feed them.
 *
 * “Leaf” means that no additional parameter/block dependency is known below
 * that node in the dependency graph.
 */
class IStatSourcesProxy {
public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~IStatSourcesProxy() = default;

    /**
     * @brief Resolves a set of parameters to their leaf dependency sources.
     *
     * Implementations are expected to recursively traverse parameter and/or
     * block dependencies until only root parameters remain.
     *
     * @param param_ids Set of starting parameter identifiers.
     * @return Set of leaf/root parameter identifiers contributing to @p param_ids.
     */
    virtual std::unordered_set<ParamId> get_all_leaf_sources(const std::unordered_set<ParamId>& param_ids) 
                                      const = 0;

};

#endif