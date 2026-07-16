#ifndef PARAM_SOURCES_PROVIDER_H
#define PARAM_SOURCES_PROVIDER_H

#include "ISourcesProvider.h"

/**
 * @file ParamSourcesProvider.h
 * @brief Concrete ISourcesProvider implementation using MemoryManager.
 *
 * This header declares ParamSourcesProvider, which queries the
 * MemoryManager to retrieve ultimate source parameters (leaves of
 * the dependency graph) for a given set of ParamId objects.
 */

 /**
 * @class ParamSourcesProvider
 * @ingroup DependencyManagementModule
 * @brief Provides access to leaf source parameters via MemoryManager.
 *
 * ParamSourcesProvider is a thin adapter that forwards the call to
 * MemoryManager::get_all_source_parameters, allowing client code to
 * remain decoupled from MemoryManager directly.
 */
class ParamSourcesProvider : public ISourcesProvider {
public:
    /// Default constructor.
    ParamSourcesProvider() = default;
    
    /**
     * @brief Retrieves all leaf source parameters for the given set of ParamIds.
     *
     * Delegates to MemoryManager::get_all_source_parameters, which performs
     * the actual traversal of the full parameter space (across all
     * ParameterType instances and inputs).
     *
     * @param param_ids Set of parameter IDs to start from.
     * @return Set of leaf source ParamIds.
     */
    std::unordered_set<ParamId> get_all_leaf_sources(const std::unordered_set<ParamId>& param_ids) 
                                      const override;

};

#endif // PARAM_SOURCES_PROVIDER_H