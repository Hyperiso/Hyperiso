#ifndef ISOURCES_PROVIDER_H
#define ISOURCES_PROVIDER_H

#include <unordered_set>

#include "Include.h"

/**
 * @file ISourcesProvider.h
 * @brief Interface for retrieving ultimate source parameters behind dependencies.
 *
 * This header declares ISourcesProvider, an interface used to query the set of
 * “leaf” parameters from a collection of (possibly dependent) ParamId objects.
 *
 * A leaf source parameter is typically one that does not depend on any other
 * parameter, or is considered a primitive input in the dependency graph.
 */

 /**
 * @class ISourcesProvider
 * @ingroup DependencyManagementModule
 * @brief Interface for retrieving leaf source parameters from a dependency graph.
 *
 * Implementations of this interface provide a way to:
 *  - traverse the parameter dependency graph,
 *  - extract the ultimate “leaf” source parameters that feed a given set
 *    of dependent parameters.
 */
class ISourcesProvider {
public:
    /// Virtual destructor for interface.
    virtual ~ISourcesProvider() = default;

    /**
     * @brief Retrieves all leaf (ultimate) source parameters for the given set of parameters.
     *
     * The implementation is expected to:
     *  - follow all dependency links starting from each ParamId in param_ids,
     *  - collect the parameters that do not depend on anything else (or are
     *    considered primitive / external inputs),
     *  - return the union of those leaf parameters.
     *
     * @param param_ids Set of parameter IDs to start from.
     * @return Set of all leaf source ParamIds.
     */
    virtual std::unordered_set<ParamId> get_all_leaf_sources(const std::unordered_set<ParamId>& param_ids) 
                                      const = 0;

};

#endif // ISOURCES_PROVIDER_H