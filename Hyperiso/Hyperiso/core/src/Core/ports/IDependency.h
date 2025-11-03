#ifndef IDEPENDECY_H
#define IDEPENDECY_H

#include "Include.h"
#include "Block.h"
#include "Parameter.h"
#include "DependentParameter.h"

/**
 * @example dependency_management_example.cpp
 * 
 * @defgroup DependencyManagementModule Dependency Management System
 * @brief Provides interfaces and classes to manage dependencies between blocks and parameters.
 *
 * ## Related Classes
 * - @ref IDependency
 * - @ref CompositeParamAdapter
 */

 /**
 * @class IDependency
 * @ingroup DependencyManagementModule
 * @brief Interface for managing dependencies between parameter blocks or parameters.
 *
 * Provides an abstraction layer for adding, removing, and updating dependencies.
 */
class IDependency {
public:
    virtual ~IDependency() = default;

    /**
     * @brief Adds a dependent block.
     * @param name Name of the dependent block.
     * @param source_names Map of source blocks per parameter type.
     * @param dest Destination parameter type where the new block will be stored.
     * @param recalculateFunc Function used to recalculate the dependent block.
     */
    virtual void add_block_dependency(const BlockName& name,
                                      const std::unordered_map<ParameterType, std::vector<std::string>>& source_names,
                                      ParameterType dest,
                                      DepUpdateFunc recalculateFunc) 
                                      = 0;

    /**
     * @brief Adds a dependent parameter.
     * @param pid ID of the new dependent parameter.
     * @param source_pids Set of source parameter IDs the new parameter depends on.
     * @param recalculateFunc Function used to recalculate the dependent parameter.
     */
    virtual void add_param_dependency(const ParamId& pid,
                                      const std::unordered_set<ParamId>& source_pids,
                                      DepParamUpdateFunc recalculateFunc) 
                                      = 0;

    /**
     * @brief Removes a dependency by name and parameter type.
     * @param name Name of the dependency block.
     * @param src Source parameter type.
     */
    virtual void remove_dependency(const BlockName& name, ParameterType src) = 0;

    /**
     * @brief Updates a dependency by name and parameter type.
     * @param name Name of the dependency block.
     * @param src Source parameter type.
     */
    virtual void update_dependency(const BlockName& name, ParameterType src) = 0;
};

#endif // IDEPENDECY_H
