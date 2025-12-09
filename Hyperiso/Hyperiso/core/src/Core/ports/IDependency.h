#ifndef IDEPENDECY_H
#define IDEPENDECY_H

#include "Include.h"
#include "Block.h"
#include "Parameter.h"
#include "DependentParameter.h"

/**
 * @file IDependency.h
 * @brief Interface for managing dependencies between blocks and parameters.
 *
 * This header defines the IDependency interface, which abstracts the creation,
 * update, and removal of:
 *  - dependent blocks (based on one or more source blocks), and
 *  - dependent parameters (based on one or more source parameters).
 *
 * Implementations typically delegate to concrete mechanisms such as
 * DependentBlock / DependentParameter and their managers.
 */

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
 * Typical implementations (such as CompositeParamAdapter) internally rely on
 * DependentBlock and DependentParameter, and orchestrate them via helper
 * managers (e.g. DependentBlockManager).
 */
class IDependency {
public:
    /// Virtual destructor for proper polymorphic cleanup.
    virtual ~IDependency() = default;

    /**
     * @brief Adds a dependent block.
     *
     * Creates a new logical block whose contents are derived from other blocks.
     *
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
     *
     * Registers a new parameter whose value is computed from one or more
     * source parameters.
     *
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
     *
     * Usually targets a dependent block that should be removed (and unlinked
     * from its sources/observers).
     *
     * @param name Name of the dependency block.
     * @param src Source parameter type.
     */
    virtual void remove_dependency(const BlockName& name, ParameterType src) = 0;

    /**
     * @brief Updates a dependency by name and parameter type.
     *
     * Typically forces a recomputation of a dependent block.
     *
     * @param name Name of the dependency block.
     * @param src Source parameter type.
     */
    virtual void update_dependency(const BlockName& name, ParameterType src) = 0;
};

#endif // IDEPENDECY_H
