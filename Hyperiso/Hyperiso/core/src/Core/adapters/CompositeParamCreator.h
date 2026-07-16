#ifndef COMPOSITEPARAMCREATOR_H
#define COMPOSITEPARAMCREATOR_H

#include "IDependency.h"
#include "DependentBlockManager.h"
#include "DependentParameter.h"

/**
 * @file CompositeParamAdapter.h
 * @brief Adapter that delegates dependency management to DependentBlockManager.
 *
 * This file defines the CompositeParamAdapter class, a concrete implementation of
 * the IDependency interface. It provides a thin abstraction layer that forwards
 * all dependency operations (creation, removal, update) to the underlying
 * DependentBlockManager utilities.
 *
 * @ingroup DependencyManagementModule
 */

/**
 * @class CompositeParamAdapter
 * @ingroup DependencyManagementModule
 * @brief Concrete implementation of IDependency using DependentBlockManager.
 *
 * This adapter glues high-level dependency management code (e.g. configuration,
 * model setup) to the underlying mechanisms based on:
 * - DependentBlock (for block-level dependencies), and
 * - DependentParameter (for parameter-level dependencies).
 *
 * All methods simply forward to the static helpers in DependentBlockManager:
 * - add_block_dependency() → DependentBlockManager::addDependentBlock()
 * - add_param_dependency() → DependentBlockManager::addDependentParameter()
 * - remove_dependency()    → DependentBlockManager::removeDependentBlock()
 * - update_dependency()    → DependentBlockManager::update()
 *
 * This indirection allows you to:
 * - mock or replace the dependency manager in tests,
 * - centralize all dependency creation logic behind a single interface (IDependency).
 *
 * @see IDependency
 * @see DependentBlockManager
 * @see DependentBlock
 * @see DependentParameter
 */
class CompositeParamAdapter : public IDependency {
public:
    /**
     * @brief Adds a dependent block using DependentBlockManager.
     *
     * Creates a new DependentBlock whose content is computed from a set of
     * source blocks, possibly coming from different Parameters instances.
     *
     * Internally, this calls:
     * @code
     * DependentBlockManager::addDependentBlock(name, source_names, dest, recalculateFunc);
     * @endcode
     *
     * @param name Name of the dependent block to create.
     * @param source_names Map associating each source ParameterType to the list
     *                     of block names that will serve as inputs.
     * @param dest Destination ParameterType in which the new block will be registered.
     * @param recalculateFunc Lambda used to recompute the block from its sources.
     */
    void add_block_dependency(
        const BlockName& name,
        const std::unordered_map<ParameterType, std::vector<std::string>>& source_names,
        ParameterType dest,
        DepUpdateFunc recalculateFunc
    ) override;

    /**
     * @brief Adds a dependent parameter using DependentBlockManager.
     *
     * Creates a DependentParameter whose value is computed from one or several
     * source parameters, possibly across different ParameterType instances.
     *
     * Internally, this calls:
     * @code
     * DependentBlockManager::addDependentParameter(pid, source_pids, recalculateFunc);
     * @endcode
     *
     * @param pid ID of the dependent parameter to create.
     * @param source_pids Set of IDs for the source parameters.
     * @param recalculateFunc Lambda used to recompute the parameter value.
     */
    void add_param_dependency(
        const ParamId& pid,
        const std::unordered_set<ParamId>& source_pids,
        DepParamUpdateFunc recalculateFunc
    ) override;

    /**
     * @brief Removes a dependent block using DependentBlockManager.
     *
     * Internally, this calls:
     * @code
     * DependentBlockManager::removeDependentBlock(name, src);
     * @endcode
     *
     * @param name Name of the dependent block to remove.
     * @param src  ParameterType (model) where the block is registered.
     */
    void remove_dependency(const BlockName& name, ParameterType src) override;

    /**
     * @brief Manually triggers an update on a dependent block using DependentBlockManager.
     *
     * Internally, this calls:
     * @code
     * DependentBlockManager::update(name, src);
     * @endcode
     *
     * @param name Name of the dependent block to update.
     * @param src  ParameterType (model) where the block is located.
     */
    void update_dependency(const BlockName& name, ParameterType src) override;

};

#endif // COMPOSITEPARAMCREATOR_H
