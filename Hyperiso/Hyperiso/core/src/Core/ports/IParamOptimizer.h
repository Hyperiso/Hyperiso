#ifndef IPARAM_OPTIMIZER_H
#define IPARAM_OPTIMIZER_H

#include "Include.h"
#include "Parameter.h"

/**
 * @file IParamOptimizer.h
 * @brief Interface for batched parameter updates and optimization.
 *
 * This header declares IParamOptimizer, an abstract interface intended
 * for classes that perform batched modifications of parameter values.
 *
 * Typical implementations:
 *  - accumulate a list of operations (set/remove),
 *  - apply them atomically or in a controlled way via commit(),
 *  - offer a clear() method to reset the internal operation queue.
 */

 /**
 * @defgroup ParamOptimizationModule Parameter Optimization System
 * @brief Interfaces and adapters for batched parameter updates.
 *
 * This module provides:
 *  - IParamOptimizer: abstract interface to define a parameter optimization API.
 *  - ParamOptimizer: concrete engine working directly on BlockAccessor instances.
 *  - ParamOptimizerAdapter: high-level adapter constructing ParamOptimizer from
 *    ParameterType scopes.
 */

/**
 * @class IParamOptimizer
 * @ingroup ParamOptimizationModule
 * @brief Interface for batched parameter modifications.
 *
 * IParamOptimizer defines an abstract API for:
 *  - staging updates of individual parameters (by value or by full Parameter),
 *  - removing parameters,
 *  - committing staged operations (optionally coalescing multiple changes),
 *  - clearing the queued operations without applying them.
 */
class IParamOptimizer {
public:
    /// Virtual destructor for interface.
    virtual ~IParamOptimizer() = default;
    
    /**
     * @brief Stage a change of central value for a parameter.
     *
     * @param block Block name where the parameter resides.
     * @param id    LhaID of the parameter.
     * @param v     New value to set.
     */
    virtual void set_value(const BlockName& block, const LhaID& id, scalar_t v) = 0;

    /**
     * @brief Stage a full replacement of a parameter object.
     *
     * @param block Block name where the parameter resides.
     * @param id    LhaID of the parameter.
     * @param p     Shared pointer to the new Parameter.
     */
    virtual void set_param(const BlockName& block, const LhaID& id, std::shared_ptr<Parameter> p)  = 0;

    /**
     * @brief Stage the removal of a parameter.
     *
     * @param block Block name where the parameter resides.
     * @param id    LhaID of the parameter to remove.
     */
    virtual void remove(const BlockName& block, const LhaID& id)  = 0;

    /**
     * @brief Commit all staged operations.
     *
     * Implementations may:
     *  - coalesce multiple operations on the same (block, id) into a single
     *    effective update,
     *  - freeze/unfreeze affected blocks during the update,
     *  - notify observers only once per block.
     *
     * @param coalesce If true, coalesce multiple operations on the same parameter.
     */
    virtual void commit(bool coalesce = true)  = 0;

    /**
     * @brief Clears all currently staged operations without applying them.
     */
    virtual void clear()  = 0;

};

#endif