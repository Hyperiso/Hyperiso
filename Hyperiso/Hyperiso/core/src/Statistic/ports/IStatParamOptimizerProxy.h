#ifndef ISTAT_PARAM_OPTIMIZER_PROXY_H
#define ISTAT_PARAM_OPTIMIZER_PROXY_H

#include "Include.h"
#include "Parameter.h"

/**
 * @file IStatParamOptimizerProxy.h
 * @brief Statistical-layer port for batched parameter edits and optimized commits.
 *
 * This interface exposes a small write-oriented API used by the statistics
 * module to stage parameter modifications before applying them to the live
 * parameter graph.
 *
 * The intended workflow is:
 * - enqueue updates with @ref set_value() or @ref set_param(),
 * - optionally enqueue removals with @ref remove(),
 * - apply all pending changes with @ref commit(),
 * - or discard/reset the staging area with @ref clear().
 *
 * This is useful when many parameter updates must be applied together,
 * for example in scans, fits, profiling, or repeated likelihood evaluations,
 * where immediate propagation after each elementary change would be too costly.
 *
 * Concrete implementations typically forward to an optimizer / transaction-like
 * adapter that knows how to:
 * - batch updates,
 * - reduce redundant edits,
 * - and apply them in a dependency-friendly order.
 *
 * @see StatParamOptimizerProxy
 * @see ParamOptimizerAdapter
 */

/**
 * @class IStatParamOptimizerProxy
 * @brief Abstract interface for batched parameter mutation in the statistics layer.
 *
 * This interface is intentionally minimal and focuses on staging operations
 * over parameters identified by:
 * - a block name,
 * - an LHA-like identifier,
 * - and either a raw scalar value or a full @ref Parameter object.
 *
 * Changes are not required to take effect immediately: implementations may
 * defer them until @ref commit() is called.
 */
class IStatParamOptimizerProxy {
public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~IStatParamOptimizerProxy() = default;
    
    /**
     * @brief Stages a scalar value update for a parameter.
     *
     * The parameter identified by (@p block, @p id) will later be assigned
     * the value @p v when @ref commit() is called.
     *
     * @param block Block containing the parameter.
     * @param id    LHA-like identifier inside the block.
     * @param v     New scalar value to assign.
     */
    virtual void set_value(const BlockName& block, const LhaID& id, scalar_t v) = 0;

    /**
     * @brief Stages a full parameter replacement/update.
     *
     * The parameter identified by (@p block, @p id) will later be updated
     * from the provided @ref Parameter object when @ref commit() is called.
     *
     * @param block Block containing the parameter.
     * @param id    LHA-like identifier inside the block.
     * @param p     Shared pointer to the parameter payload to apply.
     */
    virtual void set_param(const BlockName& block, const LhaID& id, std::shared_ptr<Parameter> p)  = 0;

    /**
     * @brief Stages the removal of a parameter.
     *
     * The parameter identified by (@p block, @p id) will be removed on the
     * next @ref commit().
     *
     * @param block Block containing the parameter.
     * @param id    LHA-like identifier of the parameter to remove.
     */
    virtual void remove(const BlockName& block, const LhaID& id)  = 0;

    /**
     * @brief Applies all staged operations.
     *
     * Implementations may use @p coalesce to merge or compact multiple pending
     * edits before applying them. The exact coalescing behavior is adapter-defined.
     *
     * @param coalesce If true, allows the implementation to merge/reduce pending operations.
     */
    virtual void commit(bool coalesce = true)  = 0;

    /**
     * @brief Clears all staged operations without applying them.
     */
    virtual void clear()  = 0;

};

#endif