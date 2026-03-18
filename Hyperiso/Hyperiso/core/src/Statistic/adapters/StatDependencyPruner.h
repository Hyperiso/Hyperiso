#ifndef STATDEPENDENCYPRUNER_H
#define STATDEPENDENCYPRUNER_H

#include "IStatDependencyPruner.h"
#include "DependencyPruner.h"

/**
 * @file StatDependencyPruner.h
 * @brief Statistics-layer adapter over the core @ref DependencyPruner service.
 *
 * This class adapts the core dependency-pruning API to the needs of the
 * statistics module.
 *
 * Compared to @ref DependencyPruner:
 * - it exposes block identifiers as plain `std::string`,
 * - it keeps the same high-level detach / reattach semantics,
 * - it delegates all actual work to an internal @ref DependencyPruner instance.
 *
 * This allows the statistics module to interact with the dependency graph
 * without depending directly on the core block abstraction.
 *
 * @see IStatDependencyPruner
 * @see DependencyPruner
 * @see Parameters
 */

/**
 * @class StatDependencyPruner
 * @brief Adapter that bridges the statistics module to the core dependency-pruning backend.
 *
 * This class is a thin façade over @ref DependencyPruner.
 *
 * It forwards all requests to the underlying core adapter, converting the
 * string block name interface of the statistics layer into the corresponding
 * core calls.
 */
class StatDependencyPruner : public IStatDependencyPruner {
public:
    /**
     * @brief Reattaches a previously detached block.
     *
     * Delegates to:
     * @code
     * dp.reattach_block(tp, block_name);
     * @endcode
     *
     * @param tp         Parameter namespace containing the block.
     * @param block_name Name of the block to reattach.
     */
    void reattach_block(ParameterType tp, const std::string& block_name) override;

    /**
     * @brief Detaches a block from its upstream dependencies.
     *
     * Delegates to:
     * @code
     * dp.detach_block(tp, block_name);
     * @endcode
     *
     * @param tp         Parameter namespace containing the block.
     * @param block_name Name of the block to detach.
     */
    void detach_block(ParameterType tp, const std::string& block_name) override;

    /**
     * @brief Reattaches a previously detached dependent parameter.
     *
     * Delegates to:
     * @code
     * dp.reattach_parameter(tp, block_name, id);
     * @endcode
     *
     * @param tp         Parameter namespace containing the parameter.
     * @param block_name Name of the parent block.
     * @param id         Identifier of the parameter to reattach.
     */
    void reattach_parameter(ParameterType tp, const std::string& block_name, const LhaID& id) override;

    /**
     * @brief Detaches a dependent parameter from its upstream dependencies.
     *
     * Delegates to:
     * @code
     * dp.detach_parameter(tp, block_name, id);
     * @endcode
     *
     * @param tp         Parameter namespace containing the parameter.
     * @param block_name Name of the parent block.
     * @param id         Identifier of the parameter to detach.
     */
    void detach_parameter(ParameterType tp, const std::string& block_name, const LhaID& id) override;

private:
    DependencyPruner dp;    /// Core dependency-pruning adapter used to perform the actual operations.
};

#endif // STATDEPENDENCYPRUNER_H
