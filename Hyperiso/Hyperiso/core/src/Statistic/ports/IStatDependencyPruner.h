#ifndef ISTATDEPENDENCYPRUNER_H
#define ISTATDEPENDENCYPRUNER_H

#include "Include.h"

/**
 * @file IStatDependencyPruner.h
 * @brief Statistical-layer port for pruning and restoring parameter dependencies.
 *
 * This interface is the statistics-module counterpart of the core dependency
 * pruning API. It exposes a simplified string-based interface to:
 * - detach a dependent block,
 * - reattach a dependent block,
 * - detach a dependent parameter,
 * - reattach a dependent parameter.
 *
 * The statistical layer uses raw block names as `std::string` instead of the
 * richer @ref BlockName abstraction used in the core module. Concrete adapters
 * are expected to bridge this interface to the underlying dependency backend.
 *
 * Typical usage:
 * @code
 *   std::shared_ptr<IStatDependencyPruner> pruner = std::make_shared<StatDependencyPruner>();
 *   pruner->detach_block(ParameterType::SM, "VCKM");
 *   ...
 *   pruner->reattach_block(ParameterType::SM, "VCKM");
 * @endcode
 *
 * @see StatDependencyPruner
 * @see DependencyPruner
 * @see Parameters
 */

/**
 * @class IStatDependencyPruner
 * @brief Abstract port used by the statistics module to prune dependency links.
 *
 * This interface allows statistical workflows to temporarily disable automatic
 * dependency propagation on:
 * - whole blocks,
 * - individual dependent parameters.
 *
 * It is intentionally lightweight and string-based so that the statistics
 * module does not depend directly on the full block-identifier abstraction
 * of the core parameter system.
 */
class IStatDependencyPruner {
public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~IStatDependencyPruner() = default;

    /**
     * @brief Reattaches a previously detached block.
     *
     * @param type       Parameter namespace/domain containing the block.
     * @param block_name Name of the block to reattach.
     */
    virtual void reattach_block(ParameterType, const std::string&) = 0;

    /**
     * @brief Detaches a block from its upstream dependency graph.
     *
     * @param type       Parameter namespace/domain containing the block.
     * @param block_name Name of the block to detach.
     */
    virtual void detach_block(ParameterType, const std::string&) = 0;

    /**
     * @brief Reattaches a previously detached dependent parameter.
     *
     * @param type       Parameter namespace/domain containing the parameter.
     * @param block_name Name of the parent block.
     * @param id         Identifier of the parameter to reattach.
     */
    virtual void reattach_parameter(ParameterType, const std::string&, const LhaID&) = 0;

    /**
     * @brief Detaches a dependent parameter from its upstream dependency graph.
     *
     * @param type       Parameter namespace/domain containing the parameter.
     * @param block_name Name of the parent block.
     * @param id         Identifier of the parameter to detach.
     */
    virtual void detach_parameter(ParameterType, const std::string&, const LhaID&) = 0;
};

#endif // ISTATDEPENDENCYPRUNER_H
