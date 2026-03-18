#ifndef DEPENDENCYPRUNER_H
#define DEPENDENCYPRUNER_H

#include "IDependencyPruner.h"
#include "BlockName.h"
#include "Include.h"
#include "Parameters.h"

/**
 * @file DependencyPruner.h
 * @brief Adapter implementing @ref IDependencyPruner on top of @ref Parameters.
 *
 * This class is a thin adapter between the dependency-pruning port
 * (@ref IDependencyPruner) and the concrete HyperISO parameter backend
 * (@ref Parameters).
 *
 * For each operation, it:
 * - retrieves the appropriate @ref Parameters instance via
 *   `Parameters::GetInstance(tp)`,
 * - delegates the detach / reattach request to the corresponding repository
 *   method.
 *
 * Supported operations:
 * - block-level detach / reattach,
 * - parameter-level detach / reattach.
 *
 * This adapter is useful when higher layers should depend only on a port
 * interface, without directly coupling to @ref Parameters.
 *
 * @see IDependencyPruner
 * @see Parameters
 * @see DependentBlock
 * @see DependentParameter
 */

/**
 * @class DependencyPruner
 * @brief Concrete dependency-pruning adapter backed by @ref Parameters.
 *
 * Template specialization:
 * - Type   = @ref ParameterType
 * - Block  = @ref BlockName
 * - Id     = @ref LhaID
 *
 * In practice, this class provides a simple façade to:
 * - `Parameters::reattach_block`
 * - `Parameters::detach_block`
 * - `Parameters::reattach_param`
 * - `Parameters::detach_param`
 */
class DependencyPruner : public IDependencyPruner<ParameterType, BlockName, LhaID> {
public:
    /**
     * @brief Reattaches a previously detached dependent block.
     *
     * Delegates to:
     * @code
     * Parameters::GetInstance(tp)->reattach_block(block_name);
     * @endcode
     *
     * @param tp         Parameter namespace containing the block.
     * @param block_name Name of the block to reattach.
     */
    void reattach_block(ParameterType tp, const BlockName& block_name) override;

    /**
     * @brief Detaches a dependent block from its upstream dependencies.
     *
     * Delegates to:
     * @code
     * Parameters::GetInstance(tp)->detach_block(block_name);
     * @endcode
     *
     * @param tp         Parameter namespace containing the block.
     * @param block_name Name of the block to detach.
     */
    void detach_block(ParameterType tp, const BlockName& block_name) override;

    /**
     * @brief Reattaches a previously detached dependent parameter.
     *
     * Delegates to:
     * @code
     * Parameters::GetInstance(tp)->reattach_param(block_name, id);
     * @endcode
     *
     * @param tp         Parameter namespace containing the parameter.
     * @param block_name Name of the parent block.
     * @param id         Identifier of the parameter to reattach.
     */
    void reattach_parameter(ParameterType tp, const BlockName& block_name, const LhaID& id) override;

    /**
     * @brief Detaches a dependent parameter from its upstream dependencies.
     *
     * Delegates to:
     * @code
     * Parameters::GetInstance(tp)->detach_param(block_name, id);
     * @endcode
     *
     * @param tp         Parameter namespace containing the parameter.
     * @param block_name Name of the parent block.
     * @param id         Identifier of the parameter to detach.
     */
    void detach_parameter(ParameterType tp, const BlockName& block_name, const LhaID& id) override;
};

#endif // DEPENDENCYPRUNER_H
