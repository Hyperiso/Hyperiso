#ifndef IDEPENDENCYMOUNTER_H
#define IDEPENDENCYMOUNTER_H

/**
 * @file IDependencyPruner.h
 * @brief Interface for detaching and reattaching parameter/block dependencies.
 *
 * This interface defines a small port used to manipulate the dependency graph
 * of parameter repositories.
 *
 * It abstracts operations that temporarily prune dependency links:
 * - detach a whole block from its upstream dependencies,
 * - reattach a whole block,
 * - detach a single dependent parameter,
 * - reattach a single dependent parameter.
 *
 * The interface is generic over:
 * - the parameter repository / namespace identifier type,
 * - the block identifier type,
 * - the parameter identifier type.
 *
 * Typical concrete implementations delegate these operations to a backend
 * parameter repository such as @ref Parameters.
 *
 * @tparam Type_T  Type identifying the parameter namespace/domain
 *                 (e.g. ParameterType).
 * @tparam Block_T Type identifying a block
 *                 (e.g. BlockName).
 * @tparam Id_T    Type identifying a parameter inside a block
 *                 (e.g. LhaID).
 *
 * @see Parameters
 * @see BlockAccessor
 * @see DependentBlock
 * @see DependentParameter
 */

/**
 * @class IDependencyPruner
 * @brief Port interface for pruning and restoring dependency links.
 *
 * This interface is meant to be used by higher-level services that need to
 * temporarily disable automatic dependency propagation, for example:
 * - batch updates,
 * - silent mutations,
 * - optimization / fitting workflows,
 * - staged recomputation.
 *
 * The exact semantics of detach / reattach are backend-dependent, but in the
 * current HyperISO design they usually map to:
 * - block-level detach/reattach of @ref DependentBlock instances,
 * - parameter-level detach/reattach of @ref DependentParameter instances.
 *
 * @tparam Type_T  Parameter namespace/domain identifier.
 * @tparam Block_T Block identifier.
 * @tparam Id_T    Parameter identifier inside a block.
 */
template <typename Type_T, typename Block_T, typename Id_T> 
class IDependencyPruner {
public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~IDependencyPruner() = default;

    /**
     * @brief Reattaches a previously detached block.
     *
     * @param type       Parameter namespace/domain containing the block.
     * @param block_name Block identifier.
     */
    virtual void reattach_block(Type_T, const Block_T&) = 0;

    /**
     * @brief Detaches a block from its upstream dependency graph.
     *
     * After detachment, the block typically keeps its current cached content
     * but no longer reacts to source updates until reattached.
     *
     * @param type       Parameter namespace/domain containing the block.
     * @param block_name Block identifier.
     */
    virtual void detach_block(Type_T, const Block_T&) = 0;

    /**
     * @brief Reattaches a previously detached parameter.
     *
     * @param type       Parameter namespace/domain containing the parameter.
     * @param block_name Block identifier.
     * @param id         Parameter identifier inside the block.
     */
    virtual void reattach_parameter(Type_T, const Block_T&, const Id_T&) = 0;

    /**
     * @brief Detaches a parameter from its upstream dependency graph.
     *
     * After detachment, the parameter typically keeps its current cached value
     * but no longer reacts to source updates until reattached.
     *
     * @param type       Parameter namespace/domain containing the parameter.
     * @param block_name Block identifier.
     * @param id         Parameter identifier inside the block.
     */
    virtual void detach_parameter(Type_T, const Block_T&, const Id_T&) = 0;
};

#endif // IDEPENDENCYMOUNTER_H
