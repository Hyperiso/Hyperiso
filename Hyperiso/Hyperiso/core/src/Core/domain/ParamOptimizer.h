#ifndef PARAM_OPTIMIZER_H
#define PARAM_OPTIMIZER_H

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <variant>
#include <memory>
#include <stdexcept>
#include <string>
#include <algorithm>

#include "BlockAccessor.h"
#include "Block.h"
#include "Parameter.h"

/**
 * @file ParamOptimizer.h
 * @brief Defines the ParamOptimizer class to apply batched updates to parameters.
 *
 * This file declares:
 * - BAKeyHash: hash functor for (block, id) pairs.
 * - ParamOptimizer: a helper class to queue and commit parameter modifications
 *   across one or several BlockAccessor scopes, with coalescing and freeze/unfreeze
 *   of dependent structures.
 *
 * @ingroup ParametersModule
 * @see BlockAccessor
 * @see Block
 * @see Parameter
 */

 /**
 * @struct BAKeyHash
 * @brief Hash functor for (block name, parameter id string) pairs.
 *
 * Used internally by ParamOptimizer to keep track of operations
 * keyed by (BlockName, LhaID.to_string()).
 */
struct BAKeyHash {
    /**
     * @brief Computes the hash of a pair of strings.
     * @param k Pair (block name, id string).
     * @return Combined hash value.
     */
    std::size_t operator()(const std::pair<std::string,std::string>& k) const noexcept {
        std::size_t h1 = std::hash<std::string>()(k.first);
        std::size_t h2 = std::hash<std::string>()(k.second);
        return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1<<6) + (h1>>2));
    }
};

/**
 * @class ParamOptimizer
 * @brief Helper class to batch parameter updates on one or several BlockAccessor scopes.
 *
 * ParamOptimizer is designed to:
 * - Queue operations on parameters (set value, set full Parameter, remove).
 * - Optionally coalesce multiple operations on the same (block, id) into the last one.
 * - Apply all operations in a single commit, with:
 *   - freezing of all involved blocks before changes,
 *   - performing assignments / stores / removals,
 *   - a single notification per block when new parameters are stored,
 *   - unfreezing at the end.
 *
 * Typical usage:
 * @code
 * auto ba = Parameters::GetInstance(ParameterType::SM)->get_block_accessor();
 * ParamOptimizer opt(ba);
 * opt.set_value("MASS", LhaID(5), 4.8);
 * opt.set_value("MASS", LhaID(6), 173.0);
 * opt.commit();  // applies changes, notifies dependents
 * @endcode
 */
class ParamOptimizer {
public:
    /**
     * @brief Constructs a ParamOptimizer operating on a single BlockAccessor scope.
     * @param scope Shared pointer to the BlockAccessor to be modified.
     */
    explicit ParamOptimizer(std::shared_ptr<BlockAccessor> scope);

    /**
     * @brief Constructs a ParamOptimizer operating on multiple BlockAccessor scopes.
     *
     * This is useful when parameters are distributed over several
     * BlockAccessor instances (for example, several Parameters objects)
     * but blocks must remain unambiguous (a given block name must resolve
     * to the same Block in all scopes where it appears).
     *
     * @param scopes Vector of shared pointers to BlockAccessor scopes.
     */
    explicit ParamOptimizer(std::vector<std::shared_ptr<BlockAccessor>> scopes);

    /**
     * @brief Queue a simple numerical value assignment for a parameter.
     *
     * If the parameter already exists, it will be updated. Otherwise, it will be created
     * as a new Parameter with zero uncertainties.
     *
     * @param block Name of the block.
     * @param id LHA identifier of the parameter.
     * @param v New value to set.
     */
    void set_value(const BlockName& block, const LhaID& id, scalar_t v);

    /**
     * @brief Queue a full Parameter object assignment.
     *
     * If the parameter already exists, it will be overwritten (via Block::assign).
     * Otherwise, it will be stored (via Block::store).
     *
     * @param block Name of the block.
     * @param id LHA identifier of the parameter.
     * @param p Shared pointer to the Parameter object to use.
     */
    void set_param(const BlockName& block, const LhaID& id, std::shared_ptr<Parameter> p);

    /**
     * @brief Queue a removal operation for a parameter.
     *
     * If the parameter does not exist initially in the block, the operation is ignored.
     *
     * @param block Name of the block.
     * @param id LHA identifier of the parameter to remove.
     */
    void remove(const BlockName& block, const LhaID& id);

    /**
     * @brief Applies all queued operations.
     *
     * The algorithm is:
     * - Freeze all blocks in all scopes.
     * - Coalesce operations if requested (only the last op per (block, id) is kept).
     * - For each operation:
     *   - If it is a SetValue / SetParam and the parameter already existed, use Block::assign.
     *   - If it is a SetValue / SetParam on a new parameter, use Block::store and mark the
     *     block for notification.
     *   - If it is a Remove and the parameter existed initially, remove it (Block::remove).
     * - Notify observers once per block that received new parameters.
     * - Unfreeze all blocks.
     * - Clear the internal operation queue.
     *
     * @param coalesce If true (default), only the last operation per (block, id) is applied.
     */
    void commit(bool coalesce = true);

    /**
     * @brief Clears the pending operation queue without applying changes.
     */
    void clear();

private:
    /// @brief Internal operation: set a numeric value.
    struct OpSetValue { BlockName block; LhaID id; scalar_t value; };

    /// @brief Internal operation: set a full Parameter object.
    struct OpSetParam { BlockName block; LhaID id; std::shared_ptr<Parameter> param; };

    /// @brief Internal operation: remove a parameter.
    struct OpRemove   { BlockName block; LhaID id; };

    /// @brief Variant type collecting all possible operations.
    using Op = std::variant<OpSetValue, OpSetParam, OpRemove>;

    /**
     * @brief Coalesces the internal operation list by (block, id).
     *
     * Only the last operation per (block, id) pair is kept, preserving
     * the relative order of the surviving operations.
     *
     * @return A vector of coalesced operations.
     */
    std::vector<Op> coalesce_ops_() const;

    /**
     * @brief Finds the Block corresponding to the given block name in the configured scopes.
     *
     * The lookup rules are:
     * - Search all scopes in order.
     * - If no scope contains the block, throws std::invalid_argument.
     * - If more than one scope contains the block, and they point to different Block
     *   instances, throws std::invalid_argument (ambiguous block).
     *
     * @param name Name of the block to find.
     * @return Shared pointer to the resolved Block.
     */
    std::shared_ptr<Block> find_block_(const BlockName& name) const;

    /**
     * @brief Freezes all blocks in all scopes.
     *
     * Calls Block::freeze() on every block reachable from all BlockAccessor scopes.
     */
    void freeze_all_();

    /**
     * @brief Unfreezes all blocks in all scopes.
     *
     * Calls Block::unfreeze() on every block reachable from all BlockAccessor scopes.
     */
    void unfreeze_all_();

    /**
     * @brief BlockAccessor scopes on which operations will be applied.
     *
     * Each scope must provide consistent access to blocks; if a given block name
     * is present in multiple scopes, they must resolve to the same Block instance,
     * otherwise an ambiguity error is raised.
     */
    std::vector<std::shared_ptr<BlockAccessor>> scopes_;

    /**
     * @brief Internal list of queued operations.
     */
    std::vector<Op> ops_;
};

#endif