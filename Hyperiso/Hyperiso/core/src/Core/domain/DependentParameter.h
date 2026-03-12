#ifndef DEPENDENT_PARAMETER_H
#define DEPENDENT_PARAMETER_H

#include <functional>
#include <memory>
#include <unordered_map>

#include "Parameter.h"

class ParamSrc;
class DependentParameter;

/**
 * @file DependentParameter.h
 * @brief Defines parameters whose values are lazily computed from other parameters.
 *
 * This file declares:
 * - @ref DepParamUpdateFunc: callback type used to recompute a dependent parameter
 * - @ref DependentParameter: a @ref Parameter whose content is derived from source parameters
 *
 * A DependentParameter participates in the parameter dependency graph:
 * - it observes a set of source parameters,
 * - it becomes dirty when one of them changes,
 * - and it recomputes its cached value only on demand.
 *
 * This class is the parameter-level counterpart of @ref DependentBlock.
 *
 * @see Parameter
 * @see ParamSrc
 * @see Block
 */

/**
 * @typedef DepParamUpdateFunc
 * @brief Function signature used to recompute a DependentParameter.
 *
 * The callback receives:
 * - a @ref ParamSrc view over all source parameters,
 * - a shared pointer to the dependent parameter being updated.
 *
 * The callback is expected to modify the dependent parameter in-place
 * (for example with `set_expected_silent()`, `set_std()`, `set_scale()`, etc.).
 *
 * @note In the current implementation, the callback is called from
 *       @ref DependentParameter::ensure_up_to_date().
 */
typedef std::function<void(const ParamSrc&, std::shared_ptr<DependentParameter>)> DepParamUpdateFunc;

/**
 * @class DependentParameter
 * @brief Parameter whose value is computed from other parameters and cached lazily.
 *
 * A DependentParameter stores:
 * - a map of source parameters,
 * - a @ref ParamSrc helper wrapping these sources,
 * - a recomputation lambda,
 * - dirty/frozen flags controlling lazy updates.
 *
 * Lifecycle:
 * - after construction, @ref init() must be called once the object is managed by
 *   `std::shared_ptr`,
 * - @ref init() registers the parameter as observer of all its sources,
 * - when a source changes, @ref update() marks the parameter dirty and propagates,
 * - actual recomputation happens later in @ref ensure_up_to_date(), usually triggered
 *   by @ref get_val().
 *
 * Freezing:
 * - while frozen, recomputation is disabled,
 * - updates are deferred through `update_at_unfreeze`,
 * - once unfrozen, a pending update is replayed.
 *
 * Detach / reattach:
 * - @ref detach() removes all active dependencies while keeping the last computed value,
 * - @ref reattach() restores the saved dependencies and recomputation rule.
 *
 * @warning Unlike @ref Parameter::get_val(), this override currently returns the cached
 *          `expected` value directly after `ensure_up_to_date()`. It does **not**
 *          apply `mode/shift` logic from the base class.
 */
class DependentParameter : public Parameter {
public:
    /**
     * @brief Constructs a dependent parameter from its id, sources, and recomputation rule.
     *
     * The base @ref Parameter is initialized with:
     * - expected = 0
     * - statistical uncertainty = 0
     * - systematic uncertainty = 0
     *
     * The constructor stores:
     * - the active source map,
     * - a @ref ParamSrc wrapper over that map,
     * - the recomputation lambda,
     * - and a saved copy of sources/lambda used by @ref detach() / @ref reattach().
     *
     * The parameter is initially marked dirty.
     *
     * @note This constructor does **not** register observers. Call @ref init()
     *       once the object is owned by a `std::shared_ptr`.
     *
     * @param pid              Identifier of the dependent parameter.
     * @param sources          Source parameters keyed by @ref ParamId.
     * @param recalculateFunc  Callback used to recompute this parameter.
     */
    explicit DependentParameter(ParamId pid, std::unordered_map<ParamId, std::shared_ptr<Parameter>> sources, DepParamUpdateFunc recalculateFunc);

    /**
     * @brief Checks whether this parameter depends on a given source parameter.
     *
     * @param pid Identifier of the candidate source parameter.
     * @return True if @p pid is present in the current source map.
     */
    bool dependsOn(const ParamId& pid);

    /**
     * @brief Initializes observer registration on all source parameters.
     *
     * This method:
     * - stores a weak self-reference,
     * - registers this parameter as observer of each source parameter.
     *
     * @warning Must only be called when this object is already managed by
     *          `std::shared_ptr`, because it uses `shared_from_this()`.
     */
    void init();

    /**
     * @brief Marks the parameter as dirty and propagates the update downstream.
     *
     * Current behavior:
     * - if frozen:
     *   - sets `update_at_unfreeze = true`,
     *   - still calls @ref notifyObservers() immediately,
     *   - does not recompute now;
     * - otherwise:
     *   - sets `dirty = true`,
     *   - calls @ref notifyObservers().
     *
     * Actual recomputation is deferred to @ref ensure_up_to_date().
     */
    void update() override;

    /**
     * @brief Freezes the dependent parameter.
     *
     * While frozen:
     * - @ref ensure_up_to_date() becomes a no-op,
     * - @ref update() only records that an update is pending.
     */
    void freeze() override;

    /**
     * @brief Unfreezes the dependent parameter.
     *
     * If an update was requested while frozen, this method triggers one call to
     * @ref update() after unfreezing.
     */
    void unfreeze() override;

    /**
     * @brief Detaches this parameter from all upstream source parameters.
     *
     * If the parameter has been initialized, this removes this dependent parameter
     * from each source observer list.
     *
     * This does not clear local downstream observers.
     */
    void clear_above() override;

    /**
     * @brief Clears this parameter and all downstream dependent parameters.
     *
     * Current behavior:
     * - moves the current observer list to a local snapshot,
     * - detaches from all source parameters via @ref clear_above(),
     * - if an owning block exists, removes this parameter entry locally from that block
     *   using `owner_block->erase_local(id.code)`,
     * - recursively calls `clear_below()` on downstream observers.
     *
     * @note The erase is performed with `id.code`, consistent with block storage being
     *       keyed by @ref LhaID inside a block.
     */
    void clear_below() override;

    /**
     * @brief Returns the source parameters this dependent parameter currently uses.
     *
     * @return Map of source @ref ParamId to source @ref Parameter shared pointers.
     */
    std::unordered_map<ParamId, std::shared_ptr<Parameter>> get_source_parameters() const override {
        return this->sources_raw;
    }

    /**
     * @brief Returns the current cached value, recomputing it on demand if needed.
     *
     * This method:
     * - calls @ref ensure_up_to_date(),
     * - then returns the cached `expected` value directly.
     *
     * @note In the current implementation, this override does **not** apply
     *       `ParameterMode::SHIFTABLE` / `shift`. It returns `expected` only.
     *
     * @return Current cached dependent value.
     */
    scalar_t get_val() const override;
    
    /**
     * @brief Marks the cached value as dirty.
     *
     * This only sets the internal `dirty` flag to true.
     * No notification is emitted here.
     */
    void mark_dirty();

    /**
     * @brief Recomputes the parameter if its cached value is stale.
     *
     * Behavior:
     * - if not dirty: returns immediately,
     * - if frozen: returns immediately,
     * - otherwise:
     *   - forces all source parameters to be up to date by calling `src->get_val()`,
     *   - checks that the recomputation lambda exists,
     *   - clears the dirty flag,
     *   - calls the recomputation lambda on this parameter.
     *
     * @note The lambda is expected to update this parameter in-place.
     * @note Errors are logged if a source is null or if the lambda is missing.
     */
    void ensure_up_to_date();

    /**
     * @brief Replaces the dependency set and recomputation rule at runtime.
     *
     * This method:
     * - detaches from current sources,
     * - replaces active and saved sources,
     * - rebuilds the @ref ParamSrc helper,
     * - replaces the recomputation lambda,
     * - re-registers this parameter as observer of the new sources if already initialized,
     * - marks the parameter dirty,
     * - resets detached state,
     * - notifies downstream observers.
     *
     * @param new_sources New source parameter map.
     * @param new_lambda  New recomputation lambda.
     */
    void rebind(std::unordered_map<ParamId, std::shared_ptr<Parameter>> new_sources,
            DepParamUpdateFunc new_lambda);
        
    /**
     * @brief Temporarily detaches all dependencies while keeping the current cached value.
     *
     * Current behavior:
     * - returns immediately if already detached,
     * - forces one last recomputation via @ref ensure_up_to_date(),
     * - unregisters from all current sources,
     * - clears the active source map,
     * - rebuilds an empty @ref ParamSrc wrapper,
     * - resets the active recomputation lambda,
     * - clears the dirty flag,
     * - marks the dependency state as detached,
     * - notifies downstream observers.
     *
     * This is useful when one wants to “freeze” the dependency graph structure
     * while preserving the last computed numeric value.
     *
     * @see reattach()
     */
    void detach();

    /**
     * @brief Restores previously detached dependencies and recomputation rule.
     *
     * Current behavior:
     * - returns immediately if not detached,
     * - restores the saved source map,
     * - rebuilds the @ref ParamSrc wrapper,
     * - restores the saved recomputation lambda,
     * - re-registers this parameter as observer of restored sources,
     * - marks the parameter dirty,
     * - clears detached state,
     * - notifies downstream observers.
     *
     * After reattachment, the value will be recomputed lazily on the next access.
     *
     * @see detach()
     */
    void reattach();

    /**
     * @brief Destructor.
     *
     * If the parameter is still initialized and can lock its self-reference,
     * it unregisters itself from all currently active source parameters.
     *
     * This prevents stale observer pointers from remaining in source parameters.
     */
    ~DependentParameter();
    
private:
    std::weak_ptr<DependentParameter> self;                                 ///< Weak self-reference used for safe unregister/rebind logic.
    std::unordered_map<ParamId, std::shared_ptr<Parameter>> sources_raw;    ///< Active source parameters.
    std::unique_ptr<ParamSrc> sources;                                      ///< Helper view over active sources.
    DepParamUpdateFunc recalculateLambda;                                   ///< Active recomputation callback.
    bool frozen {false};                                                    ///< True if recomputation is currently frozen.
    bool update_at_unfreeze {false};                                        ///< True if an update was requested while frozen.
    bool dirty = true;                                                      ///< True if cached value must be recomputed before use.

    /**
     * @brief Saved copy of the source map used by @ref detach() / @ref reattach().
     *
     * This stores the “restorable” upstream dependency configuration.
     */
    std::unordered_map<ParamId, std::shared_ptr<Parameter>> saved_sources_raw;

    /**
     * @brief Saved recomputation callback used by @ref detach() / @ref reattach().
     */
    DepParamUpdateFunc saved_recalculateLambda;

    /**
     * @brief True when the active dependencies have been detached.
     */
    bool dependency_detached = false;
};

/**
 * @brief Adds the payload of one shared parameter into another.
 *
 * Behavior:
 * - if @p rhs is null: do nothing,
 * - if both are non-null and have the same @ref ParamId:
 *   performs `*lhs += *rhs`,
 * - otherwise replaces @p lhs with a new copied @ref Parameter built from `*rhs`.
 *
 * This helper is intended for convenient accumulation of parameter values stored
 * as `std::shared_ptr<Parameter>`.
 *
 * @param lhs Left-hand side parameter pointer.
 * @param rhs Right-hand side parameter pointer.
 * @return Reference to @p lhs after modification/replacement.
 */
inline std::shared_ptr<Parameter>& operator+=(std::shared_ptr<Parameter>& lhs, const std::shared_ptr<Parameter>& rhs) {
    if (rhs) {
        if (lhs && lhs->get_id() == rhs->get_id()) {
            *lhs += *rhs;
        }
        else
            lhs = std::make_shared<Parameter>(*rhs);
    }
    return lhs;
}

/**
 * @brief Scales a shared parameter in place.
 *
 * Current behavior:
 * - if @p rhs is non-zero, applies `*lhs *= rhs`,
 * - if @p rhs is zero, leaves @p lhs unchanged,
 * - if @p lhs is null, behavior is undefined because it is dereferenced.
 *
 * @warning Despite the operator spelling, this function mutates @p lhs in place
 *          and returns it by reference.
 *
 * @param lhs Shared parameter pointer to scale.
 * @param rhs Multiplicative factor.
 * @return Reference to @p lhs after scaling.
 */
inline std::shared_ptr<Parameter>& operator*(std::shared_ptr<Parameter>& lhs, const scalar_t& rhs) {
    if (rhs) {
        *lhs *= rhs;
    }
    return lhs;
}

#endif