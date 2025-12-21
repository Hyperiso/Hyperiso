#ifndef DEPENDENT_PARAMETER_H
#define DEPENDENT_PARAMETER_H

#include <functional>
#include <memory>
#include <unordered_map>

#include "Parameter.h"

class ParamSrc;
class DependentParameter;

/**
 * @typedef DepParamUpdateFunc
 * @brief Function signature for updating a DependentParameter.
 *
 * The callback receives:
 *   - a ParamSrc view giving access to all source parameters
 *   - a shared pointer to the DependentParameter to update
 *
 * The lambda is expected to *modify* the dependent parameter in-place
 * (e.g. via set_expected / set_std / etc.).
 */
typedef std::function<void(const ParamSrc&, std::shared_ptr<DependentParameter>)> DepParamUpdateFunc;

/**
 * @class DependentParameter
 * @brief Parameter whose value is computed from other parameters (lazy / cached).
 *
 * A DependentParameter stores a set of *source* parameters and a user-provided
 * recomputation lambda (@ref DepParamUpdateFunc). It implements an observer
 * mechanism:
 *  - In @ref init(), it registers itself as observer of each source parameter.
 *  - When a source changes, @ref update() is called:
 *      - it does NOT recompute immediately,
 *      - it marks the parameter as dirty and notifies its own observers.
 *  - The actual recomputation happens on demand in @ref ensure_up_to_date(),
 *    typically triggered by calling @ref get_val().
 *
 * Freezing:
 *  - When frozen, updates are deferred (a flag is set) and no recomputation is
 *    performed until @ref unfreeze() is called.
 *
 * Lifetime / safety notes:
 *  - @ref init() must be called once after construction, when managed by a
 *    std::shared_ptr (so shared_from_this() is valid).
 *  - The class stores a weak self-reference to safely unregister from sources
 *    in destructor / clear logic.
 */
class DependentParameter : public Parameter {
public:
    /**
     * @brief Constructs a dependent parameter with given sources and recomputation rule.
     *
     * Internally initializes the base @ref Parameter with:
     *  - expected = 0, stat = 0, syst = 0
     * and stores:
     *  - the raw source map,
     *  - a @ref ParamSrc view wrapper (used for convenient access / error context),
     *  - the recomputation lambda.
     *
     * @note Dependency registration is NOT performed here. Call @ref init()
     *       once the object is owned by a std::shared_ptr.
     *
     * @param pid             ID of the dependent parameter.
     * @param sources         Source parameters keyed by ParamId.
     * @param recalculateFunc Lambda that recomputes this parameter in-place.
     */
    explicit DependentParameter(ParamId pid, std::unordered_map<ParamId, std::shared_ptr<Parameter>> sources, DepParamUpdateFunc recalculateFunc);

    /**
     * @brief Checks if this dependent parameter depends on the given parameter id.
     *
     * @param pid ID of the source parameter.
     * @return True if pid is present in the source map, false otherwise.
     */
    bool dependsOn(const ParamId& pid);

    /**
     * @brief Initializes dependency tracking and observer registration.
     *
     * This method:
     *  - stores a weak self-reference (for safe unregistration later),
     *  - registers this parameter as an observer of each non-null source parameter
     *    using @ref Parameter::addObserver().
     *
     * @warning Must be called exactly once after construction, and only when this
     *          object is already managed by std::shared_ptr (shared_from_this valid).
     */
    void init();

    /**
     * @brief Marks the parameter as dirty and propagates change downstream.
     *
     * This override implements a *lazy* update strategy:
     *  - If frozen: sets a flag so an update will be triggered at unfreeze time,
     *    then notifies observers immediately.
     *  - Otherwise: marks the parameter dirty and notifies observers.
     *
     * Actual recomputation is deferred to @ref ensure_up_to_date(), typically
     * called by @ref get_val().
     */
    void update() override;

    /**
     * @brief Freezes this parameter (disables recomputation).
     *
     * While frozen:
     *  - @ref ensure_up_to_date() will not recompute,
     *  - @ref update() will only set a pending-update flag and notify observers.
     */
    void freeze() override;

    /**
     * @brief Unfreezes this parameter and triggers a pending update if any.
     *
     * If updates occurred while frozen, a call to @ref update() is triggered once,
     * which marks the parameter dirty and propagates to observers.
     */
    void unfreeze() override;

    /**
     * @brief Unregisters this DependentParameter from all its source parameters.
     *
     * Removes this parameter from each source observer list.
     * Safe even if already detached.
     */
    void clear_above() override;

    /**
     * @brief Clears this parameter and all downstream dependents.
     *
     * Steps:
     *  - detaches from all sources via @ref clear_above(),
     *  - asks the owning Block (if any) to erase this local parameter entry
     *    via Block::erase_local(code),
     *  - recursively calls clear_below() on its observers.
     */
    void clear_below() override;

    /**
     * @brief Retrieves all the source blocks of the dependent block.
     * @return A map of source block names to their shared pointers.
     */
    std::unordered_map<ParamId, std::shared_ptr<Parameter>> get_source_parameters() const override {
        return this->sources_raw;
    }

    /**
     * @brief Returns the current value, recomputing it on demand if needed.
     *
     * This method forces recomputation through @ref ensure_up_to_date() if the
     * parameter is marked dirty and not frozen.
     *
     * @return Current expected value (cached after recomputation).
     */
    scalar_t get_val() const override;
    
    /**
     * @brief Marks the cached value as invalid.
     *
     * The next call to @ref ensure_up_to_date() (or @ref get_val()) will recompute
     * the value (unless frozen).
     */
    void mark_dirty();

    /**
     * @brief Ensures the cached value is up-to-date by recomputing if needed.
     *
     * If not dirty, does nothing.
     * If frozen, does nothing.
     *
     * Otherwise:
     *  - touches each source via src->get_val() to ensure sources are up-to-date,
     *  - calls the recomputation lambda to update this parameter in-place
     *    (expected/std/etc.),
     *  - clears the dirty flag.
     *
     * @note If a source pointer is null or the lambda is not set, an error is logged.
     */
    void ensure_up_to_date();

    /**
     * @brief Replaces the dependency set and recomputation rule at runtime.
     *
     * This method:
     *  - detaches from current sources (clear_above),
     *  - replaces the source map and the ParamSrc wrapper,
     *  - replaces the recomputation lambda,
     *  - re-registers as observer of the new sources (if already initialized),
     *  - marks the parameter dirty and notifies observers.
     *
     * @param new_sources New source parameter map.
     * @param new_lambda  New recomputation lambda.
     */
    void rebind(std::unordered_map<ParamId, std::shared_ptr<Parameter>> new_sources,
            DepParamUpdateFunc new_lambda);
            
    /**
     * @brief Destructor: unregisters from sources if still initialized.
     *
     * If @ref init() has been called and the self weak pointer can be locked,
     * removes this dependent parameter from each source observer list.
     *
     * This prevents dangling observer entries in source parameters.
     */
    ~DependentParameter();
    
private:
    std::weak_ptr<DependentParameter> self;                                 ///< Self-reference used for safe observer management in destructor / clear_*.
    std::unordered_map<ParamId, std::shared_ptr<Parameter>> sources_raw;    ///< Raw source parameters, keyed by ParamId.
    std::unique_ptr<ParamSrc> sources;                                      ///< View wrapper around sources_raw (with context string for nicer errors).
    DepParamUpdateFunc recalculateLambda;                                   ///< Lambda used to recompute the value of this parameter.
    bool frozen {false};                                                    ///< If true, update is delayed.
    bool update_at_unfreeze {false};                                        ///< If true, an update is triggered upon unfreezing.
    bool dirty = true;                                                      ///< True if cached value is invalid and must be recomputed on demand.
};

/**
 * @brief Overloads the += operator for shared_ptr<Parameter>.
 *
 * Behaviour:
 *   - If rhs is null: do nothing.
 *   - If lhs is non-null and has the same ParamId as rhs:
 *       *lhs += *rhs;
 *   - Else:
 *       lhs is replaced by a new Parameter copied from *rhs.
 *
 * @param lhs Shared pointer to the Parameter to be modified.
 * @param rhs Shared pointer to the Parameter to add.
 * @return Reference to lhs after the operation.
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
 * @brief Overloads the * operator for shared_ptr<Parameter> and a scalar factor.
 *
 * Behaviour:
 *   - If rhs (scale factor) is non-zero, scales *lhs in-place using
 *     Parameter::operator*=.
 *   - If lhs is null, behaviour is undefined (same as dereferencing a null
 *     shared_ptr elsewhere).
 *
 * @param lhs Shared pointer to the Parameter to be scaled.
 * @param rhs Scale factor.
 * @return Reference to lhs after scaling.
 */
inline std::shared_ptr<Parameter>& operator*(std::shared_ptr<Parameter>& lhs, const scalar_t& rhs) {
    if (rhs) {
        *lhs *= rhs;
    }
    return lhs;
}

#endif