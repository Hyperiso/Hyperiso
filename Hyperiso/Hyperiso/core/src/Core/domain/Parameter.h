#ifndef PARAMETER_H
#define PARAMETER_H

#include <string>
#include <iostream>
#include <functional>

#include "Logger.h"
#include "Include.h"
#include "Math.h"

/**
 * @file Parameter.h
 * @brief Defines the Parameter class used to store individual physical/model parameters.
 *
 * This file declares @ref Parameter, the core scalar container used throughout the framework.
 *
 * A Parameter stores:
 * - a unique identifier (@ref ParamId),
 * - a central/expected value,
 * - statistical and systematic uncertainties,
 * - an optional additive shift,
 * - an operation mode (@ref ParameterMode),
 * - optional metadata such as scale and binning,
 * - and dependency/observer links to other parameters and to an owning @ref Block.
 *
 * The class participates in the dependency graph in two directions:
 * - upstream, through its owner block (if any),
 * - downstream, through its observer parameters.
 *
 * @see Block
 * @see ParamId
 */


class Block;

/**
 * @enum ParameterMode
 * @brief Defines the modes in which a parameter can operate.
 *
 * - FIXED: the parameter remains constant.
 * - SHIFTABLE: the parameter value can be shifted dynamically.
 */
enum class ParameterMode {
    FIXED,      ///< The parameter remains fixed at its expected value.
    SHIFTABLE   ///< The parameter value is shifted by an additive offset.
};

/**
 * @class Parameter
 * @brief Represents a single parameter with value, uncertainties, and dependency links.
 *
 * This is the fundamental data container for model, flavor, observable, and Wilson parameters.
 *
 * A Parameter stores:
 * - its identifier (@ref ParamId),
 * - its expected value,
 * - statistical and systematic standard deviations,
 * - an optional additive shift,
 * - optional metadata (scale, bin),
 * - a list of observer parameters,
 * - and optionally the block that owns it.
 *
 * Notification model:
 * - when the parameter changes, it can notify downstream parameter observers,
 * - if attached to a block, it can also trigger block-level observer propagation
 *   through the owning block.
 *
 * Read model:
 * - `get_val()` first ensures that the owning block (if any) is up to date,
 *   so accessing a parameter transparently triggers lazy recomputation of
 *   dependent blocks when needed.
 *
 * @see Block
 * @see ParameterMode
 */
class Parameter : public std::enable_shared_from_this<Parameter> {
protected:
    ParamId id;                                         ///< Unique identifier for the parameter.
    scalar_t expected;                                  ///< Expected value of the parameter.    
    scalar_t deviation_stat;                            ///< Statistical standard deviation.
    scalar_t deviation_syst;                            ///< Systematic standard deviation.
    scalar_t shift;                                     ///< Additive shift applied when the parameter is SHIFTABLE.
    ParameterMode mode;                                 ///< Current operating mode.
    std::vector<std::shared_ptr<Parameter>> observers;  ///< Observers notified when this parameter changes.
    std::weak_ptr<Block> owner_block;                   ///< Owning Block (if any) that this parameter belongs to.             
    std::optional<double> scale;                        ///< Optional renormalization scale associated with the parameter.
    std::optional<std::pair<double, double>> binning;   ///< Optional energy bin [low, high] associated with the parameter.

public:
    /**
     * @brief Default constructor.
     *
     * Creates a neutral parameter:
     * - id = (SM, "NullBlock", 0)
     * - expected = 0
     * - zero uncertainties
     * - zero shift
     * - mode = FIXED
     */
    inline Parameter() : id({ParameterType::SM, "NullBlock", 0}), expected(0), deviation_stat(0), deviation_syst(0), shift(0), mode(ParameterMode::FIXED) {}
    
    /**
     * @brief Constructs a parameter from its identifier, value, and uncertainties.
     *
     * The shift is initialized to zero and the mode to @ref ParameterMode::FIXED.
     *
     * @param id        Unique identifier of the parameter.
     * @param mean      Central value.
     * @param std_stat  Statistical uncertainty.
     * @param std_syst  Systematic uncertainty.
     */
    Parameter(ParamId id, scalar_t mean, scalar_t std_stat, scalar_t std_syst);

    /**
     * @brief Sets the parameter mode.
     *
     * This does not trigger notifications by itself.
     *
     * @param mode New mode to assign.
     */
    void set_mode(ParameterMode mode);

    /**
     * @brief Sets the statistical and systematic uncertainties.
     *
     * @param stat Statistical standard deviation.
     * @param syst Systematic standard deviation.
     */
    void set_std(scalar_t stat, scalar_t syst);

    /**
     * @brief Sets the parameter scale metadata.
     *
     * This is purely metadata storage; no notification is emitted.
     *
     * @param scale Scale value.
     */
    void set_scale(double scale);

    /**
     * @brief Sets the binning metadata of the parameter.
     *
     * @param bin Bin interval as `(low, high)`.
     */
    void set_bin(std::pair<double, double> bin);

    /**
     * @brief Replaces the identifier of the parameter.
     *
     * @param id New identifier.
     */
    void set_id(ParamId id);

    /**
     * @brief Returns the current value of the parameter.
     *
     * Behavior:
     * - if the parameter belongs to a block, the owning block is first asked to
     *   `ensure_up_to_date()`,
     * - then the returned value is:
     *   - `expected` if mode is @ref ParameterMode::FIXED,
     *   - `expected + shift` if mode is @ref ParameterMode::SHIFTABLE.
     *
     * @return Current parameter value.
     */

    virtual scalar_t get_val() const;

    /**
     * @brief Sets the expected value and notifies dependents.
     *
     * This updates the central value, then calls @ref notifyObservers().
     *
     * @param val New expected value.
     */
    void set_expected(scalar_t val);

    /**
     * @brief Sets the expected value without notifying observers.
     *
     * This is intended for internal/batched updates where propagation is handled separately.
     *
     * @warning This can leave dependents stale until notifications are manually triggered.
     *
     * @param val New expected value.
     */
    void set_expected_silent(scalar_t val);

    /**
     * @brief Returns the combined standard deviation.
     *
     * The combined uncertainty is computed as:
     * \f[
     * \sqrt{\sigma_{\mathrm{stat}}^2 + \sigma_{\mathrm{syst}}^2}
     * \f]
     *
     * @return Combined uncertainty.
     */
    scalar_t get_combined_std() const;

    /**
     * @brief Returns statistical and systematic uncertainties separately.
     *
     * @return Pair `(stat, syst)`.
     */

    std::pair<scalar_t, scalar_t> get_std() const;

    /**
     * @brief Returns the parameter scale metadata.
     *
     * If no scale has been assigned, returns `-1.0`.
     *
     * @return Scale value, or `-1.0` if undefined.
     */
    double get_scale();

    /**
     * @brief Returns the parameter bin metadata.
     *
     * If no bin has been assigned, returns `(-1.0, -1.0)`.
     *
     * @return Bin interval, or sentinel pair if undefined.
     */
    std::pair<double, double> get_bin();


    /**
     * @brief Returns the identifier of the parameter.
     *
     * @return Parameter id.
     */
    ParamId get_id() const;

    /**
     * @brief Changes the owner ParameterType part of the identifier.
     *
     * This is typically used when re-tagging a parameter as SM, BSM, WILSON, etc.
     *
     * @param type New owner type.
     */
    void set_owner(ParameterType type);

    /**
     * @brief Sets the additive shift of the parameter.
     *
     * This is only allowed if the parameter is in @ref ParameterMode::SHIFTABLE mode.
     *
     * @param shift New shift value.
     *
     * @throws std::runtime_error if called while the parameter is FIXED.
     */
    void set_shift(scalar_t shift);

    /**
     * @brief Adds an observer parameter.
     *
     * Observers are stored uniquely by pointer identity.
     *
     * @param observer Downstream observer parameter.
     */
    void addObserver(std::shared_ptr<Parameter> observer);

    /**
     * @brief Removes an observer parameter.
     *
     * Removal is performed by pointer identity.
     *
     * @param observer Observer to remove.
     */
    void removeObserver(std::shared_ptr<Parameter> observer);

    /**
     * @brief Notifies downstream parameter observers and then the owner block.
     *
     * Current behavior:
     * - first calls @ref notifyParamObserversOnly(),
     * - then, if an owning block exists, calls `owner_block->notifyObservers()`.
     *
     * This is the standard propagation entry point after a value-changing operation.
     */
    void notifyObservers();

    /**
     * @brief Notifies only parameter observers.
     *
     * This method:
     * - calls `update()` on each non-null observer parameter,
     * - removes null observer entries,
     * - does not notify the owning block.
     *
     * This is useful when propagating parameter-level dirtiness without re-triggering
     * block-level cascades multiple times.
     */
    void notifyParamObserversOnly();

    /**
     * @brief Hook called when an upstream dependency changes.
     *
     * Base implementation does nothing.
     * Derived classes may override this to recompute internal state.
     */
    virtual void update() {}

    /**
     * @brief Hook to freeze the parameter.
     *
     * Base implementation does nothing.
     * Derived classes may override this if they support explicit freezing.
     */
    virtual void freeze() {}

    /**
     * @brief Hook to unfreeze the parameter.
     *
     * Base implementation does nothing.
     */
    virtual void unfreeze() {}

    /**
     * @brief Hook to clear upstream dependencies.
     *
     * Base implementation does nothing.
     * Derived classes may override this to detach from source parameters.
     */
    virtual void clear_above() {}

    /**
     * @brief Recursively clears downstream dependencies below this parameter.
     *
     * Current behavior:
     * - takes ownership of the current observer list,
     * - clears the local observer list,
     * - calls @ref clear_above() on this parameter,
     * - recursively calls `clear_below()` on each former observer.
     *
     * This is used to tear down or invalidate a parameter dependency subtree.
     */
    virtual void clear_below();

    /**
     * @brief Returns the source parameters this parameter depends on.
     *
     * A plain Parameter has no explicit upstream parameter sources, so the default
     * implementation returns an empty map.
     *
     * Derived classes such as computed/dependent parameters should override this.
     *
     * @return Map of source `ParamId -> Parameter`.
     */
    virtual std::unordered_map<ParamId, std::shared_ptr<Parameter>> get_source_parameters() const;

    /**
     * @brief Sets the owning block of this parameter.
     *
     * @param owner Weak pointer to the owning block.
     */
    void set_owner_block(std::weak_ptr<Block> owner);

    /**
     * @brief Returns the owning block of this parameter.
     *
     * @return Weak pointer to the owning block.
     */
    std::weak_ptr<Block> get_owner_block() const;
    
    /**
     * @brief Overwrites the numerical/configuration payload from another parameter.
     *
     * Copied fields:
     * - expected value
     * - statistical uncertainty
     * - systematic uncertainty
     * - shift
     * - mode
     * - scale
     * - binning
     *
     * Preserved fields:
     * - identifier
     * - observers
     * - owning block
     *
     * @param other Source parameter payload.
     */
    void overwrite_payload_from(const Parameter& other);
    
    /**
     * @brief Assignment operator.
     *
     * Copies:
     * - identifier
     * - expected value
     * - uncertainties
     * - mode
     * - shift
     * - scale
     * - binning
     *
     * Does not copy:
     * - observer list
     * - owning block
     *
     * @param other Source parameter.
     * @return Reference to this parameter.
     */
    Parameter& operator=(const Parameter& other);

    /**
     * @brief Stream output operator for a parameter.
     *
     * Prints a compact human-readable representation including block/code,
     * optional binning, central value, and uncertainties.
     *
     * @param os Output stream.
     * @param p  Parameter to print.
     * @return Output stream.
     */
    friend std::ostream& operator<<(std::ostream& os, const Parameter& p);

    /**
     * @brief Adds another parameter payload to this one.
     *
     * Update rules:
     * - expected values are summed,
     * - statistical/systematic uncertainties are combined in quadrature,
     * - shifts are summed.
     *
     * Identifier, mode, metadata, and observer relationships are unchanged.
     *
     * @param other Parameter to add.
     * @return Reference to this parameter.
     */
    Parameter& operator+=(const Parameter& other);

    /**
     * @brief Multiplies the parameter by a scalar factor.
     *
     * Update rules:
     * - expected value is scaled by `factor`,
     * - uncertainties are scaled by `abs(factor)`,
     * - shift is scaled by `factor`.
     *
     * Metadata and observer relationships are unchanged.
     *
     * @param scale Multiplicative factor.
     * @return Reference to this parameter.
     */
    Parameter& operator*=(const scalar_t& scale);

    /**
     * @brief Virtual destructor.
     */
    virtual ~Parameter() = default;
};



#endif // PARAMETER_H

