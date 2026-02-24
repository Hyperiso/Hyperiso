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
 * @brief Defines the Parameter class for managing individual parameter values and the DependentParameter class for computed parameters.
 *
 * This file declares:
 * - Parameter: represents a basic parameter with expected value, uncertainties, and operational mode.
 * - DependentParameter: represents a parameter dynamically computed from other parameters.
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
    FIXED,      ///< The parameter remains constant.
    SHIFTABLE   ///< The parameter value can be modified.
};

/**
 * @class Parameter
 * @brief Represents a parameter with an ID, central value, uncertainties, mode, and optional observers.
 *
 * This class is the core representation of a physical or model parameter:
 *   - it stores an identifier (@ref ParamId),
 *   - a central (expected) value,
 *   - statistical and systematic uncertainties,
 *   - an optional shift (for nuisance / fit parameters),
 *   - and optional metadata such as renormalization scale and energy binning.
 *
 * It also implements a simple observer pattern: a parameter can depend on other
 * parameters and be notified when they change via @ref notifyObservers() and @ref update().
 */
class Parameter : public std::enable_shared_from_this<Parameter> {
protected:
    ParamId id;                                         ///< Unique identifier for the parameter.
    scalar_t expected;                                  ///< Expected value of the parameter.    
    scalar_t deviation_stat;                            ///< Statistical standard deviation.
    scalar_t deviation_syst;                            ///< Systematic standard deviation.
    scalar_t shift;                                     ///< Current shift applied to the parameter (0 if fixed).
    ParameterMode mode;                                 ///< Mode of operation (fixed or shiftable).
    std::vector<std::shared_ptr<Parameter>> observers;  ///< Observers notified when this parameter changes.
    std::weak_ptr<Block> owner_block;                   ///< Owning Block (if any) that this parameter belongs to.             
    std::optional<double> scale;                        ///< Optional renormalization scale associated with the parameter.
    std::optional<std::pair<double, double>> binning;   ///< Optional energy bin [low, high] associated with the parameter.

public:
    /**
     * @brief Default constructor initializes a "null" parameter.
     *
     * The default parameter:
     *  - has ID (SM, "NullBlock", 0),
     *  - central value 0,
     *  - zero statistical and systematic uncertainties,
     *  - zero shift,
     *  - and is FIXED.
     */
    inline Parameter() : id({ParameterType::SM, "NullBlock", 0}), expected(0), deviation_stat(0), deviation_syst(0), shift(0), mode(ParameterMode::FIXED) {}
    
    /**
     * @brief Constructs a Parameter with specified ID, mean value, and standard deviations.
     * @param id        The ParamId object for this parameter.
     * @param mean      Expected (central) value.
     * @param std_stat  Statistical standard deviation.
     * @param std_syst  Systematic standard deviation.
     */
    Parameter(ParamId id, scalar_t mean, scalar_t std_stat, scalar_t std_syst);

    /**
     * @brief Sets the operation mode of the parameter.
     * @param mode The new ParameterMode to set (FIXED or SHIFTABLE).
     */
    void set_mode(ParameterMode mode);

    /**
     * @brief Sets both statistical and systematic standard deviations.
     * @param stat New statistical standard deviation.
     * @param syst New systematic standard deviation.
     */
    void set_std(scalar_t stat, scalar_t syst);

    /**
     * @brief Sets the renormalization scale associated with this parameter.
     * @param scale Renormalization scale (e.g. in GeV).
     */
    void set_scale(double scale);

    /**
     * @brief Sets the energy bin associated with this parameter.
     * @param bin Pair (low, high) defining the bin edges.
     */
    void set_bin(std::pair<double, double> bin);

    /**
     * @brief Sets the identifier of the parameter.
     * @param id New ParamId to assign.
     */
    void set_id(ParamId id);

    /**
     * @brief Retrieves the current value of the parameter.
     *
     * If the mode is:
     *  - FIXED     → returns the expected value,
     *  - SHIFTABLE → returns expected + shift.
     *
     * @return The current parameter value.
     */
    virtual scalar_t get_val() const;

    /**
     * @brief Sets the expected (central) value of the parameter and notifies observers.
     * @param val New expected value.
     */
    void set_expected(scalar_t val);

        /**
     * @brief Sets the expected (central) value without notifying observers.
     *
     * This method updates the internal expected value but deliberately
     * skips observer notification. It is intended for internal use
     * (e.g. batch updates, reconstruction, deserialization).
     *
     * @warning Using this method may leave dependent parameters in an
     *          inconsistent state if observers are not updated manually.
     *
     * @param val New expected value.
     */
    void set_expected_silent(scalar_t val);

    /**
     * @brief Retrieves the combined standard deviation.
     *
     * The combined uncertainty is defined as the quadratic sum:
     * \f$\sqrt{\sigma_\text{stat}^2 + \sigma_\text{syst}^2}\f$.
     *
     * @return The total standard deviation.
     */
    scalar_t get_combined_std() const;

    /**
     * @brief Retrieves statistical and systematic standard deviations separately.
     * @return A pair (stat, syst).
     */
    std::pair<scalar_t, scalar_t> get_std() const;

    /**
     * @brief Retrieves the renormalization scale of the parameter.
     *
     * If no scale has been set, returns -1.0 as a sentinel.
     *
     * @return The renormalization scale (or -1.0 if undefined).
     */
    double get_scale();

    /**
     * @brief Retrieves the energy bin associated with the parameter.
     *
     * If no binning has been set, returns (-1.0, -1.0) as a sentinel.
     *
     * @return The energy bin [low, high].
     */
    std::pair<double, double> get_bin();


    /**
     * @brief Retrieves the parameter's ID.
     * @return The ParamId object.
     */
    ParamId get_id() const;

    /**
     * @brief Changes the parameter's owner type (ParameterType field of ParamId).
     *
     * This is typically used to re-tag a parameter as SM, BSM, FLAVOR, etc.
     *
     * @param type New ParameterType to set.
     */
    void set_owner(ParameterType type);

    /**
     * @brief Shifts the value of the parameter (only if SHIFTABLE).
     *
     * The shift is additive on top of the expected value:
     *   current value = expected + shift.
     *
     * @param shift Amount of shift to apply.
     * @throws std::runtime_error If called on a FIXED parameter.
     */
    void set_shift(scalar_t shift);

    /**
     * @brief Adds an observer parameter.
     *
     * The observer will be notified via @ref update() when this parameter
     * changes (for example after @ref set_expected()).
     *
     * @param observer Shared pointer to the observer.
     */
    void addObserver(std::shared_ptr<Parameter> observer);

    /**
     * @brief Removes an observer from the observer list.
     *
     * The observer is matched by pointer identity.
     * If the observer is not present, this function has no effect.
     *
     * @param observer Shared pointer to the observer to remove.
     */
    void removeObserver(std::shared_ptr<Parameter> observer);

    /**
     * @brief Notifies all observer parameters of a change.
     *
     * For each registered observer:
     *  - calls observer->update(),
     *  - ignores null observers,
     *  - removes expired/null entries from the observer list.
     *
     * After notifying parameter observers, the owning Block (if any)
     * is also notified via Block::notifyObservers().
     *
     * This method is automatically called after state-changing operations
     * such as @ref set_expected().
     */
    void notifyObservers();

    /**
     * @brief Notify only parameter observers (does NOT notify owner block).
     *
     * Use this to propagate dirtiness to DependentParameters without causing
     * DependentBlocks (which observe the owner Block) to be spammed N times.
     */
    void notifyParamObserversOnly();

    /**
     * @brief Virtual hook called when a source parameter has changed.
     *
     * Base implementation does nothing. Derived classes can override this
     * to recompute their value based on dependencies.
     */
    virtual void update() {}

    /**
     * @brief Virtual hook to freeze the parameter (e.g. disable updates).
     *
     * Base implementation does nothing. Derived classes can override this
     * to implement freezing logic.
     */
    virtual void freeze() {}

    /**
     * @brief Virtual hook to unfreeze the parameter.
     *
     * Base implementation does nothing. Derived classes can override this.
     */
    virtual void unfreeze() {}

    /**
     * @brief Virtual hook to clear dependent parameters above this node.
     *
     * Base implementation does nothing. Derived classes can override this
     * to clear or reset state "upstream" in a dependency graph.
     */
    virtual void clear_above() {}

    /**
     * @brief Recursively clears this parameter and all downstream dependents.
     *
     * The clearing procedure is:
     *  - move the current observer list to a local container,
     *  - call @ref clear_above() on this parameter,
     *  - recursively call clear_below() on each observer.
     *
     * This is typically used to invalidate cached or derived quantities
     * in a dependency graph.
     */
    virtual void clear_below();

    /**
     * @brief Returns the parameters this parameter depends on.
     *
     * For a plain Parameter, there are no dependencies and an empty
     * map is returned.
     *
     * Derived classes representing computed or composite parameters
     * should override this method to expose their dependency graph.
     *
     * @return Map of ParamId → source Parameter.
     */
    virtual std::unordered_map<ParamId, std::shared_ptr<Parameter>> get_source_parameters() const;

    /**
     * @brief Sets the owning block of this parameter.
     * @param owner Weak pointer to the owning Block.
     */
    void set_owner_block(std::weak_ptr<Block> owner);

    /**
     * @brief Returns the owning block of this parameter, if any.
     * @return Weak pointer to the owning Block.
     */
    std::weak_ptr<Block> get_owner_block() const;
    
    /**
     * @brief Overwrites the numerical payload from another Parameter.
     *
     * Copies only the physical content:
     *  - expected value,
     *  - statistical and systematic uncertainties,
     *  - shift,
     *  - mode,
     *  - scale,
     *  - binning.
     *
     * The following are NOT modified:
     *  - parameter ID,
     *  - observers,
     *  - owning Block.
     *
     * @param other Parameter to copy the payload from.
     */
    void overwrite_payload_from(const Parameter& other);
    
    /**
     * @brief Assignment operator.
     *
     * Copies all physical and configuration fields from @p other:
     *  - id,
     *  - expected value,
     *  - uncertainties,
     *  - mode,
     *  - shift,
     *  - scale,
     *  - binning.
     *
     * Observer relationships and owning Block are intentionally
     * NOT copied.
     *
     * @param other Parameter to copy from.
     * @return Reference to this instance.
     */
    Parameter& operator=(const Parameter& other);

    /**
     * @brief Prints a human-readable representation of the parameter.
     *
     * Output format:
     *   Parameter BLOCK,ID [bin_low,bin_high] = value +- syst +- stat
     *
     * The bin information is omitted if no binning is defined.
     *
     * @param os Output stream.
     * @param p  Parameter to print.
     * @return Output stream.
     */
    friend std::ostream& operator<<(std::ostream& os, const Parameter& p);

    /**
     * @brief Adds another parameter to this one.
     *
     * Numerical fields are combined as:
     *  - expected       ← expected + other.expected
     *  - deviation_stat ← sqrt(stat^2 + other.stat^2)
     *  - deviation_syst ← sqrt(syst^2 + other.syst^2)
     *  - shift          ← shift + other.shift
     *
     * Metadata (ID, mode, scale, binning, observers) is unchanged.
     *
     * @param other Parameter to add.
     * @return Reference to this parameter.
     */
    Parameter& operator+=(const Parameter& other);

    /**
     * @brief Scales the parameter by a scalar factor.
     *
     * Scaling rules:
     *  - expected       ← expected * factor
     *  - deviation_stat ← deviation_stat * |factor|
     *  - deviation_syst ← deviation_syst * |factor|
     *  - shift          ← shift * factor
     *
     * @param scale Scalar factor.
     * @return Reference to this parameter.
     */
    Parameter& operator*=(const scalar_t& scale);

    virtual ~Parameter() = default;
};



#endif // __PARAMETER_H__

