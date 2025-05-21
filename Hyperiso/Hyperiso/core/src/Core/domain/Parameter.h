/**
 * @file Parameter.h
 * @brief Defines the Parameter class for managing individual parameter values and the DependentParameter class for computed parameters.
 *
 * This file declares:
 * - Parameter: represents a basic parameter with expected value, uncertainties, and operational mode.
 * - DependentParameter: represents a parameter dynamically computed from other parameters.
 */

#ifndef PARAMETER_H
#define PARAMETER_H

#include <string>
#include <iostream>
#include <functional>

#include "Logger.h"
#include "General.h"
#include "Math.h"

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
 * @brief Represents a parameter with an ID, expected value, deviations, mode, and optional observers.
 *
 * Supports operations such as shifting, freezing, and notifying dependent parameters.
 */
class Parameter {
private:
    ParamId id;                                         ///< Unique identifier for the parameter.
    scalar_t expected;                                  ///< Expected value of the parameter.    
    scalar_t deviation_stat;                            ///< Statistical standard deviation.
    scalar_t deviation_syst;                            ///< Systematic standard deviation.
    scalar_t shift;                                     ///< Current shift applied to the parameter (0 if fixed).
    ParameterMode mode;                                 ///< Mode of operation (fixed or shiftable).
    std::vector<std::shared_ptr<Parameter>> observers;  ///< Observers notified when this parameter changes.

public:
    /**
     * @brief Default constructor initializes a null parameter.
     */
    inline Parameter() : id({ParameterType::SM, "NullBlock", 0}), expected(0), deviation_stat(0), deviation_syst(0), shift(0), mode(ParameterMode::FIXED) {}
    
    /**
     * @brief Constructs a Parameter with specified ID, mean value, and standard deviations.
     * @param id The ParamId object.
     * @param mean Expected value.
     * @param std_stat Statistical standard deviation.
     * @param std_syst Systematic standard deviation.
     */
    Parameter(ParamId id, scalar_t mean, scalar_t std_stat, scalar_t std_syst);

    /**
     * @brief Sets the operation mode of the parameter.
     * @param mode The new ParameterMode to set.
     */
    void set_mode(ParameterMode mode);

    /**
     * @brief Sets the standard deviations of the parameter.
     * @param stat New statistical standard deviation.
     * @param syst New systematic standard deviation.
     */
    void set_std(scalar_t stat, scalar_t syst);

    /**
     * @brief Sets the identifier of the parameter.
     * @param stat New identifier.
     */
    void set_id(ParamId id);

    /**
     * @brief Retrieves the current value of the parameter (shifted if allowed).
     * @return The current parameter value.
     */
    scalar_t get_val() const;

    /**
     * @brief Sets the expected value of the parameter and notifies observers.
     * @param val New expected value.
     */
    void set_expected(scalar_t val);

    /**
     * @brief Retrieves the combined standard deviation (quadratic sum of stat and syst).
     * @return The total standard deviation.
     */
    scalar_t get_combined_std() const;

    /**
     * @brief Retrieves the separated standard deviations (stat and syst).
     * @return A pair containing the statistical and systematic standard deviations.
     */
    std::pair<scalar_t, scalar_t> get_std() const;


    /**
     * @brief Retrieves the parameter's ID.
     * @return The ParamId object.
     */
    ParamId get_id() const;

    /**
     * @brief Changes the parameter's owner (model type).
     * @param type New ParameterType to set.
     */
    void set_owner(ParameterType type);

    /**
     * @brief Shifts the value of the parameter (only if SHIFTABLE).
     * @param shift Amount of shift to apply.
     * @throws std::runtime_error If called on a FIXED parameter.
     */
    void set_shift(scalar_t shift);

    /**
     * @brief Adds an observer parameter.
     * @param observer Shared pointer to the observer.
     */
    void addObserver(std::shared_ptr<Parameter> observer) { observers.push_back(observer); }

    /**
     * @brief Removes an observer parameter.
     * @param observer Shared pointer to the observer.
     */
    void removeObserver(std::shared_ptr<Parameter> observer) { observers.erase(std::find(observers.begin(), observers.end(), observer)); }

    /**
     * @brief Notifies all observers of a change.
     */
    void notifyObservers() {
        for (auto& observer : observers) {
            observer->update();
        }
    }

    /**
     * @brief Virtual method to update the parameter (default: no operation).
     */
    virtual void update() {}

    /**
     * @brief Virtual method to freeze the parameter (default: no operation).
     */
    virtual void freeze() {}

    /**
     * @brief Virtual method to unfreeze the parameter (default: no operation).
     */
    virtual void unfreeze() {}

    /**
     * @brief Assignment operator.
     * @param other The parameter to copy.
     * @return Reference to this parameter after copy.
     */
    Parameter& operator=(const Parameter& other);

    /**
     * @brief Overloaded stream insertion operator for printing the parameter details.
     * @param os Output stream.
     * @param p The parameter to print.
     * @return The output stream.
     */
    friend std::ostream& operator<<(std::ostream& os, const Parameter& p) {
        os << "Parameter " << p.id.block << "," << p.id.code << "=" << p.expected << "+-" << p.deviation_syst << "+-" << p.deviation_stat << std::endl;
        return os;
    }

    /**
     * @brief Overloads the += operator to add another Parameter to this one.
     *
     * Adds the numerical fields (expected, deviation_stat, deviation_syst, shift)
     * from the given Parameter to this instance.
     *
     * @param other The Parameter to add.
     * @return Reference to this Parameter after addition.
     */
    Parameter& operator+=(const Parameter& other);

};

class DependentParameter;

/**
 * @typedef DepParamUpdateFunc
 * @brief Function signature for updating a DependentParameter.
 *
 * Takes a map of source parameters and a shared pointer to the dependent parameter.
 */
typedef std::function<void(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>&, std::shared_ptr<DependentParameter>)> DepParamUpdateFunc;

/**
 * @class DependentParameter
 * @brief Represents a parameter whose value is dynamically computed from other parameters.
 *
 * Automatically updates itself when any of its source parameters changes.
 */
class DependentParameter : public Parameter, public std::enable_shared_from_this<DependentParameter> {
public:
    /**
     * @brief Constructs a DependentParameter from a set of source parameters and a recalculation function.
     * @param sources Map of source parameters.
     * @param recalculateFunc Lambda function to compute the value based on sources.
     */
    explicit DependentParameter(std::unordered_map<ParamId, std::shared_ptr<Parameter>> sources, DepParamUpdateFunc recalculateFunc) 
    : sources(std::move(sources)), recalculateLambda(std::move(recalculateFunc)), frozen(false) {}

    /**
     * @brief Checks if the dependent parameter depends on a given parameter.
     * @param pid ID of the parameter to check.
     * @return True if dependency exists, false otherwise.
     */
    bool dependsOn(const ParamId& pid);

    /**
     * @brief Initializes dependency tracking and observer registration.
     *
     * Must be called after construction.
     */
    void init();

    /**
     * @brief Updates the dependent parameter by recomputing its value.
     */
    void update() override;

    /**
     * @brief Freezes the dependent parameter (prevents update until unfreeze).
     */
    void freeze() override;

    /**
     * @brief Unfreezes the dependent parameter (updates if needed).
     */
    void unfreeze() override;

    /**
     * @brief Destructor. Cleans up dependency links.
     */
    ~DependentParameter();

private:
    std::shared_ptr<DependentParameter> self;                           ///< Self-reference used for observer management.
    std::unordered_map<ParamId, std::shared_ptr<Parameter>> sources;    ///< Source parameters.
    DepParamUpdateFunc recalculateLambda;                               ///< Recalculation function.
    bool frozen {false};                                                ///< If true, update is delayed.
    bool update_at_unfreeze {false};                                    ///< If true, an update is triggered upon unfreezing.
};

/**
 * @brief Overloads the += operator for shared pointers to Parameter.
 *
 * If both shared pointers are non-null, this function adds the contents
 * of rhs to lhs using the Parameter's += operator.
 *
 * @param lhs Shared pointer to the Parameter to be modified.
 * @param rhs Shared pointer to the Parameter to add.
 * @return Reference to lhs after the addition.
 */
inline std::shared_ptr<Parameter>& operator+=(std::shared_ptr<Parameter>& lhs, const std::shared_ptr<Parameter>& rhs) {
    if (rhs) {
        if (lhs)
            *lhs += *rhs;
        else
            lhs = std::make_shared<Parameter>(*rhs);
    }
    return lhs;
}

#endif // __PARAMETER_H__

