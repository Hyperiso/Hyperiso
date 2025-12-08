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
#include "Include.h"
#include "Math.h"

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
 * @brief Represents a parameter with an ID, expected value, deviations, mode, and optional observers.
 *
 * Supports operations such as shifting, freezing, and notifying dependent parameters.
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
    std::weak_ptr<Block> owner_block;                   
    std::optional<double> scale;
    std::optional<std::pair<double, double>> binning;

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
     * @brief Sets the renormalization scale of the parameter.
     * @param stat New statistical standard deviation.
     */
    void set_scale(double scale);

    /**
     * @brief Sets the energy binning of the parameter.
     * @param stat New statistical standard deviation.
     */
    void set_bin(std::pair<double, double> bin);

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
     * @brief Retrieves the renormalization scale of the parameter.
     * @return The renormalization scale.
     */
    double get_scale();

    /**
     * @brief Retrieves the energy binning of the parameter.
     * @return The energy bin.
     */
    std::pair<double, double> get_bin();


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
    void removeObserver(std::shared_ptr<Parameter> observer);

    /**
     * @brief Notifies all observers of a change.
     */
    void notifyObservers();

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
     * @brief Virtual method to unfreeze the parameter (default: no operation).
     */
    virtual void clear_above() {}

    /**
     * @brief Virtual method to unfreeze the parameter (default: no operation).
     */
    virtual void clear_below();

    virtual std::unordered_map<ParamId, std::shared_ptr<Parameter>> get_source_parameters() const {
        return {};
    }

    void set_owner_block(std::weak_ptr<Block> owner) { owner_block = std::move(owner); }
    std::weak_ptr<Block> get_owner_block() const { return owner_block; }

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
        os << "Parameter " << p.id.block << "," << p.id.code;
        if (p.binning.has_value()) {
            os << " [" << p.binning.value().first << "," << p.binning.value().second << "] ";
        }
        
        os << "=" << p.expected << "+-" << p.deviation_syst << "+-" << p.deviation_stat << std::endl;
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

    /**
     * @brief Overloads the *= operator to scale a parameter by a scalar.
     *
     * Performs the scaling of the different fields (expected, deviation_stat, deviation_syst, shift)
     * from the given Parameter to this instance.
     *
     * @param factor The scalar to scale by.
     * @return Reference to this Parameter after scaling.
     */
    Parameter& operator*=(const scalar_t& scale);

};



#endif // __PARAMETER_H__

