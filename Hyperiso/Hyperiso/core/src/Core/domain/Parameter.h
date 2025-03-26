/**
 * @file Parameter.h
 * @brief Defines the Parameter class for managing individual parameter values.
 * 
 * This file declares the Parameter class, which encapsulates a parameter's ID,
 * expected value, deviation, and operational mode.
 */

#ifndef __PARAMETER_H__
#define __PARAMETER_H__

#include <string>
#include <iostream>
#include "Logger.h"
#include "General.h"

/**
 * @enum ParameterMode
 * @brief Defines the modes in which a parameter can operate.
 */
enum class ParameterMode {
    FIXED,      ///< The parameter remains constant.
    SHIFTABLE   ///< The parameter value can be modified.
};

/**
 * @class Parameter
 * @brief Represents a parameter with an ID, expected value, deviation, and mode.
 */
class Parameter {
private:
    ParamId id;              ///< Unique identifier for the parameter.
    double expected;         ///< Expected value of the parameter.
    double deviation_stat;   ///< Statistical standard deviation of the parameter.
    double deviation_syst;   ///< Systematic standard deviation of the parameter.
    double value;            ///< Current value of the parameter.
    ParameterMode mode;      ///< Mode of operation.

    std::vector<std::shared_ptr<Parameter>> observers;

public:
    /**
     * @brief Default constructor initializes a null parameter.
     */
    inline Parameter() : id({ParameterType::SM, "NullBlock", 0}), expected(0), deviation_stat(0), deviation_syst(0), mode(ParameterMode::FIXED) {}
    
    /**
     * @brief Constructs a Parameter with specified ID, mean value, and standard deviation.
     */
    Parameter(ParamId id, double mean, double std_stat, double std_syst);

    /**
     * @brief Sets the mode of the parameter.
     */
    void set_mode(ParameterMode mode);

    /**
     * @brief Sets the standard deviation of the parameter.
     */
    void set_std(double stat, double syst);

    /**
     * @brief Retrieves the current value of the parameter.
     */
    double get_val() const;

    /**
     * @brief Sets the standard deviation of the parameter.
     */
    void set_expected(double val);

    /**
     * @brief Retrieves the standard deviation of the parameter.
     */
    double get_std() const;

    /**
     * @brief Retrieves the ID of the parameter.
     */
    ParamId get_id() const;

    /**
     * @brief Retrieves the ID of the parameter.
     */
    void set_owner(ParameterType type);

    /**
     * @brief Shifts the parameter value if it is shiftable.
     * @throws std::runtime_error If the parameter is fixed.
     */
    void shift(double shift);

    void addObserver(std::shared_ptr<Parameter> observer) { observers.push_back(observer); }
    void removeObserver(std::shared_ptr<Parameter> observer) { observers.erase(std::find(observers.begin(), observers.end(), observer)); }
    void notifyObservers() {
        for (auto& observer : observers) {
            observer->update();
        }
    }
    virtual void update() {}

    /**
     * @brief Assignment operator.
     */
    Parameter& operator=(const Parameter& other) {
        this->id = other.id;
        this->expected = other.expected;
        this->deviation_stat = other.deviation_stat;
        this->deviation_syst = other.deviation_syst;
        this->mode = other.mode;
        this->value = other.value;
        return *this;
    }

    /**
     * @brief Overloaded stream insertion operator for printing parameter details.
     */
    friend std::ostream& operator<<(std::ostream& os, const Parameter& p) {
        os << "Parameter " << p.id.block << "," << p.id.code << "=" << p.expected << "+-" << p.deviation_syst << "+-" << p.deviation_stat << std::endl;
        return os;
    }
};

class DependentParameter;
typedef std::function<void(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>&, std::shared_ptr<DependentParameter>)> DepParamUpdateFunc;

class DependentParameter : public Parameter, public std::enable_shared_from_this<DependentParameter> {
public:
    explicit DependentParameter(std::unordered_map<ParamId, std::shared_ptr<Parameter>> sources, DepParamUpdateFunc recalculateFunc) 
    : sources(std::move(sources)), recalculateLambda(std::move(recalculateFunc)) {}

    bool dependsOn(const ParamId& pid) {
        return sources.contains(pid);
    }

    void init() {
        self = shared_from_this();
        if (self) {
            for (auto src : sources){
                src.second->addObserver(self);   
            }
        } else {
            std::cerr << "Error: DependentBlock must be created with std::make_shared!" << std::endl;
        }
    }

    void update() override {
        if (recalculateLambda 
            && std::all_of(sources.begin(), sources.end(), 
                        [](std::pair<ParamId, std::shared_ptr<Parameter>> block) { return block.second; })) 
        {
            if (auto self = shared_from_this()) { 
                recalculateLambda(sources, self);
            } else {
                std::cerr << "Error: shared_from_this() failed in update()" << std::endl;
            }
        }
    }

    ~DependentParameter() {
        LOG_INFO("Destruct DependentParameter at", self.get());
        if (self) {
            for (auto src : sources){
                src.second->removeObserver(self);   
            }
        }
    }

private:
    std::shared_ptr<DependentParameter> self;
    std::unordered_map<ParamId, std::shared_ptr<Parameter>> sources;
    DepParamUpdateFunc recalculateLambda;
};


#endif // __PARAMETER_H__

