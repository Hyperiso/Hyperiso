#ifndef DEPENDENT_PARAMETER_H
#define DEPENDENT_PARAMETER_H

#include "Parameter.h"

class ParamSrc;

class DependentParameter;

/**
 * @typedef DepParamUpdateFunc
 * @brief Function signature for updating a DependentParameter.
 *
 * Takes a map of source parameters and a shared pointer to the dependent parameter.
 */
// using DepParamUpdateFunc = std::function<scalar_t(const ParamSrc&)>;
typedef std::function<void(const ParamSrc&, std::shared_ptr<DependentParameter>)> DepParamUpdateFunc;

/**
 * @class DependentParameter
 * @brief Represents a parameter whose value is dynamically computed from other parameters.
 *
 * Automatically updates itself when any of its source parameters changes.
 */
class DependentParameter : public Parameter {
public:
    /**
     * @brief Constructs a DependentParameter from a set of source parameters and a recalculation function.
     * @param pid ID of the dependent parameter.
     * @param sources Map of source parameters.
     * @param recalculateFunc Lambda function to compute the value based on sources.
     */
    explicit DependentParameter(ParamId pid, std::unordered_map<ParamId, std::shared_ptr<Parameter>> sources, DepParamUpdateFunc recalculateFunc);

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

    void clear_above() override;

    void clear_below() override;

    /**
     * @brief Destructor. Cleans up dependency links.
     */
    ~DependentParameter();

private:
    std::weak_ptr<DependentParameter> self;                           ///< Self-reference used for observer management.
    std::unordered_map<ParamId, std::shared_ptr<Parameter>> sources_raw;    ///< Source parameters.
    std::unique_ptr<ParamSrc> sources;                                      ///< Source parameters.
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
        if (lhs && lhs->get_id() == rhs->get_id()) {
            *lhs += *rhs;
        }
        else
            lhs = std::make_shared<Parameter>(*rhs);
    }
    return lhs;
}

/**
 * @brief Overloads the * operator for a shared pointer to Parameter.
 *
 * If the shared pointer is non-null, this function scales the contents
 * of rhs by lhs using the Parameter's *= operator.
 *
 * @param lhs Shared pointer to the Parameter to be modified.
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