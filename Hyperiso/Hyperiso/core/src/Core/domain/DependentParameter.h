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
 * @brief Parameter whose value is dynamically computed from other parameters.
 *
 * This class implements a simple observer pattern:
 *   - It registers itself as an observer of its source parameters in init().
 *   - When any source parameter changes, it receives update(), and calls
 *     the provided DepParamUpdateFunc to recompute its own value.
 *   - It can be frozen/unfrozen to postpone updates.
 */
class DependentParameter : public Parameter {
public:
    /**
     * @brief Constructs a DependentParameter.
     *
     * @param pid   ID of the dependent parameter.
     * @param sources Map of source parameters (by ParamId).
     * @param recalculateFunc Lambda used to recompute the parameter's value.
     *                        It will be called from update().
     */
    explicit DependentParameter(ParamId pid, std::unordered_map<ParamId, std::shared_ptr<Parameter>> sources, DepParamUpdateFunc recalculateFunc);

    /**
     * @brief Checks if this dependent parameter depends on a given source.
     * @param pid ID of the parameter to check.
     * @return True if pid is among the sources, false otherwise.
     */
    bool dependsOn(const ParamId& pid);

    /**
     * @brief Initializes dependency tracking and observer registration.
     *
     * Must be called once after construction, when the DependentParameter is
     * already owned by a std::shared_ptr (so shared_from_this() is valid).
     */
    void init();

    /**
     * @brief Recomputes the dependent parameter based on its sources.
     *
     * If frozen, the update is delayed until unfreeze() is called.
     * If some source is null or the recalculate lambda is not set,
     * logs an error and does nothing.
     */
    void update() override;

    /**
     * @brief Freezes the dependent parameter (defers updates).
     *
     * Any update requested while frozen will set a flag and be executed
     * once unfreeze() is called.
     */
    void freeze() override;

    /**
     * @brief Unfreezes the parameter and triggers a pending update if any.
     */
    void unfreeze() override;

    /**
     * @brief Breaks links from this parameter to its sources.
     *
     * Removes this DependentParameter as an observer from all its sources.
     */
    void clear_above() override;

    /**
     * @brief Recursively clears this parameter and all its dependents.
     *
     * - Unregisters from sources (clear_above())
     * - Requests its owning Block (if any) to erase this parameter
     * - Recursively calls clear_below() on its own observers.
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
     * @brief Returns the raw map of source parameters.
     */
    ~DependentParameter();
    
private:
    std::weak_ptr<DependentParameter> self;                                 ///< Self-reference used for safe observer management in destructor / clear_*.
    std::unordered_map<ParamId, std::shared_ptr<Parameter>> sources_raw;    ///< Raw source parameters, keyed by ParamId.
    std::unique_ptr<ParamSrc> sources;                                      ///< View wrapper around sources_raw (with context string for nicer errors).
    DepParamUpdateFunc recalculateLambda;                                   ///< /// Lambda used to recompute the value of this parameter.
    bool frozen {false};                                                    ///< If true, update is delayed.
    bool update_at_unfreeze {false};                                        ///< If true, an update is triggered upon unfreezing.
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