#ifndef PARAMETERSHIFTER_H
#define PARAMETERSHIFTER_H

#include "IDataMutator.h"
#include "Parameters.h"

/**
 * @file ParameterShifter.h
 * @brief Concrete mutator that shifts parameter values.
 *
 * This header declares the ParameterShifter class, a concrete implementation
 * of ::IDataMutator that:
 *  - applies additive shifts to parameter values (rather than overwriting),
 *  - can also change the parameter mode via Parameters::changeParameterMode.
 *
 * It is particularly useful in fit / scan contexts where parameters are
 * updated relative to some initial value.
 */

/**
 * @class ParameterShifter
 * @ingroup DataMutationModule
 * @brief Concrete mutator that shifts parameter values rather than directly setting them.
 *
 * ParameterShifter differs from ParameterSetter in the semantics of mutate():
 *  - instead of assigning an absolute value,
 *  - it calls Parameters::shiftParameter(), adding the given amount to the
 *    current central value.
 *
 * Modes can be adjusted via `change_mode()` in the same way as ParameterSetter.
 */
class ParameterShifter : public IDataMutator<ParamId, scalar_t, ParameterMode> {
public:
    /**
     * @brief Shifts the value of a parameter by a given amount.
     *
     * If the ParamId has no type set, an error is logged, as the target
     * Parameters instance cannot be resolved.
     *
     * Internally this calls `Parameters::shiftParameter(pid, value)`,
     * so the new value is:
     * \f[
     *     p_{\text{new}} = p_{\text{old}} + \text{value}
     * \f]
     *
     * @param pid   The parameter ID.
     * @param value The amount by which to shift the parameter value.
     */
    void mutate(const ParamId& pid, scalar_t value) override;

    /**
     * @brief Changes the mode of a parameter (e.g., from FIXED to SHIFTABLE).
     *
     * Delegates to `Parameters::changeParameterMode` on the appropriate
     * Parameters instance, as determined by @p pid.type.
     *
     * @param pid  The parameter ID.
     * @param mode The new parameter mode.
     */
    void change_mode(const ParamId& pid, ParameterMode mode) override;
};


#endif // PARAMETERSHIFTER_H
