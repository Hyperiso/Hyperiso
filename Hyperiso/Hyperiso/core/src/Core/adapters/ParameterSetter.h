#ifndef PARAMETERMODIFIER_H
#define PARAMETERMODIFIER_H

#include "IDataMutator.h"
#include "Parameters.h"
#include "MemoryManager.h"
#include "Math.h"

/**
 * @file ParameterModifier.h
 * @brief Concrete mutator that sets parameter values and modes.
 *
 * This header declares the ParameterSetter class, a concrete implementation
 * of ::IDataMutator that:
 *  - directly sets the central value of a parameter in Parameters,
 *  - updates its mode via Parameters::changeParameterMode.
 *
 * It is designed for use in steering scripts, APIs, or higher-level
 * user interfaces that want to enforce explicit values on parameters.
 */

/**
 * @class ParameterSetter
 * @ingroup DataMutationModule
 * @brief Concrete mutator that directly sets the value of parameters.
 *
 * ParameterSetter is a thin adapter on top of the Parameters singleton:
 *  - `mutate()` calls Parameters::setBlockValue(),
 *  - `change_mode()` calls Parameters::changeParameterMode().
 *
 * It additionally enforces some high-level consistency rules, such as
 * preventing changes to EW_SCALE when Wilson coefficients are provided
 * by the user (HAS_WILSON_INPUT).
 */
class ParameterSetter : public IDataMutator<ParamId, scalar_t, ParameterMode> {
public:
    /**
     * @brief Sets the value of a parameter.
     *
     * If the ParamId has no type set, an error is logged since the
     * Parameters instance to target cannot be determined.
     *
     * Special case:
     *  - If the block is `"EW_SCALE"` and the flag
     *    ExternalFlag::HAS_WILSON_INPUT is true in MemoryManager's cache,
     *    the change is refused and a warning is logged.
     *
     * @param pid   The parameter ID.
     * @param value The new central value to set.
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

#endif // PARAMETERMODIFIER_H
