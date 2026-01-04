#ifndef OBSPARAMETERMUTATOR_H
#define OBSPARAMETERMUTATOR_H

#include "Include.h"
#include "Math.h"
#include "IObsParameterMutator.h"
#include "ParameterSetter.h"
#include "ParameterShifter.h"

/**
 * @file ObsParameterMutator.h
 * @brief High-level mutator for parameters used by the observable layer.
 *
 * This class provides a unified interface to:
 *  - set a parameter value absolutely (@ref set),
 *  - apply a shift relative to the current value (@ref shift),
 *  - switch the parameter mode (value / gaussian nuisance / flat scan, etc.) (@ref change_mode).
 *
 * It delegates the work to:
 *  - @ref ParameterSetter  for absolute assignments
 *  - @ref ParameterShifter for relative shifts and mode changes
 *
 * The identifier type is @ref ParamId and the stored value type is @ref scalar_t.
 *
 * Typical usage:
 * @code
 *   ObsParameterMutator mut;
 *   mut.set(ParamId{ParameterType::SM, "MASS", 5}, 4.18);
 *   mut.shift(ParamId{ParameterType::SM, "MASS", 5}, 0.02);
 *   mut.change_mode(ParamId{ParameterType::SM, "MASS", 5}, ParameterMode::GAUSSIAN);
 * @endcode
 *
 * @see IObsParameterMutator
 * @see ParameterSetter
 * @see ParameterShifter
 */
class ObsParameterMutator: public IObsParameterMutator<ParamId, scalar_t, ParameterMode> {
public:
    /**
     * @brief Assigns an absolute value to a parameter.
     * @param pid   Parameter identifier.
     * @param value New value.
     */
    void set(const ParamId&, scalar_t) override;

    /**
     * @brief Shifts a parameter value (relative update).
     * @param pid   Parameter identifier.
     * @param value Shift amount (interpretation depends on ParameterShifter policy).
     */
    void shift(const ParamId&, scalar_t) override;

    /**
     * @brief Changes the “mode” of a parameter (e.g. fixed, nuisance, scan...).
     * @param pid Parameter identifier.
     * @param mod New mode.
     */
    void change_mode(const ParamId&, ParameterMode) override;

private:
    /// Backend for absolute assignments.
    ParameterSetter p_set;

    /// Backend for relative shifts and mode management.
    ParameterShifter p_shift;
};

#endif // OBSPARAMETERMUTATOR_H