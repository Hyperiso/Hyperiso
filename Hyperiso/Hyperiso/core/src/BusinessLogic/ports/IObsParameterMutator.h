#ifndef IOBS_PARAMETER_MUTATOR_H
#define IOBS_PARAMETER_MUTATOR_H

#include "Include.h"
#include "Math.h"

/**
 * @file IObsParameterMutator.h
 * @brief Interface for mutating parameters from the observable/business layer.
 *
 * This interface abstracts how parameters are updated:
 *  - set: absolute assignment
 *  - shift: relative update
 *  - change_mode: change statistical/scan mode for the parameter
 *
 * Template parameters:
 *  - T: parameter identifier type (e.g. ParamId)
 *  - U: parameter value type (e.g. scalar_t)
 *  - V: mode enum type (e.g. ParameterMode)
 *
 * @see ObsParameterMutator
 */
template <typename T, typename U, typename V>
class IObsParameterMutator {
public:
    virtual ~IObsParameterMutator() = default;
    
    /// Assign an absolute value.
    virtual void set(const T&, U) = 0;

    /// Apply a relative shift.
    virtual void shift(const T&, U) = 0;

    /// Change the mode (fixed/nuisance/scan...).
    virtual void change_mode(const T&, V) = 0;
};

#endif // IOBS_PARAMETER_MUTATOR_H
