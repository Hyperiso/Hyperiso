#ifndef PARAMETERSHIFTER_H
#define PARAMETERSHIFTER_H

#include "IDataMutator.h"
#include "Parameters.h"

/**
 * @class ParameterShifter
 * @ingroup DataMutationModule
 * @brief Concrete mutator that shifts parameter values rather than directly setting them.
 */
class ParameterShifter : public IDataMutator<ParamId, scalar_t, ParameterMode> {
public:

    /**
     * @brief Shifts the value of a parameter by a given amount.
     * @param pid The parameter ID.
     * @param value The amount by which to shift the parameter value.
     */
    void mutate(const ParamId& pid, scalar_t value) override;

    /**
     * @brief Changes the mode of a parameter (e.g., from FIXED to SHIFTABLE).
     * @param pid The parameter ID.
     * @param mode The new parameter mode.
     */
    void change_mode(const ParamId& pid, ParameterMode mode) override;
};


#endif // PARAMETERSHIFTER_H
