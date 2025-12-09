#ifndef PARAMETERMODIFIER_H
#define PARAMETERMODIFIER_H

#include "IDataMutator.h"
#include "Parameters.h"
#include "MemoryManager.h"
#include "Math.h"

/**
 * @class ParameterSetter
 * @ingroup DataMutationModule
 * @brief Concrete mutator that directly sets the value of parameters.
 */
class ParameterSetter : public IDataMutator<ParamId, scalar_t, ParameterMode> {
public:

    /**
     * @brief Sets the value of a parameter.
     * @param pid The parameter ID.
     * @param value The new value to set.
     */
    void mutate(const ParamId& pid, scalar_t value) override;

    /**
     * @brief Changes the mode of a parameter (e.g., from FIXED to SHIFTABLE).
     * @param pid The parameter ID.
     * @param mode The new parameter mode.
     */
    void change_mode(const ParamId& pid, ParameterMode mode) override;
};

#endif // PARAMETERMODIFIER_H
