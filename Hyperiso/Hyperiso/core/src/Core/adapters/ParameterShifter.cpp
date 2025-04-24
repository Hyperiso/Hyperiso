#include "ParameterShifter.h"

void ParameterShifter::mutate(const ParamId &pid, scalar_t value) {
    if (!pid.type.has_value()) {
        LOG_ERROR("LogicError", "Use of incomplete ParamId to mutate parameter. Please specify a ParameterType.");
    }

    Parameters::GetInstance(pid.type.value())->shiftParameter(pid, value);
}