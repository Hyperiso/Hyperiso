#include "ParameterSetter.h"

void ParameterSetter::mutate(const ParamId &pid, scalar_t value) {
    if (!pid.type.has_value()) {
        LOG_ERROR("LogicError", "Use of incomplete ParamId to mutate parameter. Please specify a ParameterType.");
    }

    Parameters::GetInstance(pid.type.value())->setBlockValue(pid.block, pid.code, value);
}

void ParameterSetter::change_mode(const ParamId& pid, ParameterMode mode) {
    if (!pid.type.has_value()) {
        LOG_ERROR("LogicError", "Use of incomplete ParamId to mutate parameter. Please specify a ParameterType.");
    }

    Parameters::GetInstance(pid.type.value())->changeParameterMode(pid, mode);
}