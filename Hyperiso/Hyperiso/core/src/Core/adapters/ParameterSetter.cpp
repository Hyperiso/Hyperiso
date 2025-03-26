#include "ParameterSetter.h"

void ParameterSetter::mutate(const ParamId &pid, double value) {
    if (!pid.type.has_value()) {
        LOG_ERROR("LogicError", "Use of incomplete ParamId to mutate parameter. Please specify a ParameterType.");
    }

    Parameters::GetInstance(pid.type.value())->get_parameter(pid.block, pid.code)->set_expected(value);
}