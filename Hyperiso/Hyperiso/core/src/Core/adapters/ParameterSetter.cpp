#include "ParameterSetter.h"

void ParameterSetter::mutate(const ParamId &pid, scalar_t value) {
    if (!pid.type.has_value()) {
        LOG_ERROR("LogicError", "Use of incomplete ParamId to mutate parameter. Please specify a ParameterType.");
    }

    if (pid.block == "EW_SCALE" && MemoryManager::GetInstance()->getMemoryCache().config.flags.at(ExternalFlag::HAS_WILSON_INPUT)) {
        LOG_WARN("Cannot set matching scale when Wilson coefficients are user input.");
        return;
    }

    Parameters::GetInstance(pid.type.value())->setBlockValue(pid.block, pid.code, value);
}

void ParameterSetter::change_mode(const ParamId& pid, ParameterMode mode) {
    if (!pid.type.has_value()) {
        LOG_ERROR("LogicError", "Use of incomplete ParamId to mutate parameter. Please specify a ParameterType.");
    }

    Parameters::GetInstance(pid.type.value())->changeParameterMode(pid, mode);
}