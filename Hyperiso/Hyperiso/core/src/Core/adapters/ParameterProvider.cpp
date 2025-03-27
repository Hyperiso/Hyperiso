#include "ParameterProvider.h"

double ParameterProvider::operator()(const ParamId &pid, DataType d_type) {
    if (this->p_type.has_value()) {
        LOG_WARN("LogicError", "This ParameterProvider already has a type.");
    }

    if (!pid.type.has_value()) {
        LOG_WARN("LogicError", "Use of untyped ParamId in ParameterProvider.");
    }

    return (*Parameters::GetInstance(pid.type.value()))(pid.block, pid.code);
}

double ParameterProvider::operator()(const std::string &block, const LhaID &id, DataType d_type) const {
    if (!this->p_type.has_value()) {
        LOG_ERROR("LogicError", "Please specify a parameter type for the ParameterProvider.");
    }

    return (*Parameters::GetInstance(this->p_type.value()))(block, id);
}

bool ParameterProvider::exists(const ParamId &pid) const {
    return Parameters::GetInstance(pid.type.value())->exist(pid.block, pid.code);
}

bool ParameterProvider::exists(const std::string &block, const LhaID &id) const {
    return Parameters::GetInstance(this->p_type.value())->exist(block, id);
}

ParameterType ParameterProvider::get_type() const {
    if (!p_type.has_value()) {
        LOG_ERROR("LogicError", "ParameterProvider has no type.");
    }

    return p_type.value();
}
