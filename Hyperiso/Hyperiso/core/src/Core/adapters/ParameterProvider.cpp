#include "ParameterProvider.h"

double ParameterProvider::operator()(const ParamId &pid, DataType d_type) {
    if (this->p_type.has_value()) {
        LOG_WARN("LogicError", "This ParameterProvider already has a type.");
    }

    return (*Parameters::GetInstance(pid.type))(pid.block, pid.code);
}

double ParameterProvider::operator()(const std::string &block, const LhaID &id, DataType d_type) {
    if (!this->p_type.has_value()) {
        LOG_ERROR("LogicError", "Please specify a parameter type for the ParameterProvider.");
    }

    return (*Parameters::GetInstance(this->p_type.value()))(block, id);
}
