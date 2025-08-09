#include "ParameterProvider.h"

scalar_t ParameterProvider::operator()(const ParamId &pid, DataType d_type) {
    if (this->p_type.has_value()) {
        LOG_WARN("LogicError", "This ParameterProvider already has a type.");
    }

    return this->get_value(pid, d_type);
}

scalar_t ParameterProvider::operator()(const std::string &block, const LhaID &id, DataType d_type) const {
    if (!this->p_type.has_value()) {
        LOG_ERROR("LogicError", "Please specify a parameter type for the ParameterProvider.");
    }

    ParamId pid {this->p_type.value(), block, id};
    return this->get_value(pid, d_type);
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

double ParameterProvider::get_scale(const std::string& block) const {
    if (!p_type.has_value()) {
        LOG_ERROR("LogicError", "ParameterProvider has no type.");
    }

    return Parameters::GetInstance(this->p_type.value())->get_block_scale(block);

}

std::shared_ptr<Parameter> ParameterProvider::get_parameter(const ParamId &pid) const {
    return (*Parameters::GetInstance(pid.type.value())).get_parameter(pid.block, pid.code);
}

scalar_t ParameterProvider::get_value(const ParamId& pid, DataType d_type) const {
    switch (d_type) {
    case DataType::VALUE:
        return (*Parameters::GetInstance(pid.type.value()))(pid.block, pid.code);
        break;
    case DataType::STD_STAT:
        return get_parameter(pid)->get_std().first;
        break;
    case DataType::STD_SYST:
        return get_parameter(pid)->get_std().second;
        break;
    case DataType::STD_COMBINED:
        return get_parameter(pid)->get_combined_std();
        break;
    default:
        LOG_ERROR("LogicError", "ParameterProvider does not support this data type.");
    }
}