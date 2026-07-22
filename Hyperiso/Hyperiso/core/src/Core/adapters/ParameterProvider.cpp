#include "ParameterProvider.h"

ParamId ParameterProvider::resolve_param_id(const ParamId& pid) const {
    ParamId resolved = pid;

    if (this->p_type.has_value()) {
        if (resolved.type.has_value() && resolved.type.value() != this->p_type.value()) {
            LOG_ERROR(
                "LogicError",
                "ParameterProvider type does not match ParamId type."
            );
        }
        resolved.type = this->p_type.value();
    } else if (!resolved.type.has_value()) {
        LOG_ERROR(
            "LogicError",
            "Please specify a parameter type either on ParameterProvider or ParamId."
        );
    }

    return resolved;
}

scalar_t ParameterProvider::operator()(const ParamId &pid, DataType d_type) const {
    return this->get_value(resolve_param_id(pid), d_type);
}

scalar_t ParameterProvider::operator()(const std::string &block, const LhaID &id, DataType d_type) const {
    if (!this->p_type.has_value()) {
        LOG_ERROR("LogicError", "Please specify a parameter type for the ParameterProvider.");
    }

    ParamId pid {this->p_type.value(), block, id};
    return this->get_value(pid, d_type);
}

bool ParameterProvider::exists(const ParamId &pid) const {
    const ParamId resolved = resolve_param_id(pid);
    return Parameters::GetInstance(resolved.type.value())->exist(resolved.block, resolved.code);
}

bool ParameterProvider::exists(const std::string &block, const LhaID &id) const {
    if (!this->p_type.has_value()) {
        LOG_ERROR("LogicError", "Please specify a parameter type for the ParameterProvider.");
    }
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
    const ParamId resolved = resolve_param_id(pid);
    return (*Parameters::GetInstance(resolved.type.value())).get_parameter(resolved.block, resolved.code);
}

scalar_t ParameterProvider::get_value(const ParamId& pid, DataType d_type) const {
    const ParamId resolved = resolve_param_id(pid);
    switch (d_type) {
    case DataType::VALUE:
        return (*Parameters::GetInstance(resolved.type.value()))(resolved.block, resolved.code);
    case DataType::STD_STAT:
        return (*Parameters::GetInstance(resolved.type.value())).get_parameter(resolved.block, resolved.code)->get_std().first;
    case DataType::STD_SYST:
        return (*Parameters::GetInstance(resolved.type.value())).get_parameter(resolved.block, resolved.code)->get_std().second;
    case DataType::STD_COMBINED:
        return (*Parameters::GetInstance(resolved.type.value())).get_parameter(resolved.block, resolved.code)->get_combined_std();
    default:
        LOG_ERROR("LogicError", "ParameterProvider does not support this data type.");
    }
}
