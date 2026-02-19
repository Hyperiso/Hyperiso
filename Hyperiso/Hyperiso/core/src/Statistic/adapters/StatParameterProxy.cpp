#include "StatParameterProxy.h"

std::shared_ptr<Parameter> StatParameterProxy::get_param(const ParamId& pid) const {
    return pp.get_parameter(pid);
}
std::shared_ptr<Parameter> StatParameterProxy::get_param(const std::string& block, const LhaID& id) const {
    return pp_with_type.get_parameter(ParamId(pp_with_type.get_type(), block, id));
}

scalar_t StatParameterProxy::operator()(const ParamId& pid, DataType d_type) const {
    if (!pid.type.has_value()) {
        LOG_WARN("LogicError", "Use of untyped ParamId in ParameterProvider.");
    }

    if (!StatParameterProxy::ALLOWED.contains(pid.type.value())) {
        LOG_ERROR("ValueError", "BusinessLogic cannot access parameter type", ParameterTypeMapper::str(pid.type.value()));
    }
    if (pid.type == ParameterType::WILSON) {
        return pp.exists(pid) ? pp(pid, d_type) : scalar_t();
    } 
    return pp(pid, d_type); 

}
double StatParameterProxy::operator()(const ObservableId& id, DataType d_type) const {
    return pp(ParamId(ParameterType::OBSERVABLE, "FOBS", ObservableMapper::flha(id)), d_type);
}

double StatParameterProxy::operator()(const BinnedObservableId &id, DataType d_type) const {
    return pp(ParamId(ParameterType::OBSERVABLE, "FOBS", id.flha()), d_type);
}

scalar_t StatParameterProxy::operator()(const std::string& block, const LhaID& id, DataType d_type) const {
    if (pp_with_type.get_type() == ParameterType::WILSON) {
        return pp_with_type.exists(block, id) ? pp_with_type(block, id, d_type) : scalar_t();
    } 
    scalar_t value = pp_with_type(block, id, d_type);
    return value;
}
std::shared_ptr<Parameter> StatParameterProxy::get_obs_param(const BinnedObservableId& id) const {
    return pp.get_parameter(ParamId(ParameterType::OBSERVABLE, "FOBS", id.flha()));
}

StatParameterProxy::StatParameterProxy(ParameterType type) { 
    if (!StatParameterProxy::ALLOWED.contains(type)) {
        LOG_ERROR("ValueError", "BusinessLogic cannot access parameter type", ParameterTypeMapper::str(type));
    }
    this->pp_with_type = ParameterProvider(type);
    this->pp = ParameterProvider();
} 