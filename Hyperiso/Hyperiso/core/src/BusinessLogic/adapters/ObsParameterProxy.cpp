#include "ObsParameterProxy.h"

ObsParameterProxy::ObsParameterProxy(ParameterType type) { 
    if (!ObsParameterProxy::ALLOWED.contains(type)) {
        LOG_ERROR("ValueError", "PhysicalModel cannot access parameter type", ParameterTypeMapper::str(type));
    }
    this->pp_with_type = ParameterProvider(type);
    this->pp = ParameterProvider();
} 

scalar_t ObsParameterProxy::operator()(const std::string& block, const LhaID& id, ParameterProvider::DataType d_type) const {
    LOG_INFO(block, id);
    if (pp_with_type.get_type() == ParameterType::WILSON) {
        return pp_with_type.exists(block, id) ? pp_with_type(block, id, d_type) : scalar_t();
    } 
    scalar_t value = pp_with_type(block, id, d_type);
    LOG_INFO("coucou");
    return value;
};

scalar_t ObsParameterProxy::operator()(const ParamId& pid, ParameterProvider::DataType d_type) { 

    if (!pid.type.has_value()) {
        LOG_WARN("LogicError", "Use of untyped ParamId in ParameterProvider.");
    }

    return pp(pid, d_type); 
};

std::shared_ptr<Parameter> ObsParameterProxy::get_parameter(const ParamId& pid) const{
    return pp.get_parameter(pid);
}