#include "ObsParameterProxy.h"

ObsParameterProxy::ObsParameterProxy(ParameterType type) { 
    if (!ObsParameterProxy::ALLOWED.contains(type)) {
        LOG_ERROR("ValueError", "PhysicalModel cannot access parameter type", ParameterTypeMapper::str(type));
    }
    this->pp_with_type = ParameterProvider(type);
    this->pp = ParameterProvider();
} 

scalar_t ObsParameterProxy::operator()(const std::string& block, const LhaID& id) const { 
    if (pp.get_type() == ParameterType::WILSON) {
        return pp.exists(block, id) ? pp(block, id) : scalar_t();
    } 
    return pp(block, id); 
};

scalar_t ObsParameterProxy::operator()(const ParamId& pid, ParameterProvider::DataType d_type=ParameterProvider::DataType::VALUE) { 

    if (!pid.type.has_value()) {
        LOG_WARN("LogicError", "Use of untyped ParamId in ParameterProvider.");
    }
    if (pp.get_type() == ParameterType::WILSON) {
        return pp.exists(pid) ? pp(pid) : scalar_t();
    }
    pp(pid, d_type);

    return pp(pid); 
};

std::shared_ptr<Parameter> ObsParameterProxy::get_parameter(const ParamId& pid) const{
    return pp.get_parameter(pid);
}