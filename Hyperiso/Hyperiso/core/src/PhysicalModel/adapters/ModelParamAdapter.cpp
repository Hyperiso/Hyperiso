#include "ModelParamAdapter.h"

ParameterProxy::ParameterProxy(ParameterType type) { 
    if (!ParameterProxy::ALLOWED.contains(type)) {
        LOG_ERROR("ValueError", "PhysicalModel cannot access parameter type", ParameterTypeMapper::str(type));
    }
    this->pp = ParameterProvider(type);
} 

scalar_t ParameterProxy::operator()(const std::string& block, const LhaID& id) const { 
    if (pp.get_type() == ParameterType::WILSON) {
        return pp.exists(block, id) ? pp(block, id) : scalar_t();
    } 
    return pp(block, id); 
};