#include "ModelParamAdapter.h"

ParameterProxy::ParameterProxy(ParameterType type) { 
    if (!ParameterProxy::ALLOWED.contains(type)) {
        LOG_ERROR("ValueError", "PhysicalModel cannot access parameter type", ParameterTypeMapper::str(type));
    }
    this->pp = ParameterProvider(type);
} 

double ParameterProxy::operator()(const std::string& block, const LhaID& id) const { 
    return pp(block, id); 
};