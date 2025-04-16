#include "MartyParameterProxy.h"

MartyParameterProxy::MartyParameterProxy(ParameterType type) { 
    if (!MartyParameterProxy::ALLOWED.contains(type)) {
        LOG_ERROR("ValueError", "MartyInterface cannot access parameter type", ParameterTypeMapper::str(type));
    }
    this->pp = ParameterProvider(type);
} 

scalar_t MartyParameterProxy::operator()(const std::string& block, const LhaID& id) const { 
    return pp(block, id); 
};