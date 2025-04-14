#include "ModelParamAdapter.h"

ParameterProxy::ParameterProxy(ParameterType type) { 
    if (!ParameterProxy::ALLOWED.contains(type)) {
        LOG_ERROR("ValueError", "PhysicalModel cannot access parameter type", ParameterTypeMapper::str(type));
    }
    this->pp = ParameterProvider(type);
} 

double ParameterProxy::operator()(const std::string& block, const LhaID& id) const { 
    if (pp.get_type() == ParameterType::WILSON) {
        std::cout << block << " : "<< id << std::endl;
        return pp.exists(block, id) ? pp(block, id) : 0.;
    } 
    std::cout << "v2" << block << id << std::endl;
    return pp(block, id); 
};