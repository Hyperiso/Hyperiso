#include "ModelParamAdapter.h"

ParameterProxy::ParameterProxy(ParameterType type) { 
    if (!ParameterProxy::ALLOWED.contains(type)) {
        LOG_ERROR("ValueError", "PhysicalModel cannot access parameter type", ParameterTypeMapper::str(type));
    }
    std::cout << "fuuuuck proxy" << std::endl;
    this->pp = ParameterProvider(type);
    std::cout << "end fuuuuck proxy" << std::endl;
} 

scalar_t ParameterProxy::operator()(const std::string& block, const LhaID& id) const { 
    if (pp.get_type() == ParameterType::WILSON) {
        std::cout << block << " : "<< id << std::endl;
        return pp.exists(block, id) ? pp(block, id) : scalar_t();
    } 
    std::cout << "v2" << block << id << std::endl;
    return pp(block, id); 
};