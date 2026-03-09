#include "StatParamOptimizerProxy.h"

StatParamOptimizerProxy::StatParamOptimizerProxy() : poa({ParameterType::SM, ParameterType::FLAVOR, ParameterType::DECAY, ParameterType::WILSON}) { //TODO : add the right parameterType
    
}


void StatParamOptimizerProxy::set_value(const BlockName& block, const LhaID& id, scalar_t v) {
    poa.set_value(block, id, v);
}
void StatParamOptimizerProxy::set_param(const BlockName& block, const LhaID& id, std::shared_ptr<Parameter> p) {
    poa.set_param(block,id, p);
}
void StatParamOptimizerProxy::remove(const BlockName& block, const LhaID& id) {
    poa.remove(block, id);
}

void StatParamOptimizerProxy::commit(bool coalesce) {
    poa.commit(coalesce);
}

void StatParamOptimizerProxy::clear() {
    poa.clear();
}