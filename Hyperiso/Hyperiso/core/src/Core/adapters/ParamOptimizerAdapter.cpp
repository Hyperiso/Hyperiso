#include "ParamOptimizerAdapter.h"

ParamOptimizerAdapter::ParamOptimizerAdapter() {

}


void ParamOptimizerAdapter::set_value(const BlockName& block, const LhaID& id, scalar_t v) {
}
void ParamOptimizerAdapter::set_param(const BlockName& block, const LhaID& id, std::shared_ptr<Parameter> p) {
}
void ParamOptimizerAdapter::remove(const BlockName& block, const LhaID& id) {

}

void ParamOptimizerAdapter::commit(bool coalesce = true) {

}

void ParamOptimizerAdapter::clear() { }