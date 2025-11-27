#include "ParamOptimizerAdapter.h"
#include "Parameters.h"

ParamOptimizerAdapter::ParamOptimizerAdapter(std::vector<ParameterType> scopes) {
    std::vector<std::shared_ptr<BlockAccessor>> input;
    for (auto scope : scopes) {
        input.push_back(Parameters::GetInstance(scope)->get_block_accessor());
    }
    this->po = std::make_shared<ParamOptimizer>(input);
}


void ParamOptimizerAdapter::set_value(const BlockName& block, const LhaID& id, scalar_t v) {
    po->set_value(block, id, v);
}
void ParamOptimizerAdapter::set_param(const BlockName& block, const LhaID& id, std::shared_ptr<Parameter> p) {
    po->set_param(block,id, p);
}
void ParamOptimizerAdapter::remove(const BlockName& block, const LhaID& id) {
    po->remove(block, id);
}

void ParamOptimizerAdapter::commit(bool coalesce) {
    po->commit(coalesce);
}

void ParamOptimizerAdapter::clear() {
    po->clear();
}