#include "CompositeParamCreator.h"

void CompositeParamAdapter::add_block_dependency(
    const std::string &name,
    const std::unordered_map<ParameterType, std::vector<std::string>>
        &source_names,
    ParameterType dest,
    DepUpdateFunc recalculateFunc)
{
    DependentBlockManager::addDependentBlock(name, source_names, dest, recalculateFunc);
}

void CompositeParamAdapter::add_param_dependency(
    const ParamId &pid,
    const std::unordered_set<ParamId> &source_pids,
    DepParamUpdateFunc recalculateFunc) 
{
    DependentBlockManager::addDependentParameter(pid, source_pids, recalculateFunc);
}

void CompositeParamAdapter::remove_dependency(const std::string &name, ParameterType src) {
    DependentBlockManager::removeDependentBlock(name, src);
}

void CompositeParamAdapter::update_dependency(const std::string &name, ParameterType src) {
    DependentBlockManager::update(name, src);
}
