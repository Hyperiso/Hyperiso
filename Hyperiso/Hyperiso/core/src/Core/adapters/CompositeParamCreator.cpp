#include "CompositeParamCreator.h"

void CompositeParamAdapter::add_block_dependency(
    BlockName name,
    std::unordered_map<ParameterType, std::vector<std::string>> source_names,
    ParameterType dest,
    DepUpdateFunc recalculateFunc)
{
    DependentBlockManager::addDependentBlock(name, source_names, dest, recalculateFunc);
}

void CompositeParamAdapter::add_param_dependency(
    ParamId pid,
    std::unordered_set<ParamId> source_pids,
    DepParamUpdateFunc recalculateFunc) 
{
    DependentBlockManager::addDependentParameter(pid, source_pids, recalculateFunc);
}

void CompositeParamAdapter::remove_dependency(const BlockName &name, ParameterType src) {
    DependentBlockManager::removeDependentBlock(name, src);
}

void CompositeParamAdapter::update_dependency(const BlockName &name, ParameterType src) {
    DependentBlockManager::update(name, src);
}
