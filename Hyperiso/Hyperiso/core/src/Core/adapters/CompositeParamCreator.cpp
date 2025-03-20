#include "CompositeParamCreator.h"

void CompositeParamAdapter::add_dependency(
    const std::string &name,
    const std::unordered_map<ParameterType, std::vector<std::string>>
        &source_names,
    ParameterType dest,
    DepUpdateFunc recalculateFunc)
{
    DependentBlockManager::addDependentBlock(name, source_names, dest, recalculateFunc);
}

void CompositeParamAdapter::remove_dependency(const std::string &name, ParameterType src) {
    DependentBlockManager::removeDependentBlock(name, src);
}

void CompositeParamAdapter::update_dependency(const std::string &name, ParameterType src) {
    DependentBlockManager::update(name, src);
}
