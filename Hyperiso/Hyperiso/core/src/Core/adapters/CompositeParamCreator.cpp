#include "CompositeParamCreator.h"

void CompositeParamCreator::add_dependency(
    const std::string &name,
    const std::unordered_map<ParameterType, std::vector<std::string>>
        &source_names,
    ParameterType dest,
    std::function<
        void(const std::unordered_map<std::string, std::shared_ptr<Block>> &,
             std::shared_ptr<DependentBlock>)> recalculateFunc)
{
    DependentBlockManager::addDependentBlock(name, source_names, dest, recalculateFunc);
}

void CompositeParamCreator::remove_dependency(const std::string &name, ParameterType src) {
    DependentBlockManager::removeDependentBlock(name, src);
}
