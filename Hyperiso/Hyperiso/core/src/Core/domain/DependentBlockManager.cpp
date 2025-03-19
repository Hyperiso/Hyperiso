#include "DependentBlockManager.h"
 
void DependentBlockManager::addDependentBlock(
    const std::string& name,
    const std::unordered_map<ParameterType, std::vector<std::string>>& source_names,
    ParameterType dest,
    std::function<void(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>)> recalculateFunc
) {
    std::unordered_map<std::string, std::shared_ptr<Block>> sources;
    for (const auto& [k, v] : source_names) {
        for (const auto& src_name : v) {
            sources.emplace(src_name, Parameters::GetInstance(k)->blockAccessor->at(src_name));
        }
    }

    auto dependentBlock = std::make_shared<DependentBlock>(sources, recalculateFunc);
    dependentBlock->init();

    Parameters::GetInstance(dest)->blockAccessor->emplace(name, dependentBlock);
}

void DependentBlockManager::removeDependentBlock(const std::string &name, ParameterType src){
    Parameters::GetInstance(src)->blockAccessor->erase(name);
}
