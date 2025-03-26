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
    dependentBlock->blockname = name;
    dependentBlock->init();
    dependentBlock->update();
    Parameters::GetInstance(dest)->blockAccessor->emplace(name, dependentBlock);
}

void DependentBlockManager::addDependentParameter(
    const ParamId &pid,
    const std::unordered_set<ParamId> &source_pids,
    DepParamUpdateFunc recalculateFunc)
{
    std::unordered_map<ParamId, std::shared_ptr<Parameter>> sources;
    
    for (const auto& id : source_pids) {
        sources.emplace(id, Parameters::GetInstance(id.type.value())->blockAccessor->getParameter(id.block, id.code));
    }
    
    auto dependentParam = std::make_shared<DependentParameter>(sources, recalculateFunc);
    dependentParam->init();
    dependentParam->update();
    auto ba = Parameters::GetInstance(pid.type.value())->blockAccessor;

    if (!ba->contains(pid.block)) {
        ba->emplace(pid.block, std::make_shared<Block>());
        ba->at(pid.block)->blockname = pid.block;
    } 
    
    ba->at(pid.block)->store(pid.code, dependentParam);
}

void DependentBlockManager::removeDependentBlock(const std::string &name,
                                                 ParameterType src)
{
    Parameters::GetInstance(src)->blockAccessor->erase(name);
}

void DependentBlockManager::update(const std::string &name, ParameterType src) {
    Parameters::GetInstance(src)->blockAccessor->at(name)->update();
}
