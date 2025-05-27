#include "DependentBlockManager.h"
#include "Parameters.h"

void DependentBlockManager::addDependentBlock(
    std::string name,
    std::unordered_map<ParameterType, std::vector<std::string>> source_names,
    ParameterType dest,
    std::function<void(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>)> recalculateFunc
) {
    std::unordered_map<std::string, std::shared_ptr<Block>> sources;

    for (const auto& [k, v] : source_names) {
        for (const auto& src_name : v) {
            auto it = std::find_if(
                Parameters::GetInstance(k)->blockAccessor->begin(),
                Parameters::GetInstance(k)->blockAccessor->end(),
                [&](const auto& pair) { return pair.first == src_name; }
            );

            if (it == Parameters::GetInstance(k)->blockAccessor->end()) {
                std::cout << "k: " << ParameterTypeMapper::str(k) << std::endl;
                std::cout << "bad: " << src_name << std::endl;
                std::cout << Parameters::GetInstance(k)->blockAccessor << std::endl;
            }
            sources.emplace(src_name, Parameters::GetInstance(k)->blockAccessor->at(src_name));
        }
    }
    auto dependentBlock = std::make_shared<DependentBlock>(sources, recalculateFunc);
    dependentBlock->blockname = name;
    dependentBlock->init();
    dependentBlock->update();
    Parameters::GetInstance(dest)->blockAccessor->emplace(name, dependentBlock);
    // std::cout << "--------------------------------------------------------------------------" << std::endl;
    // std::cout << Parameters::GetInstance(dest)->blockAccessor << std::endl;
}

void DependentBlockManager::addDependentParameter(
    ParamId pid,
    std::unordered_set<ParamId> source_pids,
    DepParamUpdateFunc recalculateFunc)
{
    std::unordered_map<ParamId, std::shared_ptr<Parameter>> sources;
    
    for (const auto& id : source_pids) {
        sources.emplace(id, Parameters::GetInstance(id.type.value())->blockAccessor->getParameter(id.block, id.code));
    }
    
    auto dependentParam = std::make_shared<DependentParameter>(pid, sources, recalculateFunc);
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
    if (!Parameters::GetInstance(src)->blockAccessor->contains(name)) {
        return;
    }

    std::shared_ptr<DependentBlock> dep_block = std::static_pointer_cast<DependentBlock>(Parameters::GetInstance(src)->blockAccessor->at(name));
    for (auto& src : dep_block->get_source_blocks()) {
        src.second->removeObserver(dep_block);
    }

    for (auto& observer : dep_block->getObservers()) {
        removeDependentBlock(observer->get_name(), src);
    }

    Parameters::GetInstance(src)->blockAccessor->erase(name);
}

void DependentBlockManager::update(const std::string &name, ParameterType src) {
    Parameters::GetInstance(src)->blockAccessor->at(name)->update();
}
