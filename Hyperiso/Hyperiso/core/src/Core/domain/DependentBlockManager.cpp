#include "DependentBlockManager.h"
#include "Parameters.h"

void DependentBlockManager::addDependentBlock(
    std::string name,
    std::unordered_map<ParameterType, std::vector<std::string>> source_names,
    ParameterType dest,
    DepUpdateFunc recalculateFunc 
) {
    std::unordered_map<std::string, std::shared_ptr<Block>> sources;
    std::vector<std::string> missing;

    for (const auto& [k, v] : source_names) {
        auto ba = Parameters::GetInstance(k)->blockAccessor;
        for (const auto& src_name : v) {
            if (!ba->contains(src_name)) {
                missing.push_back(ParameterTypeMapper::str(k) + "::" + src_name);
                continue;
            }
            sources.emplace(src_name, ba->at(src_name));
        }
    }

    if (!missing.empty()) {
        std::ostringstream oss;
        oss << "addDependentBlock '" << name << "': missing source blocks: ";
        for (size_t i = 0; i < missing.size(); ++i) { if (i) oss << ", "; oss << missing[i]; }
        oss << ". Available per model:\n";
        for (const auto& [k, _] : source_names) {
            auto ba = Parameters::GetInstance(k)->blockAccessor;
            oss << "  - " << ParameterTypeMapper::str(k) << ": ";
            bool first = true;
            for (auto& bn : ba->get_block_names()) { if (!first) oss << ", "; first = false; oss << bn; }
            oss << "\n";
        }
        throw std::invalid_argument(oss.str());
    }

    auto dependentBlock = std::make_shared<DependentBlock>(sources, recalculateFunc);
    dependentBlock->blockname = name;
    Parameters::GetInstance(dest)->blockAccessor->emplace(name, dependentBlock);
    dependentBlock->init();
    dependentBlock->update();
}

void DependentBlockManager::addDependentParameter(
    ParamId pid,
    std::unordered_set<ParamId> source_pids,
    DepParamUpdateFunc recalculateFunc 
) {
    std::unordered_map<ParamId, std::shared_ptr<Parameter>> sources;
    std::vector<std::string> missing;

    for (const auto& id : source_pids) {
        auto ba = Parameters::GetInstance(id.type.value())->blockAccessor;

        auto blk = ba->at(id.block);

        try {
            auto p = blk->retrieve(id.code);
            sources.emplace(id, p);
        } catch (...) {
            missing.push_back("Param " + ParameterTypeMapper::str(id.type.value()) +
                            "::" + id.block + "::" + id.code.to_string());
            continue;
        }
    }

    if (!missing.empty()) {
        std::ostringstream oss;
        oss << "addDependentParameter '" 
            << ParameterTypeMapper::str(pid.type.value()) << "::" << pid.block << "::" << pid.code.to_string()
            << "': missing sources: ";
        for (size_t i = 0; i < missing.size(); ++i) { if (i) oss << ", "; oss << missing[i]; }
        throw std::invalid_argument(oss.str());
    }

    auto ba = Parameters::GetInstance(pid.type.value())->blockAccessor;
    if (!ba->contains(pid.block)) {
        ba->emplace(pid.block, std::make_shared<Block>());
        ba->at(pid.block)->blockname = pid.block;
    }

    auto blk = ba->at(pid.block);

    if (blk->contains(pid.code)) {
        auto existing = blk->retrieve(pid.code);
        if (auto dep = std::dynamic_pointer_cast<DependentParameter>(existing)) {
            dep->rebind(std::move(sources), recalculateFunc);
            dep->init();
            return;
        }

    }

    auto dependentParam = std::make_shared<DependentParameter>(pid, std::move(sources), recalculateFunc);
    dependentParam->init();
    dependentParam->update();
    blk->store(pid.code, dependentParam);
}

void DependentBlockManager::removeDependentBlock(const std::string &name,
                                                 ParameterType src)
{
    auto params = Parameters::GetInstance(src);
    auto accessor = params->blockAccessor;

    if (!accessor->contains(name)) {
        return;
    }

    std::shared_ptr<Block> dep_block = accessor->at(name);
    if (dep_block) {
        dep_block->destroy();
    }

    accessor->erase_block(name);
}

void DependentBlockManager::update(const std::string &name, ParameterType src) {
    Parameters::GetInstance(src)->blockAccessor->at(name)->update();
}
