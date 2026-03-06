#include "DependentBlockManager.h"
#include "Parameters.h"

void DependentBlockManager::addDependentBlock(
    std::string name,
    std::unordered_map<ParameterType, std::vector<std::string>> source_names,
    ParameterType dest,
    DepUpdateFunc recalculateFunc   // <-- utilise le nouveau typedef
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
            for (auto& bn : ba->get_block_names()) { if (!first) oss << ", "; first = false; oss << bn; }  // :contentReference[oaicite:2]{index=2}
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
        // if (!ba->contains(id.block)) {
        //     missing.push_back("Block " + ParameterTypeMapper::str(id.type.value()) + "::" + id.block);
        //     continue;
        // }
        // auto blk = ba->at(id.block);
        // if (!blk->contains(id.code)) {
        //     missing.push_back("Param " + ParameterTypeMapper::str(id.type.value()) + "::" + id.block + "::" + id.code.to_string());
        //     continue;
        // }
        // sources.emplace(id, blk->retrieve(id.code));

        auto blk = ba->at(id.block);

        try {
            // retrieve() doit déclencher le lazy ensure_up_to_date() côté block
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

    // auto dependentParam = std::make_shared<DependentParameter>(pid, sources, recalculateFunc);
    // dependentParam->init();
    // dependentParam->update();

    auto ba = Parameters::GetInstance(pid.type.value())->blockAccessor;
    if (!ba->contains(pid.block)) {
        ba->emplace(pid.block, std::make_shared<Block>());
        ba->at(pid.block)->blockname = pid.block;
    }
    // std::cout << "setting dep param : " << pid.block << " : " << pid.code << " = "<< *dependentParam << std::endl;
    // ba->at(pid.block)->store_or_assign(pid.code, dependentParam);
    auto blk = ba->at(pid.block);

    if (blk->contains(pid.code)) {
        auto existing = blk->retrieve(pid.code);
        if (auto dep = std::dynamic_pointer_cast<DependentParameter>(existing)) {
            dep->rebind(std::move(sources), recalculateFunc);
            dep->init();
            // pas besoin de dep->update() : rebind met dirty + notify
            return;
        }
        // sinon: cas dangereux (placeholder non-dependent). On peut soit throw, soit overwrite payload.
        // throw std::logic_error("Expected placeholder DependentParameter for " + pid.code.to_string());
    }

    // sinon il n'existe pas => création normale
    auto dependentParam = std::make_shared<DependentParameter>(pid, std::move(sources), recalculateFunc);
    dependentParam->init();
    dependentParam->update();
    blk->store(pid.code, dependentParam);
}

void DependentBlockManager::removeDependentBlock(const std::string &name,
                                                 ParameterType src)
{
    if (!Parameters::GetInstance(src)->blockAccessor->contains(name)) {
        return;
    }
    std::shared_ptr<Block> dep_block = Parameters::GetInstance(src)->blockAccessor->at(name);
    dep_block->clear_above();


    Parameters::GetInstance(src)->blockAccessor->erase(name);
}

void DependentBlockManager::update(const std::string &name, ParameterType src) {
    Parameters::GetInstance(src)->blockAccessor->at(name)->update();
}
