#include "ParamSourcesProvider.h"
#include "MemoryManager.h"

std::unordered_set<ParamId> ParamSourcesProvider::get_all_leaf_sources(const std::unordered_set<ParamId>& param_ids) const {
    return MemoryManager::GetInstance()->get_all_source_parameters(param_ids);
}