#include "APIAdapter.h"

std::map<LhaID, double> APIAdapter::get_block_infos(const std::string& block, ParameterType param_type) {
    return Parameters::GetInstance(param_type)->get_block_infos(block);
}

std::unordered_set<std::string> APIAdapter::get_blocks_list(ParameterType param_type) {
    return Parameters::GetInstance(param_type)->get_blocks_list();
}

bool APIAdapter::check_flag(ExternalFlag flag) {
    return MemoryManager::GetInstance()->getMemoryCache().config.flags.at(flag);
}

fs::path APIAdapter::get_path(APIPath path_name) {
    switch (path_name) {
    case APIPath::LHA_PATH:
        return MemoryManager::GetInstance()->getMemoryCache().lha_path;
        break;
    default:
        LOG_ERROR("ValueError", "Unknown path for APIAdapter.");
    };
}

std::unordered_set<std::string> APIAdapter::get_all_blocks()
{
    std::unordered_set<std::string> all_blocks;
    for (auto p_type : MemoryManager::GetInstance()->getMemoryCache().parameter_types) {
        auto blocks = get_blocks_list(p_type);
        all_blocks.insert(blocks.begin(), blocks.end());
    }   
    return all_blocks;
}

std::vector<ParameterType> APIAdapter::get_type_of_block(const std::string& block) {
    std::vector<ParameterType> param_type;
    for (auto& elem : MemoryManager::GetInstance()->getMemoryCache().parameter_types) {
        for (auto& block_ : get_blocks_list(elem)) {
            if (block == block_) {
                param_type.push_back(elem);
                break;
            }
        }
    }
    return param_type;
}