#include "ParameterRouter.h"

std::vector<std::string> ParameterBlockRepartition::filter_custom_blocks(const std::vector<std::string> &source) {
    std::vector<std::string> custom_blocks {};

    for (const auto& block : source) {
        if (to_lowercase(block) == "mass" || to_lowercase(block) == "gauge") {
            custom_blocks.emplace_back(block);
            continue;
        }

        for (const auto& [_, known_blocks]: ParameterBlockRepartition::BLOCKS) {
            if (std::find(known_blocks.begin(), known_blocks.end(), block) == known_blocks.end()) {
                custom_blocks.emplace_back(block);
            }
        }
    }

    return custom_blocks;
}

ParameterType ParamRouter::GetType(std::string block, LhaID id) {
    auto conflict_blocks = get_keys<std::string, std::unordered_set<long>>(ParametersAccessRights::SM_RIGHTS);

    if (conflict_blocks.contains(block)) {
        switch (MemoryManager::GetInstance()->getMemoryCache().config.model) {
        case Model::THDM:
            if (ParametersAccessRights::THDM_RIGHTS.at(block).contains(id))
                return ParameterType::THDM;
            break;
        case Model::SUSY:
            if (ParametersAccessRights::SUSY_RIGHTS.at(block).contains(id))
                return ParameterType::SUSY;
            break;
        case Model::CUSTOM: // Janky
            if (!ParametersAccessRights::SM_RIGHTS.at(block).contains(id))
                return ParameterType::THDM;
            break;
        default:
            break;
        }
        if (ParametersAccessRights::SM_RIGHTS.at(block).contains(id))
            return ParameterType::SM;
    } else {
        for (auto& [type, owned_blocks] : ParameterBlockRepartition::BLOCKS) {
            if (owned_blocks.contains(block)) {
                return type;
            }
        }
    }

    LOG_ERROR("Invalid Parameter", "Parameter", block, ",", id, "is undefined.");
}

std::vector<ParameterType> ParamRouter::GetType(std::string block) {
    std::vector<ParameterType> param_types;
    for (auto& [type, owned_blocks] : ParameterBlockRepartition::BLOCKS) {
        if (owned_blocks.contains(block)) {
            param_types.emplace_back(type);
        }
    }
    return param_types;
}

std::unordered_set<std::string> ParamRouter::GetOwnedBlocks(ParameterType ptype) {
    return ParameterBlockRepartition::BLOCKS.at(ptype);
}