#include "ParameterRouter.h"

std::vector<BlockName> ParameterBlockRepartition::filter_custom_blocks(const std::vector<BlockName> &source) {
    std::unordered_set<std::string> known_norm;
    for (const auto& [_, known] : ParameterBlockRepartition::BLOCKS) {
        for (const auto& b : known) {
            known_norm.insert(to_lowercase(b));
        }
    }

    std::vector<BlockName> custom_blocks;
    custom_blocks.reserve(source.size());

    for (const auto& block : source) {
        const auto bnorm = to_lowercase(block);

        // Special cases
        if (bnorm == "mass" || bnorm == "gauge") {
            custom_blocks.emplace_back(block);
            continue;
        }

        if (!known_norm.count(bnorm)) {
            custom_blocks.emplace_back(block);
        }
    }

    return custom_blocks;
}

ParameterType ParamRouter::GetType(BlockName block, LhaID id) {
    auto conflict_blocks = get_keys<BlockName, std::unordered_set<long>>(ParametersAccessRights::SM_RIGHTS);

    if (conflict_blocks.contains(block)) {
        return ParametersAccessRights::SM_RIGHTS.at(block).contains(id) ? ParameterType::SM : ParameterType::BSM;
    } else {
        for (auto& [type, owned_blocks] : ParameterBlockRepartition::BLOCKS) {
            if (owned_blocks.contains(block)) {
                return type;
            }
        }
    }

    LOG_ERROR("Invalid Parameter", "Parameter", block, ",", id, "is undefined.");
}

std::vector<ParameterType> ParamRouter::GetType(BlockName block) {
    std::vector<ParameterType> param_types;
    for (auto& [type, owned_blocks] : ParameterBlockRepartition::BLOCKS) {
        if (owned_blocks.contains(block)) {
            param_types.emplace_back(type);
        }
    }
    return param_types;
}

std::unordered_set<BlockName> ParamRouter::GetOwnedBlocks(ParameterType ptype) {
    return ParameterBlockRepartition::BLOCKS.at(ptype);
}