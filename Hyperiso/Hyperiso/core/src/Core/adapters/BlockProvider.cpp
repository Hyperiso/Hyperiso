#include "BlockProvider.h"

bool BlockProvider::exists(const std::string& blockname, ParameterType pt) {
    auto blocks_list = Parameters::GetInstance(pt)->get_blocks_list();
    for (auto& elem: blocks_list) {
        if (elem == blockname) {
            return true;
        }
    }
    return false;
}

void BlockProvider::log_all_blocks(ParameterType type) {
    LOG_INFO(Parameters::GetInstance(type));
}

void BlockProvider::log_block(ParameterType type, const std::string& blockname) {
    Parameters::GetInstance(type)->print_block(blockname);
}

std::map<LhaID, scalar_t> BlockProvider::get_block(ParameterType type, const std::string& blockname) {
    return Parameters::GetInstance(type)->get_block_infos(blockname);
}

std::unordered_set<std::string> get_all_blocks(ParameterType type) {
    auto blocks_list = Parameters::GetInstance(type)->get_blocks_list();
    std::unordered_set<std::string> out;
    for (auto& elem: blocks_list) {
        out.insert(elem.to_string());
    }
    return out;
}