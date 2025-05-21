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