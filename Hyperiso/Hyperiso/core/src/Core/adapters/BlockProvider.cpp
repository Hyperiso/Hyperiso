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