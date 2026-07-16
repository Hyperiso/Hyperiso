#include "LhaParamsHelper.h"

const std::map<BlockName, std::vector<std::vector<long>>> LhaParamsHelper::minimal_blocks = {
    {"FMASS", {{211}, {321}, {323}, {411}, {421}, {423}, {431}, {511}, {521}, {531}}},
    {"FLIFE", {{211}, {321}, {323}, {411}, {421}, {431}, {511}, {521}, {531}}},
    {"FCONST", {{511, 1}, {521, 1}, {531, 1}, {323, 1}, {323, 2}}},
};

std::vector<std::vector<long>> LhaParamsHelper::get_minimal_content(const BlockName &block_name) {
    if (LhaParamsHelper::minimal_blocks.contains(block_name)) {
        return LhaParamsHelper::minimal_blocks.at(block_name);
    }
    LOG_ERROR("LhaParamsHelper", "Unknown block", block_name);
}