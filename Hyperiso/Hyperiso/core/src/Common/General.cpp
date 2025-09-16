#include "General.h"

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

std::ostream &operator<<(std::ostream &os, const LhaID &id) {
    os << id.to_string();
    return os;
}

LhaID::LhaID(const std::string &str_id) {
    for (const auto &num : split(str_id, '_')) {
        parts.emplace_back(std::stol(num));
    }
}

std::string LhaID::to_string() const {
    std::stringstream ss;
    if (!this->parts.empty()) {
        ss << this->parts.at(0);
        for (size_t i = 1; i < this->parts.size(); i++) {
            ss << '_' << this->parts.at(i);
        }
    }
    return ss.str();
}