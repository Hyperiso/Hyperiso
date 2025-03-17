#include "HyperisoMaster.h"

void HyperisoMaster::init(const std::string &lhaFile, Config config) {
    MemoryManager::GetInstance()->init(lhaFile, std::move(config));
}

void HyperisoMaster::init(const std::string &lhaFile) {
    init(lhaFile, Config());
}

bool HyperisoMaster::check_flag(ExternalFlag flag) {
    return MemoryManager::GetInstance()->getMemoryCache().config.flags.at(flag);
}

Model HyperisoMaster::get_model() {
    return MemoryManager::GetInstance()->getMemoryCache().config.model;
}
