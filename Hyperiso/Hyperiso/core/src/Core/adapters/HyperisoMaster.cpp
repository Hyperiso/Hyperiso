#include "HyperisoMaster.h"

void HyperisoMaster::init(const std::string &lhaFile, Config config) {
    std::shared_ptr<ParamBlockLoader> pbl = std::make_shared<ParamBlockLoader>();
    std::shared_ptr<CorrelationLoader<ParamId>> cl_param = std::make_shared<CorrelationLoader<ParamId>>();
    std::shared_ptr<CorrelationLoader<ObservableId>> cl_obs = std::make_shared<CorrelationLoader<ObservableId>>();
    MemoryManager::GetInstance()->init(lhaFile, std::move(config), pbl, cl_param, cl_obs);
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

void HyperisoMaster::switch_lha(const std::string &lhaFile, Config config) {
    MemoryManager::GetInstance()->switch_lha(lhaFile, config);
}