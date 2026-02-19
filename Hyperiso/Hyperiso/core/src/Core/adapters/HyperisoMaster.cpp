#include "HyperisoMaster.h"

void HyperisoMaster::init(const std::string &lhaFile, HyperisoConfig config) {
    std::shared_ptr<ParamBlockLoader> pbl = std::make_shared<ParamBlockLoader>();
    std::shared_ptr<CorrelationLoader<ParamId>> cl_param = std::make_shared<CorrelationLoader<ParamId>>();
    std::shared_ptr<CorrelationLoader<BinnedObservableId>> cl_obs = std::make_shared<CorrelationLoader<BinnedObservableId>>();
    std::shared_ptr<SpectrumCalculator> spectrum_c = std::make_shared<SpectrumCalculator>();
    std::shared_ptr<DefaultPathsProvider> dpp = std::make_shared<DefaultPathsProvider>();
    MemoryManager::Create(pbl, cl_param, cl_obs, spectrum_c, dpp);
    MemoryManager::GetInstance()->init(lhaFile, std::move(config));
}

void HyperisoMaster::init(const std::string &lhaFile) {
    init(lhaFile, HyperisoConfig());
}

bool HyperisoMaster::check_flag(ExternalFlag flag) {
    return MemoryManager::GetInstance()->getMemoryCache().config.flags.at(flag);
}

Model HyperisoMaster::get_model() {
    return MemoryManager::GetInstance()->getMemoryCache().config.model;
}

void HyperisoMaster::switch_lha(const std::string &lhaFile, HyperisoConfig config) {
    std::shared_ptr<ParamBlockLoader> pbl = std::make_shared<ParamBlockLoader>();
    MemoryManager::GetInstance()->switch_lha(lhaFile, config);
}