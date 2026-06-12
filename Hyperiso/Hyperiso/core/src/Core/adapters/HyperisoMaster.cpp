#include "HyperisoMaster.h"

#include <utility>

void HyperisoMaster::ensure_memory_manager_created() {
    std::shared_ptr<ParamBlockLoader> pbl = std::make_shared<ParamBlockLoader>();
    std::shared_ptr<CorrelationLoader<ParamId>> cl_param = std::make_shared<CorrelationLoader<ParamId>>();
    std::shared_ptr<CorrelationLoader<ExperimentObs>> cl_obs = std::make_shared<CorrelationLoader<ExperimentObs>>();
    std::shared_ptr<SpectrumCalculator> spectrum_c = std::make_shared<SpectrumCalculator>();
    std::shared_ptr<DefaultPathsProvider> dpp = std::make_shared<DefaultPathsProvider>();
    MemoryManager::Create(pbl, cl_param, cl_obs, spectrum_c, dpp, pbl);
}

bool HyperisoMaster::should_validate_marty_runtime(const HyperisoConfig& config) const {
    const auto flag_it = config.flags.find(ExternalFlag::HYP_AS_SM_MARTY);
    const bool hyp_as_sm_marty = flag_it != config.flags.end() && flag_it->second;
    return hyp_as_sm_marty || config.model == Model::MARTY;
}

bool HyperisoMaster::validate_marty_runtime_if_needed(const HyperisoConfig& config, const std::string& context) const {
    if (!should_validate_marty_runtime(config)) {
        return true;
    }

    const auto install = MartyRuntimeConfig::require_available(context);
    return install.valid;
}

void HyperisoMaster::init(const std::string &lhaFile, HyperisoConfig config) {
    ensure_memory_manager_created();
    if (!validate_marty_runtime_if_needed(config, "HyperisoMaster::init")) {
        return;
    }
    MemoryManager::GetInstance()->init(lhaFile, std::move(config));
}


void HyperisoMaster::pre_init_add_block(BlockName blockName,
                              size_t itemCount,
                              size_t valueIdx,
                              int scaleIdx,
                              int rgIdx,
                              int binIdx,
                              bool globalScale)
{
    ensure_memory_manager_created();

    auto* mm = MemoryManager::GetInstance();
    if (mm->is_ready()) {
        LOG_WARN("HyperisoMaster", "pre_init_add_block must be called before init; this prototype will only affect future LHA reloads.");
    }

    mm->add_lha_prototype(blockName, itemCount, valueIdx, scaleIdx, rgIdx, binIdx, globalScale);
}

void HyperisoMaster::pre_init_add_blocks(const std::vector<LhaPrototypeSpec>& prototypes)
{
    ensure_memory_manager_created();

    auto* mm = MemoryManager::GetInstance();
    if (mm->is_ready()) {
        LOG_WARN("HyperisoMaster", "pre_init_add_blocks must be called before init; these prototypes will only affect future LHA reloads.");
    }

    mm->add_lha_prototypes(prototypes);
}

void HyperisoMaster::pre_init_set_marty_path(const std::string& martyInstallPath)
{
    const auto install = MartyRuntimeConfig::set_external_install_path(martyInstallPath);
    if (!install.valid) {
        LOG_ERROR("MartyConfigError", install.error);
        return;
    }

    if (install.has_executable) {
        LOG_DEBUG("MARTY executable detected at", install.marty_executable.string());
    } else {
        LOG_WARN("MartyRuntimeConfig", "No bin/marty executable found under", install.prefix.string(), "continuing because libmarty and marty.h are present.");
    }

    LOG_INFO("MARTY runtime path registered:", install.prefix.string());
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
    ensure_memory_manager_created();
    if (!validate_marty_runtime_if_needed(config, "HyperisoMaster::switch_lha")) {
        return;
    }
    MemoryManager::GetInstance()->switch_lha(lhaFile, std::move(config));
}
