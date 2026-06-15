#include "HyperisoMaster.h"

#include <algorithm>
#include <cctype>
#include <memory>
#include <system_error>
#include <utility>

namespace {

std::string api_path_name(APIPath path_name)
{
    switch (path_name) {
    case APIPath::LHA_PATH:             return "LHA_PATH";
    case APIPath::ASSETS_ROOT:          return "ASSETS_ROOT";
    case APIPath::DEFAULT_PARAM_VALUES: return "DEFAULT_PARAM_VALUES";
    case APIPath::DEFAULT_OBS_VALUES:   return "DEFAULT_OBS_VALUES";
    case APIPath::DEFAULT_PARAM_CORR:   return "DEFAULT_PARAM_CORR";
    case APIPath::DEFAULT_OBS_CORR:     return "DEFAULT_OBS_CORR";
    case APIPath::USER_SM_PARAMS:       return "USER_SM_PARAMS";
    case APIPath::USER_FLAVOR_PARAMS:   return "USER_FLAVOR_PARAMS";
    case APIPath::USER_DECAY_PARAMS:    return "USER_DECAY_PARAMS";
    case APIPath::USER_OBS_VALUES:      return "USER_OBS_VALUES";
    case APIPath::USER_PARAM_CORR:      return "USER_PARAM_CORR";
    case APIPath::USER_OBS_CORR:        return "USER_OBS_CORR";
    case APIPath::PARAM_MAPPING_DIR:    return "PARAM_MAPPING_DIR";
    case APIPath::TEMPLATE_DIR:         return "TEMPLATE_DIR";
    case APIPath::SPECTRUM_DIR:         return "SPECTRUM_DIR";
    }
    return "UNKNOWN_API_PATH";
}

bool is_default_input_path(APIPath path_name)
{
    switch (path_name) {
    case APIPath::DEFAULT_PARAM_VALUES:
    case APIPath::DEFAULT_OBS_VALUES:
    case APIPath::DEFAULT_PARAM_CORR:
    case APIPath::DEFAULT_OBS_CORR:
        return true;
    default:
        return false;
    }
}

bool is_user_input_path(APIPath path_name)
{
    switch (path_name) {
    case APIPath::USER_SM_PARAMS:
    case APIPath::USER_FLAVOR_PARAMS:
    case APIPath::USER_DECAY_PARAMS:
    case APIPath::USER_OBS_VALUES:
    case APIPath::USER_PARAM_CORR:
    case APIPath::USER_OBS_CORR:
        return true;
    default:
        return false;
    }
}

bool is_directory_path(APIPath path_name)
{
    switch (path_name) {
    case APIPath::ASSETS_ROOT:
    case APIPath::PARAM_MAPPING_DIR:
    case APIPath::TEMPLATE_DIR:
    case APIPath::SPECTRUM_DIR:
        return true;
    default:
        return false;
    }
}

std::string lower_extension(const fs::path& path)
{
    std::string extension = path.extension().string();
    std::transform(extension.begin(), extension.end(), extension.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return extension;
}

bool validate_pre_init_path_override(APIPath path_name, const fs::path& path)
{
    if (path_name == APIPath::LHA_PATH) {
        LOG_ERROR("HyperisoMaster", "LHA_PATH is provided through init()/switch_lha() and cannot be overridden through pre_init_set_paths.");
        return false;
    }

    std::error_code ec;
    if (!fs::exists(path, ec) || ec) {
        LOG_ERROR("HyperisoMaster", "Path override for", api_path_name(path_name), "does not exist:", path.string());
        return false;
    }

    if (is_directory_path(path_name)) {
        if (!fs::is_directory(path, ec) || ec) {
            LOG_ERROR("HyperisoMaster", "Path override for", api_path_name(path_name), "must be an existing directory:", path.string());
            return false;
        }
        return true;
    }

    if (!fs::is_regular_file(path, ec) || ec) {
        LOG_ERROR("HyperisoMaster", "Path override for", api_path_name(path_name), "must be an existing file:", path.string());
        return false;
    }

    const std::string extension = lower_extension(path);
    if (is_default_input_path(path_name)) {
        if (extension != ".json") {
            LOG_ERROR("HyperisoMaster", "Path override for", api_path_name(path_name), "must use the .json extension:", path.string());
            return false;
        }
        LOG_WARN("HyperisoMaster", "Overriding a Hyperiso default path is not the expected runtime behavior. Make sure you know what you are doing:", api_path_name(path_name));
        return true;
    }

    if (is_user_input_path(path_name)) {
        if (extension != ".yaml" && extension != ".yml") {
            LOG_ERROR("HyperisoMaster", "Path override for", api_path_name(path_name), "must use the .yaml or .yml extension:", path.string());
            return false;
        }
        return true;
    }

    LOG_ERROR("HyperisoMaster", "Unsupported APIPath override:", api_path_name(path_name));
    return false;
}

} // namespace

void HyperisoMaster::ensure_memory_manager_created() {
    std::shared_ptr<ParamBlockLoader> pbl = std::make_shared<ParamBlockLoader>();
    std::shared_ptr<CorrelationLoader<ParamId>> cl_param = std::make_shared<CorrelationLoader<ParamId>>();
    std::shared_ptr<CorrelationLoader<ExperimentObs>> cl_obs = std::make_shared<CorrelationLoader<ExperimentObs>>();
    std::shared_ptr<SpectrumCalculator> spectrum_c = std::make_shared<SpectrumCalculator>();
    std::shared_ptr<IPathsProvider> dpp = path_overrides.empty()
        ? std::static_pointer_cast<IPathsProvider>(std::make_shared<DefaultPathsProvider>())
        : std::static_pointer_cast<IPathsProvider>(std::make_shared<OverridePathsProvider>(path_overrides));
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

void HyperisoMaster::pre_init_set_paths(const std::map<APIPath, std::string>& pathOverrides)
{
    ensure_memory_manager_created();

    auto* mm = MemoryManager::GetInstance();
    if (mm->is_ready()) {
        LOG_WARN("HyperisoMaster", "pre_init_set_paths must be called before init; these path overrides will only affect future reloads.");
    }

    if (pathOverrides.empty()) {
        LOG_WARN("HyperisoMaster", "pre_init_set_paths called with an empty override map; keeping existing paths.");
        return;
    }

    bool installed_any_path = false;
    for (const auto& [path_name, raw_path] : pathOverrides) {
        const fs::path candidate(raw_path);
        if (!validate_pre_init_path_override(path_name, candidate)) {
            continue;
        }

        std::error_code ec;
        const fs::path normalized = fs::weakly_canonical(candidate, ec);
        path_overrides[path_name] = ec ? fs::absolute(candidate) : normalized;
        installed_any_path = true;
    }

    if (!installed_any_path) {
        return;
    }

    mm->set_paths_provider(std::make_shared<OverridePathsProvider>(path_overrides));
    LOG_INFO("Hyperiso path overrides registered:", path_overrides.size());
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
