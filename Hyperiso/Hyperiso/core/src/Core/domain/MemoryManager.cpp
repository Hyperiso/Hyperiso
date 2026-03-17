#include "MemoryManager.h"
#include "Parameters.h"

MemoryManager* MemoryManager::instance = nullptr;

MemoryManager::MemoryManager() : memento(DBMemento()) {
    this->cache.is_ready = false;
}

void MemoryManager::check_if_ready() {
    if (!cache.is_ready) {
        LOG_ERROR("MemoryManager", "Please init the memory manager before using it.");
    }
}

std::shared_ptr<BlockAccessor> MemoryManager::extract_blocks(std::unordered_set<BlockName> block_names) {
    return (*input_cache)[block_names];
}

std::shared_ptr<BlockAccessor> MemoryManager::extract_block_accessor() {
    auto global_ba = std::make_shared<BlockAccessor>();

    for (auto type : cache.parameter_types)
    {
        auto params = Parameters::GetInstance(type);
        if (!params)
            continue;

        auto ba = params->get_block_accessor();
        if (!ba)
            continue;

        for (const auto& [block_name, block_ptr] : *ba)
        {
            global_ba->emplace(block_name, block_ptr);
        }
    }

    if (input_cache)
    {
        for (const auto& [block_name, block_ptr] : *input_cache)
        {
            if (!global_ba->contains(block_name))
            {
                global_ba->emplace(block_name, block_ptr);
            }
        }
    }
    return global_ba;
}

const CorrelationRepository &MemoryManager::get_correlation_repository() {
    check_if_ready();
    return correlation_repository;
}

void MemoryManager::save_input_cache() {
    memento.takeSnapshot(input_cache);
}

void MemoryManager::read_default_input() {
    dl_ba->load(input_cache, paths_provider->default_param_values()); 
    auto obs_blocks = std::make_shared<BlockAccessor>();
    dl_ba->load(obs_blocks, paths_provider->default_obs_values(), true); 
    input_cache = input_cache >> obs_blocks;
    save_input_cache();
    LOG_DEBUG("Default cache stored");

    auto default_param_corr = std::make_shared<CorrelationMatrixPair<ParamId>>();
    dl_cmp_p->load(default_param_corr, paths_provider->default_param_corr().string());
    LOG_DEBUG("Default param correlations loaded");
    auto default_obs_corr = std::make_shared<CorrelationMatrixPair<ExperimentObs>>();
    dl_cmp_o->load(default_obs_corr, paths_provider->default_obs_corr().string());
    LOG_DEBUG("Default observable correlations loaded");
    correlation_repository.set_correlation_matrix(default_param_corr);
    correlation_repository.set_correlation_matrix(default_obs_corr);

    LOG_INFO("Default files loaded");
}

void MemoryManager::read_user_input() {
    // ParamBlockLoader p_loader;
    fs::path ui_paths[3] = {
        paths_provider->user_sm_params(), paths_provider->user_flavor_params(),
        paths_provider->user_decay_params()
    };
    for (auto& path : ui_paths) {
        auto ui_ba = std::make_shared<BlockAccessor>();
        dl_ba->load(ui_ba, path); 
        input_cache = ui_ba >> input_cache;
    }

    auto ui_ba_obs = std::make_shared<BlockAccessor>();
    dl_ba->load(ui_ba_obs, paths_provider->user_obs_values(), true); 
    input_cache = ui_ba_obs >> input_cache;

    save_input_cache();

    auto user_param_corr = std::make_shared<CorrelationMatrixPair<ParamId>>();
    dl_cmp_p->load(user_param_corr, paths_provider->user_param_corr().string());
    auto user_obs_corr = std::make_shared<CorrelationMatrixPair<ExperimentObs>>();
    dl_cmp_o->load(user_obs_corr, paths_provider->user_obs_corr().string());
    correlation_repository.merge_correlation_matrix(user_param_corr);
    correlation_repository.merge_correlation_matrix(user_obs_corr);

    LOG_INFO("User input files loaded");
}

void MemoryManager::read_lha_input(const std::string& lhaFile, const HyperisoConfig& config) {
    fs::path lha_path = this->format_lha_path(lhaFile);
    fs::path spectrum_path = calculate_spectrum(lha_path, config);

    auto lha_ba = std::make_shared<BlockAccessor>();
    dl_ba->load(lha_ba, spectrum_path);
    input_cache = lha_ba >> input_cache;
    save_input_cache();

    LOG_INFO("LHA file loaded");
}

fs::path MemoryManager::calculate_spectrum(fs::path input_lha_path, const HyperisoConfig &config) {
    if (!(config.model == Model::THDM || config.model == Model::SUSY)) {
        return input_lha_path;
    }

    if (!sc) {
        LOG_WARN("MemoryManager", "No ISpectrumCalculator provided, skipping spectrum calculation.");
        return input_lha_path;
    }

    fs::path spectrum_path = paths_provider->spectrum_dir()/input_lha_path.filename();
    // SpectrumCalculator sc;
    sc->calculate_spectrum(input_lha_path, spectrum_path, config.model);
    return spectrum_path;
}

MemoryManager* MemoryManager::GetInstance() {
    if (!MemoryManager::instance) {
        MemoryManager::instance = new MemoryManager();
    }
    return MemoryManager::instance;
}

MemoryManager::MemoryManager(std::shared_ptr<IDataLoader<BlockAccessor>> loader, std::shared_ptr<IDataLoader<CorrelationMatrixPair<ParamId>>> param_corr, std::shared_ptr<IDataLoader<CorrelationMatrixPair<ExperimentObs>>> obs_corr, std::shared_ptr<ISpectrumCalculator> spectrum_c, std::shared_ptr<IPathsProvider> paths_provider_) : memento(DBMemento()) {
    this->sc = spectrum_c;
    
    this->dl_ba = loader;
    this->dl_cmp_p = param_corr;
    this->dl_cmp_o = obs_corr;

    this->paths_provider = paths_provider_;
    this->cache.is_ready = false;
}

MemoryManager* MemoryManager::Create(std::shared_ptr<IDataLoader<BlockAccessor>> loader, std::shared_ptr<IDataLoader<CorrelationMatrixPair<ParamId>>> param_corr, std::shared_ptr<IDataLoader<CorrelationMatrixPair<ExperimentObs>>> obs_corr, std::shared_ptr<ISpectrumCalculator> spectrum_c, std::shared_ptr<IPathsProvider> paths_provider) {
    if (!MemoryManager::instance) {
        MemoryManager::instance = new MemoryManager(loader, param_corr, obs_corr, spectrum_c, paths_provider);
    }
    return MemoryManager::instance;
}

void MemoryManager::init(const std::string& lhaFile, HyperisoConfig config) {
    if (cache.is_ready) {
        LOG_WARN("MemoryManager has already been initialized.");
        return;
    }

    input_cache = std::make_shared<BlockAccessor>();
    this->read_default_input();
    this->read_user_input();
    this->read_lha_input(lhaFile, config);
    cache.lha_path = lhaFile;
    cache.config = config;
    cache.thread_id = std::this_thread::get_id();
    this->deduce_parameter_types(config);
    cache.is_ready = true;

    LOG_DEBUG("Hyperiso successfully initialized !");
}

void MemoryManager::deduce_parameter_types(const HyperisoConfig &config) {
    cache.parameter_types = {ParameterType::SM,
                             ParameterType::FLAVOR,
                             ParameterType::DECAY,
                             ParameterType::OBSERVABLE,
                             ParameterType::PASSTHROUGH,
                             ParameterType::WILSON};
    if (config.model != Model::SM) {
        cache.parameter_types.push_back(ParameterType::BSM);
    }
}

void MemoryManager::switch_lha(const std::string& lhaFile, HyperisoConfig config) {
    ReadyGuard guard(cache.is_ready);
    memento.restore();
    this->read_lha_input(lhaFile, config);
    cache.lha_path = lhaFile; 
    cache.config   = std::move(config);
    this->deduce_parameter_types(cache.config);
    this->cache.flags[InternalFlag::PARAMS_CHANGED] = true;
}

void MemoryManager::reload_user_input(HyperisoConfig config) {
    reload_user_input(cache.lha_path, config);
}

void MemoryManager::reload_user_input(const std::string &lhaFile, HyperisoConfig config) {
    ReadyGuard guard(cache.is_ready);
    memento.restore(2);
    this->read_user_input();
    this->read_lha_input(lhaFile, config);
    
    cache.lha_path = lhaFile;
    cache.config   = std::move(config);
    this->deduce_parameter_types(cache.config);
    this->cache.flags[InternalFlag::PARAMS_CHANGED] = true;

}

fs::path MemoryManager::format_lha_path(const std::string &path) {
    const fs::path input_path(path);
    const fs::path assets_dir = paths_provider->assets_root();
    fs::path full_path;
    if (input_path.is_relative()) { // Path is specified relative to Assets/
        full_path = assets_dir/path;
    } else if (input_path.is_absolute()) {
        full_path = path;
    } else {
        LOG_ERROR("MemoryManager", "LHA File path is undefined:", path);
    }

    if (!std::filesystem::exists(full_path)) {
        LOG_ERROR("MemoryManager", "Cannot find LHA File:", full_path.string());
    }

    LOG_DEBUG("LHA File path:", full_path);
    return full_path;
}

//TODO careful with use_marty, no need for bool but path and model_name 
void MemoryManager::switch_model(Model model) {
    // this->cache.config.flags[ExternalFlag::USE_MARTY] = use_marty;
    this->cache.config.model = model;
    this->deduce_parameter_types(this->cache.config);

    this->cache.flags[InternalFlag::PARAMS_CHANGED] = true;
}

std::unordered_set<ParamId>
MemoryManager::get_all_source_parameters(const std::unordered_set<ParamId>& param_ids) const
{

    auto global_ba = std::make_shared<BlockAccessor>();

    for (auto type : cache.parameter_types)
    {
        auto params = Parameters::GetInstance(type);
        if (!params)
            continue;

        auto ba = params->get_block_accessor();
        if (!ba)
            continue;

        for (const auto& [block_name, block_ptr] : *ba)
        {
            global_ba->emplace(block_name, block_ptr);
        }
    }

    if (input_cache)
    {
        for (const auto& [block_name, block_ptr] : *input_cache)
        {
            if (!global_ba->contains(block_name))
            {
                global_ba->emplace(block_name, block_ptr);
            }
        }
    }

    return global_ba->get_all_source_parameters(param_ids);
}