#include "MemoryManager.h"

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

const CorrelationRepository &MemoryManager::get_correlation_repository() {
    check_if_ready();
    return correlation_repository;
}

void MemoryManager::save_input_cache() {
    memento.takeSnapshot(input_cache);
}

void MemoryManager::read_default_input(std::shared_ptr<IDataLoader<BlockAccessor>> loader, std::shared_ptr<IDataLoader<CorrelationMatrixPair<ParamId>>> param_corr, std::shared_ptr<IDataLoader<CorrelationMatrixPair<Observables>>> obs_corr) {
    loader->load(input_cache, FilePaths::default_param_values_path); 
    auto obs_blocks = std::make_shared<BlockAccessor>();
    loader->load(obs_blocks, FilePaths::default_obs_values_path); 
    input_cache = input_cache + obs_blocks; 
    LOG_INFO("Default input loaded");
    save_input_cache();
    LOG_INFO("Default cache stored");

    auto default_param_corr = std::make_shared<CorrelationMatrixPair<ParamId>>();
    param_corr->load(default_param_corr, FilePaths::default_param_corr_path.string());
    LOG_INFO("Default param correlations loaded");
    auto default_obs_corr = std::make_shared<CorrelationMatrixPair<Observables>>();
    obs_corr->load(default_obs_corr, FilePaths::default_obs_corr_path.string());
    LOG_INFO("Default observable correlations loaded");
    correlation_repository.set_correlation_matrix(default_param_corr);
    correlation_repository.set_correlation_matrix(default_obs_corr);

    LOG_INFO("Default files loaded");
}

void MemoryManager::read_user_input(std::shared_ptr<IDataLoader<BlockAccessor>> loader, std::shared_ptr<IDataLoader<CorrelationMatrixPair<ParamId>>> param_corr, std::shared_ptr<IDataLoader<CorrelationMatrixPair<Observables>>> obs_corr) {
    // ParamBlockLoader p_loader;
    fs::path ui_paths[4] = {FilePaths::user_sm_params_path, FilePaths::user_flavor_params_path, FilePaths::user_decay_params_path, FilePaths::user_obs_values_path};
    for (auto& path : ui_paths) {
        auto ui_ba = std::make_shared<BlockAccessor>();
        loader->load(ui_ba, path); 
        input_cache = ui_ba >> input_cache;
    }
    save_input_cache();

    auto user_param_corr = std::make_shared<CorrelationMatrixPair<ParamId>>();
    param_corr->load(user_param_corr, FilePaths::user_param_corr_path.string());
    auto user_obs_corr = std::make_shared<CorrelationMatrixPair<Observables>>();
    obs_corr->load(user_obs_corr, FilePaths::user_obs_corr_path.string());
    correlation_repository.merge_correlation_matrix(user_param_corr);
    correlation_repository.merge_correlation_matrix(user_obs_corr);

    LOG_DEBUG("User input files loaded");
}

void MemoryManager::read_lha_input(const std::string& lhaFile, const Config& config) {
    fs::path lha_path = this->format_lha_path(lhaFile);
    fs::path spectrum_path = calculate_spectrum(lha_path, config);

    ParamBlockLoader p_loader;
    auto lha_ba = std::make_shared<BlockAccessor>();
    p_loader.load(lha_ba, spectrum_path);
    input_cache = lha_ba >> input_cache;
    save_input_cache();

    LOG_DEBUG("LHA file loaded");
}

fs::path MemoryManager::calculate_spectrum(fs::path input_lha_path, const Config &config) {
    if (config.flags.at(ExternalFlag::USE_MARTY) || !(config.model == Model::THDM || config.model == Model::SUSY)) {
        return input_lha_path;
    }

    fs::path spectrum_path = DirPaths::spectrum_dir_path/input_lha_path.filename();
    SpectrumCalculator sc;
    sc.calculate_spectrum(input_lha_path, spectrum_path, config.model);
    return spectrum_path;
}

MemoryManager* MemoryManager::GetInstance() {
    if (!MemoryManager::instance) {
        MemoryManager::instance = new MemoryManager();
    }
    return MemoryManager::instance;
}

void MemoryManager::init(const std::string& lhaFile, Config config, std::shared_ptr<IDataLoader<BlockAccessor>> loader, std::shared_ptr<IDataLoader<CorrelationMatrixPair<ParamId>>> param_corr, std::shared_ptr<IDataLoader<CorrelationMatrixPair<Observables>>> obs_corr) {
    if (cache.is_ready) {
        LOG_WARN("MemoryManager has already been initialized.");
        return;
    }

    input_cache = std::make_shared<BlockAccessor>();
    this->read_default_input(loader, param_corr, obs_corr);
    this->read_user_input(loader, param_corr, obs_corr);
    this->read_lha_input(lhaFile, config);

    cache.lha_path = lhaFile;
    cache.config = config;
    cache.thread_id = std::this_thread::get_id();
    this->deduce_parameter_types(config);
    cache.is_ready = true;

    LOG_DEBUG("Hyperiso successfully initialized !");
}

void MemoryManager::deduce_parameter_types(const Config &config) {
    cache.parameter_types = {ParameterType::SM,
                             ParameterType::FLAVOR,
                             ParameterType::DECAY,
                             ParameterType::OBSERVABLE,
                             ParameterType::PASSTHROUGH,
                             ParameterType::WILSON};
    if (config.model != Model::SM)
        cache.parameter_types.push_back(ParameterType::BSM);
    if (config.flags.at(ExternalFlag::HAS_WILSON_INPUT))
        cache.parameter_types.push_back(ParameterType::WILSON);
}

void MemoryManager::switch_lha(const std::string& lhaFile, Config config) {
    this->cache.is_ready = false;
    memento.restore();
    this->read_lha_input(lhaFile, config);
    this->cache.flags.at(InternalFlag::PARAMS_CHANGED) = true;
}

void MemoryManager::reload_user_input(Config config, std::shared_ptr<IDataLoader<BlockAccessor>> loader, std::shared_ptr<IDataLoader<CorrelationMatrixPair<ParamId>>> param_corr, std::shared_ptr<IDataLoader<CorrelationMatrixPair<Observables>>> obs_corr) {
    reload_user_input(cache.lha_path, config, loader, param_corr, obs_corr);
}

void MemoryManager::reload_user_input(const std::string &lhaFile, Config config, std::shared_ptr<IDataLoader<BlockAccessor>> loader, std::shared_ptr<IDataLoader<CorrelationMatrixPair<ParamId>>> param_corr, std::shared_ptr<IDataLoader<CorrelationMatrixPair<Observables>>> obs_corr) {
    this->cache.is_ready = false;
    memento.restore(2);
    this->read_user_input(loader, param_corr, obs_corr);
    this->read_lha_input(lhaFile, config);
    this->cache.flags.at(InternalFlag::PARAMS_CHANGED) = true;
}

fs::path MemoryManager::format_lha_path(const std::string &path) {
    const fs::path input_path(path);
    const fs::path assets_dir(project_assets_root.data());
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

void MemoryManager::switch_model(Model model, bool use_marty) {
    this->cache.config.flags.at(ExternalFlag::USE_MARTY) = use_marty;
    this->cache.parameter_types.erase(std::find(this->cache.parameter_types.begin(), this->cache.parameter_types.end(), static_cast<ParameterType>(static_cast<int>(cache.config.model))));
    this->cache.parameter_types.push_back(static_cast<ParameterType>(static_cast<int>(model)));
    this->cache.config.model = model;
    this->cache.flags.at(InternalFlag::PARAMS_CHANGED) = true;
}