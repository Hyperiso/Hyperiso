#include "MemoryManager.h"

namespace fs = std::filesystem;

MemoryManager* MemoryManager::instance = nullptr;

MemoryManager::MemoryManager() {
    this->cache.is_ready = false;
}

void MemoryManager::check_if_ready() {
    if (!cache.is_ready) {
        LOG_ERROR("MemoryManager", "Please init the memory manager before using it.");
    }
}

std::shared_ptr<BlockAccessor> MemoryManager::get_blocks(std::unordered_set<std::string> block_names) {
    return (*input_cache)[block_names];
}

void MemoryManager::save_input_cache() {
    DBMemento().takeSnapshot(input_cache);
}

std::shared_ptr<BlockAccessor> MemoryManager::read_input_files(fs::path lha_path) {
    /* Default input */

    DBManager manager;
    auto default_param_values_root = manager.read_from_file(FilePaths::default_param_values_path);
    auto default_obs_values_root = manager.read_from_file(FilePaths::default_obs_values_path);

    // TODO : insert file check here
    auto input_blocks = ParamBlockAdapter::from_db_node(default_param_values_root); 
    input_blocks = input_blocks + ParamBlockAdapter::from_db_node(default_obs_values_root); 
    DBMemento memento;
    memento.takeSnapshot(input_blocks);

    // auto default_param_corr_root = json_parser->readFromFile(FilePaths::default_param_corr_path.string());
    // auto default_obs_corr_root = json_parser->readFromFile(FilePaths::default_obs_corr_path.string());
    // CorrelationMatrixPair default_param_corr = CorrelationAdapter::from_db_node<ParamId>(default_param_corr_root);
    // CorrelationMatrixPair default_obs_corr = CorrelationAdapter::from_db_node<Observables>(default_obs_corr_root);
    // CorrelationRepository cr;
    // cr.set_correlation_matrix(std::move(default_param_corr));
    // cr.set_correlation_matrix(std::move(default_obs_corr));

    /* User input */

    fs::path ui_paths[4] = {FilePaths::user_sm_params_path, FilePaths::user_flavor_params_path, FilePaths::user_decay_params_path, FilePaths::user_obs_values_path};
    for (auto& path : ui_paths) {
        auto ui_root = manager.read_from_file(path);
        auto ui_ba = ParamBlockAdapter::from_db_node(ui_root); 
        input_blocks = ui_ba >> input_blocks;
    }
    memento.takeSnapshot(input_blocks);

    // auto user_param_corr_root = yaml_parser->readFromFile(FilePaths::user_param_corr_path.string());
    // auto user_obs_corr_root = yaml_parser->readFromFile(FilePaths::user_obs_corr_path.string());
    // CorrelationMatrixPair user_param_corr = CorrelationAdapter::from_db_node<ParamId>(user_param_corr_root);
    // CorrelationMatrixPair user_obs_corr = CorrelationAdapter::from_db_node<Observables>(user_obs_corr_root);
    // cr.merge_correlation_matrix(std::move(user_param_corr));
    // cr.merge_correlation_matrix(std::move(user_obs_corr));

    /* LHA input */
    auto lha_root = manager.read_from_file(lha_path);

    auto lha_param_ba = ParamBlockAdapter::from_db_node(lha_root);
    input_blocks = lha_param_ba >> input_blocks;
    memento.takeSnapshot(input_blocks);

    return input_blocks;
}

MemoryManager* MemoryManager::GetInstance() {
    if (!MemoryManager::instance) {
        MemoryManager::instance = new MemoryManager();
    }
    return MemoryManager::instance;
}

void MemoryManager::init(const std::string& lhaFile, const Config& config) {
    if (cache.is_ready) {
        LOG_WARN("MemoryManager has already been initialized.");
        return;
    }
    fs::path lha_path = format_lha_path(lhaFile);
    fs::path spectrum_path = lha_path;
    if (!config.use_marty && (config.model == Model::THDM || config.model == Model::SUSY) && !config.is_spectrum) {
        spectrum_path = DirPaths::spectrum_dir_path/lha_path.filename();
        LOG_DEBUG("Starting spectrum calculation...");
        CalculatorType calculatorType = config.model == Model::THDM ? CalculatorType::TwoHDM : CalculatorType::Softsusy;
        GeneralCalculatorFactory::executeCommand(calculatorType, "calculateSpectrum", lha_path, spectrum_path.string());
        LOG_DEBUG("Spectrum calculation ran sucessfully");
    }
    auto block_accessor = read_input_files(spectrum_path);
    input_cache = block_accessor;

    cache.lha_path = lha_path;
    cache.config = config;
    cache.thread_id = std::this_thread::get_id();
    deduce_parameter_types(config);
    cache.is_ready = true;

    for (auto &&m : cache.parameter_types) {
        LOG_DEBUG("Initializing parameters ", (int)m);
        Parameters::GetInstance(m);
    }
}

void MemoryManager::deduce_parameter_types(const Config &config) {
    cache.parameter_types = {ParameterType::SM,
                             ParameterType::FLAVOR,
                             ParameterType::DECAY};
    if (config.model != Model::SM)
        cache.parameter_types.push_back(static_cast<ParameterType>(static_cast<int>(config.model)));
    if (config.has_wilsons)
        cache.parameter_types.push_back(ParameterType::WILSON);
}

void MemoryManager::switch_lha(const std::string& lhaFile, Config config) {
    for (auto& param : cache.parameter_types) {
        Parameters::GetInstance(param)->CleanupInstance(param);
    }
    this->cache.is_ready = false;
    this->init(lhaFile, config);
    this->cache.param_cache_okay = false;
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
    this->cache.config.model = model;
    this->cache.config.use_marty = use_marty;
    this->cache.parameter_types.push_back(static_cast<ParameterType>(static_cast<int>(model)));
    Parameters::GetInstance(static_cast<ParameterType>(static_cast<int>(model)));
    this->cache.param_cache_okay = false;
}

std::map<LhaID, double> MemoryManager::get_block_infos(const std::string& block, ParameterType param_type) {
    return Parameters::GetInstance(param_type)->get_block_infos(block);
}

std::unordered_set<std::string> MemoryManager::get_blocks_list(ParameterType param_type) {
    return Parameters::GetInstance(param_type)->get_blocks_list();
}

std::unordered_set<std::string> MemoryManager::get_all_blocks() {
    return this->input_cache->get_block_names();
}

std::vector<ParameterType> MemoryManager::get_type_of_block(const std::string& block) {
    std::vector<ParameterType> param_type;
    for (auto& elem : cache.parameter_types) {
        for (auto& block_ : get_blocks_list(elem)) {
            if (block == block_) {
                param_type.push_back(elem);
                break;
            }
        }
    }
    return param_type;
}