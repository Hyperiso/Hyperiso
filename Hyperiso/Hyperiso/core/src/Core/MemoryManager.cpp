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

std::shared_ptr<BlockAccessor> MemoryManager::get_blocks(std::vector<std::string> block_names) {
    return (*input_cache)[block_names];
}

std::shared_ptr<BlockAccessor> MemoryManager::read_params(fs::path lha_path) {
    /* Default input */

    auto json_parser = ParserFactory::createParser(ParserFactory::Type::JSON);
    auto default_param_values_root = json_parser->readFromFile(FilePaths::default_param_values_path.string());
    // TODO : insert file check here
    auto param_values_ba = BlocksCreator::from_db_node(default_param_values_root); 
    DBMemento memento;
    memento.takeSnapshot(param_values_ba);

    // TODO : read default correlations between params
    // TODO : read default observable values and correlations

    /* User input */

    auto yaml_parser = ParserFactory::createParser(ParserFactory::Type::YAML);
    fs::path ui_paths[3] = {FilePaths::user_sm_params_path, FilePaths::user_flavor_params_path, FilePaths::user_decay_params_path};
    for (auto& path : ui_paths) {
        auto ui_root = yaml_parser->readFromFile(path.string());
        auto ui_ba = BlocksCreator::from_db_node(ui_root); 
        param_values_ba = ui_ba >> param_values_ba;
    }
    memento.takeSnapshot(param_values_ba);

    // TODO : read user correlations between params
    // TODO : read user observable values and correlations

    /* LHA input */
    auto lha_reader = std::make_shared<LhaReader>(LhaReader(lha_path.string()));
    lha_reader->readAll();

    auto lha_param_ba = BlocksCreator::from_lha_reader(lha_reader);
    param_values_ba = lha_param_ba >> param_values_ba;
    memento.takeSnapshot(param_values_ba);

    return param_values_ba;
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

    auto block_accessor = read_params(spectrum_path);

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

std::vector<std::string> MemoryManager::get_blocks_list(ParameterType param_type) {
    return Parameters::GetInstance(param_type)->get_blocks_list();
}

std::vector<std::string> MemoryManager::get_all_blocks() {
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