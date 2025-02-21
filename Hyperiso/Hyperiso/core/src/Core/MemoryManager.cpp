#include "MemoryManager.h"
#include "Parameters.h"
#include <filesystem>

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


MemoryManager* MemoryManager::GetInstance() {
    if (!MemoryManager::instance) {
        MemoryManager::instance = new MemoryManager();
    }
    return MemoryManager::instance;
}

LhaReader* MemoryManager::getReader() {
    check_if_ready();
    return cache.reader.get();
}

void MemoryManager::init(const std::string& lhaFile, Model model, bool use_marty, bool is_spectrum, bool has_wilsons, bool has_obs) {
    if (cache.is_ready) {
        LOG_WARN("MemoryManager has already been initialized.");
        return;
    }
    const std::filesystem::path path(lhaFile);
    const std::filesystem::path dir_path(project_assets_root.data());
    std::filesystem::path full_path;
    if (path.is_relative()) {
        full_path = dir_path/path;
    } else if (path.is_absolute()) {
        full_path = path;
    } else {
        LOG_ERROR("PathError", "File not relative or absolute");
    }
    if (!std::filesystem::exists(full_path)) {
        LOG_ERROR("PathError", "Invalid lha path :", full_path.string());
    }
    LOG_INFO("lha path", full_path);
    std::stringstream ss;
    ss << full_path.string();
    cache.reader = std::make_shared<LhaReader>(LhaReader(ss.str()));
    cache.reader->readAll();
    cache.lha_path = std::filesystem::u8path(ss.str());
    cache.obs_cov_path = std::filesystem::u8path(project_assets_root.data() + std::string("default/observables_exp.json"));
    cache.param_cov_path = std::filesystem::u8path(project_assets_root.data() + std::string("default/parameters_exp.json"));
    cache.model = model;
    cache.is_spectrum = is_spectrum;
    cache.has_wilsons = has_wilsons;
    cache.has_obs = has_obs;
    cache.use_marty = use_marty;
    cache.thread_id = std::this_thread::get_id();
    
    cache.parameter_types = {ParameterType::SM, ParameterType::FLAVOR, ParameterType::FF};
    if (model != Model::SM)
        cache.parameter_types.push_back(static_cast<ParameterType>(static_cast<int>(model)));
    if (has_wilsons)
        cache.parameter_types.push_back(ParameterType::WILSON);

    cache.is_ready = true;

    for (auto &&m : cache.parameter_types) {
        LOG_DEBUG("Initializing parameters ", (int)m);
        Parameters::GetInstance(m);
    }
}

void MemoryManager::switch_lha(const std::string& lhaFile, Model model, bool use_marty, bool is_spectrum, bool has_wilsons, bool has_obs) {
    for (auto& param : cache.parameter_types) {
        Parameters::GetInstance(param)->CleanupInstance(param);
    }
    this->cache.is_ready = false;
    this->init(lhaFile, model, use_marty, is_spectrum, has_wilsons, has_obs);
    this->cache.param_cache_okay = false;

}

void MemoryManager::switch_model(Model model, bool use_marty) {
    this->cache.model = model;
    this->cache.use_marty = use_marty;
    this->cache.parameter_types.push_back(static_cast<ParameterType>(static_cast<int>(model)));
    Parameters::GetInstance(static_cast<ParameterType>(static_cast<int>(model)));
    this->cache.param_cache_okay = false;
}

void MemoryManager::set_observable_covariance_input_file(const std::string &path) {
    cache.obs_cov_path = path;
}

void MemoryManager::set_parameter_covariance_input_file(const std::string &path) {
    cache.param_cov_path = path;
}

std::map<int, double> MemoryManager::get_block_infos(const std::string& block, ParameterType param_type) {
    return Parameters::GetInstance(param_type)->get_block_infos(block);
}

std::vector<std::string> MemoryManager::get_blocks_list(ParameterType param_type) {
        return Parameters::GetInstance(param_type)->get_blocks_list();
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