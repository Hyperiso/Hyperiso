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

template<typename T>
std::unique_ptr<T, void(*)(void*)> MemoryManager::makeUniquePtr(T* ptr) {
    return std::unique_ptr<T, void(*)(void*)>(ptr, [](void* p) { std::free(p); });
}

template<typename T>
std::unique_ptr<T, void(*)(void*)> MemoryManager::allocate() {
    T* ptr = static_cast<T*>(std::malloc(sizeof(T)));
    if (!ptr) {
        throw std::bad_alloc();
    }
    return makeUniquePtr(ptr);
}

/**
 * @brief Search for the nearest directory containing "hyperiso" in its name.
 * 
 * This function searches for the nearest directory containing "hyperiso" in its name
 * starting from the current directory and moving upwards in the directory tree.
 * It returns the absolute path of the found directory if one is found, otherwise an empty string.
 * 
 * @return The absolute path of the nearest directory containing "hyperiso" in its name.
 */
std::string MemoryManager::findNearestHyperisoDirectory() {
    fs::path currentDir = fs::current_path();

    std::cout << "Looking for project root..." << std::endl;

    // Iterate through parent directories until "hyperiso" is found or root directory is reached
    while (!currentDir.empty()) {
        for (const auto& entry : fs::directory_iterator(currentDir)) {
            if (entry.is_directory() && entry.path().filename().string().find("Hyperiso") != std::string::npos) {
                LOG_INFO("Project root folder is " +entry.path().string());
                return entry.path().string() + "/";
            }
        }
        currentDir = currentDir.parent_path(); // Move to the parent directory
    }

    // If "hyperiso" directory is not found in any parent directory
    LOG_ERROR("OldError", "Error: Nearest directory containing 'hyperiso' not found.");
    return "";
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

void MemoryManager::init(const std::string& lhaFile, Model model, bool is_spectrum, bool has_wilsons, bool has_obs) {
    if (cache.is_ready) {
        LOG_WARN("MemoryManager has already been initialized.");
        return;
    }

    std::stringstream ss;
    ss << project_root.data() << "/" << lhaFile;
    cache.reader = std::make_unique<LhaReader>(LhaReader(ss.str()));
    cache.reader->readAll();
    cache.lha_path = std::filesystem::u8path(ss.str());
    cache.obs_cov_path = std::filesystem::u8path(project_root.data() + std::string("/DataBase/Exp/observable_covariance.json"));
    cache.param_cov_path = std::filesystem::u8path(project_root.data() + std::string("/DataBase/Exp/parameter_covariance.json"));
    cache.model = model;
    cache.is_spectrum = is_spectrum;
    cache.has_wilsons = has_wilsons;
    cache.has_obs = has_obs;
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

void MemoryManager::set_observable_covariance_input_file(const std::string &path) {
    cache.obs_cov_path = path;
}

void MemoryManager::set_parameter_covariance_input_file(const std::string &path) {
    cache.param_cov_path = path;
}
