#include "MemoryManager.h"
#include "Parameters.h"
#include <filesystem>

namespace fs = std::filesystem;

MemoryManager* MemoryManager::instance = nullptr;

// Creation pointeur unique
template<typename T>
std::unique_ptr<T, void(*)(void*)> MemoryManager::makeUniquePtr(T* ptr) {
    return std::unique_ptr<T, void(*)(void*)>(ptr, [](void* p) { std::free(p); });
}

// Méthode pour allouer de la mémoire
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
    Logger* logger = Logger::getInstance();
    fs::path currentDir = fs::current_path();

    // Iterate through parent directories until "hyperiso" is found or root directory is reached
    while (!currentDir.empty()) {
        for (const auto& entry : fs::directory_iterator(currentDir)) {
            if (entry.is_directory() && entry.path().filename().string().find("Hyperiso") != std::string::npos) {
                logger->info("Project root folder is " +entry.path().string());
                return entry.path().string() + "/";
            }
        }
        currentDir = currentDir.parent_path(); // Move to the parent directory
    }

    // If "hyperiso" directory is not found in any parent directory
    logger->error("Error: Nearest directory containing 'hyperiso' not found.");
    return "";
}

MemoryManager* MemoryManager::GetInstance(std::string lhaFile, std::vector<int> models) {
    if (!MemoryManager::instance) {
        MemoryManager::instance = new MemoryManager(findNearestHyperisoDirectory() + lhaFile, models);
    }
    return MemoryManager::instance;
}

const std::string &MemoryManager::getData(std::string key)
{
    // if (cache.find(key) == cache.end()) {
    //     cache[key] = std::make_unique<std::string>(db->readBlock(key));
    // }
    return *cache[key];
}

LhaReader* MemoryManager::getReader() {
    return reader.get();
}

void MemoryManager::init() {
    
    // Init SM parameters from SMINPUTS
    reader = std::make_unique<LhaReader>(LhaReader(lhaPath));
    reader->readAll();

    Parameters::GetInstance(0);

    // Init other parameters

}