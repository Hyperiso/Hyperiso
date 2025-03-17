/**
 * @file MemoryManager.h
 * @brief Manages memory caching and parameter storage.
 * 
 * This file defines the MemoryManager class, which is responsible for managing
 * memory caches, parameter blocks, and LHA reader instances.
 * This class need to be called before doing anything else.
 */
#if !defined(HYPERISO_MEMORY_MANAGER_H)
#define HYPERISO_MEMORY_MANAGER_H

#include <unordered_map>
#include <thread>
#include "Include.h"
#include "BlockAccessor.h"
#include "ParamBlockAdapter.h"
#include "DBMemento.h"
#include "CorrelationRepo.h"
#include "CorrelationAdapter.h"
#include "SpectrumCalculator.h"
#include "Paths.h"
#include "Config.h"

enum class InternalFlag { PARAMS_CHANGED };

/**
 * @struct MemoryCache
 * @brief Stores memory cache details for parameter management.
 */
struct MemoryCache {
    Config config;                            ///< Config struct for various flags and runtime information
    std::map<InternalFlag, bool> flags {
        {InternalFlag::PARAMS_CHANGED, false},
    };
    fs::path lha_path;                              ///< Path to LHA file
    std::vector<ParameterType> parameter_types;     ///< List of parameter types available
    std::thread::id thread_id;                      ///< ID of the thread using the cache
    bool is_ready;                                  ///< Indicates if cache is ready for use
};

/**
 * @class MemoryManager
 * @brief Singleton class for managing memory and LHA reader instances.
 * First classe to use in a script, before using Parameters.
 */
class MemoryManager {
private:
    static MemoryManager* instance;             ///< Singleton instance
    MemoryCache cache;                          ///< Internal memory cache
    std::shared_ptr<BlockAccessor> input_cache; ///< BlockAccessor filled with all the blocks read from input files
    DBMemento memento;
    CorrelationRepository correlation_repository;

    /**
     * @brief Private constructor to enforce singleton pattern.
     */
    MemoryManager();

    /**
     * @brief Checks if the memory manager is ready.
     * @throws std::logic_error If the memory manager is not initialized.
     */
    void check_if_ready();

    void read_default_input();
    void read_user_input();
    void read_lha_input(const std::string& lhaFile, const Config& config);
    fs::path format_lha_path(const std::string& path);
    fs::path calculate_spectrum(fs::path input_lha_path, const Config& config);
    void deduce_parameter_types(const Config &config);
    void save_input_cache();

public:
    /**
     * @brief Retrieves the singleton instance of MemoryManager.
     * @return Pointer to MemoryManager instance.
     */
    static MemoryManager* GetInstance();
    
    /**
     * @brief Initializes the memory manager with LHA file and model settings. First method to use, mandatory.
     */
    void init(const std::string &lhaFile, Config config);

    /**
     * @brief Switches the LHA file and reinitializes the manager.
     */
    void switch_lha(const std::string& lhaFile, Config config);

    void reload_user_input(Config config);

    void reload_user_input(const std::string &lhaFile, Config config);

    /**
     * @brief Switches the model being used. Need to be compatible with the LHA.
     */
    void switch_model(Model model = Model::SM, bool use_marty = false);

    inline const MemoryCache& getMemoryCache() { check_if_ready(); return cache; }
    std::shared_ptr<BlockAccessor> extract_blocks(std::unordered_set<std::string> block_names);

    MemoryManager(const MemoryManager&) = delete;
    MemoryManager& operator=(const MemoryManager&) = delete;
    MemoryManager(MemoryManager&&) noexcept = default;
    MemoryManager& operator=(MemoryManager&&) noexcept = default;

    friend class Parameters;
};

#endif // HYPERISO_MEMORY_MANAGER_H