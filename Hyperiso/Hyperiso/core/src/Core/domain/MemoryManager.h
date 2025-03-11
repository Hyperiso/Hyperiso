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

#include <iostream>
#include <unordered_map>
#include <memory>
#include <cstdlib>
#include <thread>
#include <filesystem>
#include "config.hpp"
#include "General.h"
#include "Block.h"
#include "Parameters.h"
#include "ParamBlockAdapter.h"
#include "DBMemento.h"
#include "CorrelationRepo.h"
#include "CorrelationAdapter.h"
#include "DBManager.h"
#include "SpectrumCalculator.h"

struct DirPaths {
    static inline const fs::path default_dir_path           = project_assets_root.data() + std::string("default");
    static inline const fs::path savestate_dir_path         = project_assets_root.data() + std::string("savestate");
    static inline const fs::path spectrum_dir_path          = project_assets_root.data() + std::string("spectrum");
    static inline const fs::path template_dir_path          = project_assets_root.data() + std::string("template");
    static inline const fs::path param_mapping_dir_path     = project_assets_root.data() + std::string("input_files/marty_mapping");
};

struct FilePaths {
    static inline const fs::path default_obs_values_path   = DirPaths::default_dir_path/"observables.json";  ///< Path to observable covariance file
    static inline const fs::path default_obs_corr_path     = DirPaths::default_dir_path/"observables_corr.json";  ///< Path to observable covariance file
    static inline const fs::path default_param_values_path = DirPaths::default_dir_path/"parameters.json";
    static inline const fs::path default_param_corr_path   = DirPaths::default_dir_path/"parameters_corr.json";           ///< Path to parameter covariance file
    static inline const fs::path user_obs_values_path      = project_assets_root.data() + std::string("input_files/observables/observables.yaml");
    static inline const fs::path user_obs_corr_path        = project_assets_root.data() + std::string("input_files/observables/correlations.yaml");
    static inline const fs::path user_sm_params_path       = project_assets_root.data() + std::string("input_files/parameters/sm.yaml");
    static inline const fs::path user_flavor_params_path   = project_assets_root.data() + std::string("input_files/parameters/flavor.yaml");
    static inline const fs::path user_decay_params_path    = project_assets_root.data() + std::string("input_files/parameters/decay.yaml");
    static inline const fs::path user_param_corr_path      = project_assets_root.data() + std::string("input_files/parameters/correlations.yaml");
};

enum class ExternalFlag { IS_LHA_SPECTRUM, HAS_WILSON_INPUT, HAS_TH_OBSERVABLE_INPUT, USE_MARTY };
enum class InternalFlag { PARAMS_CHANGED };

struct Config {
    std::map<ExternalFlag, bool> flags;
    Model model {Model::SM};            ///< Model type (current model)
    std::optional<std::string> mty_model_name;  ///< MARTY model name (name of the class in MARTY) if needed
    std::optional<fs::path> mty_model_path;     ///< Path to the MARTY model file (mty_model_name.h) if needed
};

/**
 * @struct MemoryCache
 * @brief Stores memory cache details for parameter management.
 */
struct MemoryCache {
    Config config;                            ///< Config struct for various flags and runtime information
    std::map<InternalFlag, bool> flags;
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
    MemoryCache cache;                          ///< Internal memory cache
    static MemoryManager* instance;             ///< Singleton instance
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

    std::shared_ptr<BlockAccessor> read_input_files(fs::path lha_path);

    void read_default_input();
    void read_user_input();
    void read_lha_input(const std::string& lhaFile, const Config& config);

    fs::path calculate_spectrum(fs::path input_lha_path, const Config& config);

    /**
     * @brief Ensures the LHA path is correct.
     */
    fs::path format_lha_path(const std::string& path);

    void deduce_parameter_types(const Config &config);
    void save_input_cache();

public:
    /**
     * @brief Retrieves the singleton instance of MemoryManager.
     * @return Pointer to MemoryManager instance.
     */
    static MemoryManager* GetInstance();

    inline const MemoryCache& getMemoryCache() { check_if_ready(); return cache; }
    
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

    std::shared_ptr<BlockAccessor> extract_blocks(std::unordered_set<std::string> block_names);

    MemoryManager(const MemoryManager&) = delete;
    MemoryManager& operator=(const MemoryManager&) = delete;
    MemoryManager(MemoryManager&&) noexcept = default;
    MemoryManager& operator=(MemoryManager&&) noexcept = default;

    friend class Parameters;
};



#endif // HYPERISO_MEMORY_MANAGER_H