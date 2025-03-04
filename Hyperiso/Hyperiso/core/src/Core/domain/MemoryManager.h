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

namespace fs = std::filesystem;

struct DirPaths {
    static inline const fs::path default_dir_path   = project_assets_root.data() + std::string("default");
    static inline const fs::path savestate_dir_path = project_assets_root.data() + std::string("savestate");
    static inline const fs::path spectrum_dir_path  = project_assets_root.data() + std::string("spectrum");
    static inline const fs::path template_dir_path  = project_assets_root.data() + std::string("template");
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

struct Config {
    bool is_spectrum    {false};                ///< Indicates if spectrum is available in the lha
    bool has_wilsons    {false};                ///< Indicates if Wilson coefficients are included in the lha
    bool has_obs        {false};                ///< Indicates if observables are present
    bool use_marty      {false};                ///< Indicates if Marty framework is used for Wilson calculation
    Model model         {Model::SM};            ///< Model type (current model)
    std::optional<std::string> mty_model_name;  ///< MARTY model name (name of the class in MARTY) if needed
    std::optional<fs::path> mty_model_path;     ///< Path to the MARTY model file (mty_model_name.h) if needed
};

/**
 * @struct MemoryCache
 * @brief Stores memory cache details for parameter management.
 */
struct MemoryCache {
    Config config;                            ///< Config struct for various flags and runtime information
    fs::path lha_path;                              ///< Path to LHA file
    std::vector<ParameterType> parameter_types;     ///< List of parameter types available
    std::thread::id thread_id;                      ///< ID of the thread using the cache
    bool is_ready;                                  ///< Indicates if cache is ready for use
    bool param_cache_okay;                          ///< Indicates if parameter cache is valid
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

public:
    /**
     * @brief Retrieves the singleton instance of MemoryManager.
     * @return Pointer to MemoryManager instance.
     */
    static MemoryManager* GetInstance();

    inline bool isSpectrum() { check_if_ready(); return cache.config.is_spectrum; }
    inline bool hasWilsons() { check_if_ready(); return cache.config.has_wilsons; }
    inline bool hasObservables() { check_if_ready(); return cache.config.has_obs; }
    inline bool paramCacheOkay() { check_if_ready(); return cache.param_cache_okay; }
    inline fs::path getInputLhaPath() { check_if_ready(); return cache.lha_path; }
    inline std::vector<ParameterType> getParameterTypes() { check_if_ready(); return cache.parameter_types; };
    inline Model getModel() { check_if_ready(); return cache.config.model; };
    inline bool getUseMarty() {check_if_ready(); return cache.config.use_marty;}
    
    /**
     * @brief Initializes the memory manager with LHA file and model settings. First method to use, mandatory.
     */
    void init(const std::string &lhaFile, const Config& config);

    void deduce_parameter_types(const Config &config);

    /**
     * @brief Switches the LHA file and reinitializes the manager.
     */
    void switch_lha(const std::string& lhaFile, Config config);
    
    /**
     * @brief Ensures the LHA path is correct.
     */
    fs::path format_lha_path(const std::string& path);

    /**
     * @brief Switches the model being used. Need to be compatible with the LHA.
     */
    void switch_model(Model model = Model::SM, bool use_marty = false);

    /**
     * @brief Retrieves a list of parameter blocks.
     */
    std::unordered_set<std::string> get_blocks_list(ParameterType param_type = ParameterType::SM);

    /**
     * @brief Retrieves the list of all parameter blocks stored.
     */
    std::unordered_set<std::string> get_all_blocks();

    /**
     * @brief Retrieves block information for a given block.
     */
    std::map<LhaID, double> get_block_infos(const std::string& block, ParameterType param_type = ParameterType::SM);

    /**
     * @brief Retrieves the parameter types associated with a block.
     */
    std::vector<ParameterType> get_type_of_block(const std::string& block);

    std::shared_ptr<BlockAccessor> get_blocks(std::unordered_set<std::string> block_names);
    
    void save_input_cache();

    MemoryManager(const MemoryManager&) = delete;
    MemoryManager& operator=(const MemoryManager&) = delete;
    MemoryManager(MemoryManager&&) noexcept = default;
    MemoryManager& operator=(MemoryManager&&) noexcept = default;

    friend class Parameters;
};



#endif // HYPERISO_MEMORY_MANAGER_H