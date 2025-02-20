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
#include "lha_reader.h"
#include "config.hpp"
#include "General.h"
#include "Block.h"

/**
 * @struct MemoryCache
 * @brief Stores memory cache details for parameter management.
 */
struct MemoryCache {
    std::shared_ptr<LhaReader> reader;              ///< Shared pointer to LhaReader instance
    std::filesystem::path lha_path;                 ///< Path to LHA file
    std::filesystem::path obs_cov_path;             ///< Path to observable covariance file
    std::filesystem::path param_cov_path;           ///< Path to parameter covariance file
    std::vector<ParameterType> parameter_types;     ///< List of parameter types available
    std::thread::id thread_id;                      ///< ID of the thread using the cache
    Model model;                                    ///< Model type (current model)
    bool is_spectrum;                               ///< Indicates if spectrum is available in the lha
    bool has_wilsons;                               ///< Indicates if Wilson coefficients are included in the lha
    bool has_obs;                                   ///< Indicates if observables are present
    bool use_marty;                                 ///< Indicates if Marty framework is used for Wilson calculation
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
    MemoryCache cache;              ///< Internal memory cache
    static MemoryManager* instance; ///< Singleton instance

    /**
     * @brief Private constructor to enforce singleton pattern.
     */
    MemoryManager();

    /**
     * @brief Checks if the memory manager is ready.
     * @throws std::logic_error If the memory manager is not initialized.
     */
    void check_if_ready();

public:
    /**
     * @brief Retrieves the singleton instance of MemoryManager.
     * @return Pointer to MemoryManager instance.
     */
    static MemoryManager* GetInstance();

    /**
     * @brief Retrieves the LhaReader instance.
     * @return Pointer to LhaReader instance.
     */
    LhaReader* getReader();

    inline bool isSpectrum() { check_if_ready(); return cache.is_spectrum; }
    inline bool hasWilsons() { check_if_ready(); return cache.has_wilsons; }
    inline bool hasObservables() { check_if_ready(); return cache.has_obs; }
    inline bool useMarty() { check_if_ready(); return cache.use_marty; }
    inline bool paramCacheOkay() { check_if_ready(); return cache.param_cache_okay; }
    inline std::filesystem::path getInputLhaPath() { check_if_ready(); return cache.lha_path; }
    inline std::filesystem::path getParameterCovariancePath() { check_if_ready(); return cache.param_cov_path; }
    inline std::filesystem::path getObservableCovariancePath()  { check_if_ready(); return cache.obs_cov_path; }
    inline std::vector<ParameterType> getParameterTypes() { check_if_ready(); return cache.parameter_types; };
    inline Model getModel() { check_if_ready(); return cache.model; };
    inline bool getUseMarty() {check_if_ready(); return cache.use_marty;}
    
    /**
     * @brief Initializes the memory manager with LHA file and model settings. First method to use, mandatory.
     */
    void init(const std::string& lhaFile, Model model = Model::SM, bool use_marty = false, bool is_spectrum=false, bool has_wilsons=false, bool has_obs=false);

    /**
     * @brief Switches the LHA file and reinitializes the manager.
     */
    void switch_lha(const std::string& lhaFile, Model model = Model::SM, bool use_marty = false, bool is_spectrum = false, bool has_wilson = false, bool has_obs = false);
    
    /**
     * @brief Switches the model being used. Need to be compatible with the LHA.
     */
    void switch_model(Model model = Model::SM, bool use_marty = false);

    /**
     * @brief Sets the observable covariance file path.
     */
    void set_observable_covariance_input_file(const std::string& path);

    /**
     * @brief Sets the parameter covariance file path.
     */
    void set_parameter_covariance_input_file(const std::string& path);

    /**
     * @brief Retrieves a list of parameter blocks.
     */
    std::vector<std::string> get_blocks_list(ParameterType param_type = ParameterType::SM);

    /**
     * @brief Retrieves block information for a given block.
     */
    std::map<int, double> get_block_infos(const std::string& block, ParameterType param_type = ParameterType::SM);

    /**
     * @brief Retrieves the parameter types associated with a block.
     */
    std::vector<ParameterType> get_type_of_block(const std::string& block);
    
    MemoryManager(const MemoryManager&) = delete;
    MemoryManager& operator=(const MemoryManager&) = delete;
    MemoryManager(MemoryManager&&) noexcept = default;
    MemoryManager& operator=(MemoryManager&&) noexcept = default;
};



#endif // HYPERISO_MEMORY_MANAGER_H