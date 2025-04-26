/**
 * @file MemoryManager.h
 * @brief Manages memory caching, parameter blocks, and LHA reader instances.
 *
 * This file defines the MemoryManager singleton class, responsible for managing:
 * - Memory caching (parameters, blocks, correlations)
 * - LHA file loading and reloading
 * - Spectrum calculations
 * 
 * It must be initialized before using any parameter-dependent functionality.
 */
#ifndef HYPERISO_MEMORY_MANAGER_H
#define HYPERISO_MEMORY_MANAGER_H

#include <unordered_map>
#include <thread>

#include "Include.h"
#include "BlockAccessor.h"
#include "ParamBlockLoader.h"
#include "DBMemento.h"
#include "CorrelationRepo.h"
#include "CorrelationAdapter.h"
#include "SpectrumCalculator.h"
#include "Paths.h"
#include "Config.h"

/**
 * @enum InternalFlag
 * @brief Internal flags used for memory management status.
 */
enum class InternalFlag { PARAMS_CHANGED };

/**
 * @struct MemoryCache
 * @brief Stores memory cache details for parameter management and runtime state.
 */
struct MemoryCache {
    Config config;                                                              ///< Config struct for various flags and runtime information
    std::map<InternalFlag, bool> flags {{InternalFlag::PARAMS_CHANGED, false}}; ///< Internal status flags.
    fs::path lha_path;                                                          ///< Path to the currently loaded LHA file.
    std::vector<ParameterType> parameter_types;                                 ///< List of parameter types currently managed.
    std::thread::id thread_id;                                                  ///< Thread ID associated with the current cache usage.
    bool is_ready;                                                              ///< Indicates if the memory manager is initialized and ready.

};

/**
 * @class MemoryManager
 * @brief Singleton class responsible for initializing and managing memory, input files, and parameter blocks.
 *
 * MemoryManager must be initialized first in a script or application that uses parameters.
 */
class MemoryManager {
private:
    static MemoryManager* instance;                 ///< Singleton instance pointer.
    MemoryCache cache;                              ///< Internal memory cache for session data.
    std::shared_ptr<BlockAccessor> input_cache;     ///< Cached parameters and blocks.
    DBMemento memento;                              ///< Memento to save and restore input cache.
    CorrelationRepository correlation_repository;   ///< Repository for parameter and observable correlations.

    /**
     * @brief Private constructor to enforce singleton pattern.
     */
    MemoryManager();

    /**
     * @brief Ensures the memory manager is initialized before usage.
     * @throws std::logic_error If memory manager is not ready.
     */
    void check_if_ready();

    /**
     * @brief Reads and loads default input files (parameters, observables, correlations).
     *
     * This method is used during the initialization process to load standard inputs
     * from predefined file paths into the input cache.
     */
    void read_default_input();

    /**
     * @brief Reads and loads user-specific input files (overrides or complements defaults).
     *
     * Loads user-provided parameter and observable blocks from input files, merging them
     * into the current input cache.
     */
    void read_user_input();

    /**
     * @brief Reads a LHA input file and updates the memory cache accordingly.
     *
     * @param lhaFile Path to the LHA input file.
     * @param config Configuration options for the reading process.
     */
    void read_lha_input(const std::string& lhaFile, const Config& config);

    /**
     * @brief Formats and resolves the full path to an LHA input file.
     *
     * Converts relative paths into absolute paths based on the assets directory.
     * Also verifies the existence of the file.
     *
     * @param path Relative or absolute path to the LHA file.
     * @return Resolved absolute path to the LHA file.
     */
    fs::path format_lha_path(const std::string& path);

    /**
     * @brief Calculates a new LHA file if a spectrum needs to be generated.
     *
     * If necessary (depending on the model and config flags), invokes the external spectrum calculator.
     *
     * @param input_lha_path Path to the initial LHA file.
     * @param config Configuration options to control spectrum generation.
     * @return Path to the new LHA file containing the calculated spectrum.
     */
    fs::path calculate_spectrum(fs::path input_lha_path, const Config& config);

    /**
     * @brief Determines the list of parameter types to manage based on the configuration.
     *
     * Fills the internal list of ParameterTypes (`cache.parameter_types`) based on the model
     * and flags present in the Config object.
     *
     * @param config Configuration options to interpret.
     */
    void deduce_parameter_types(const Config &config);

    /**
     * @brief Saves a snapshot of the current input cache.
     *
     * Stores the current state of the input cache to allow future restoration
     * (used when switching LHA files or reloading user input).
     */
    void save_input_cache();

public:
    /**
     * @brief Retrieves the singleton instance of MemoryManager.
     * @return Pointer to MemoryManager instance.
     */
    static MemoryManager* GetInstance();
    
    /**
     * @brief Initializes the memory manager with the provided LHA file and configuration.
     *
     * Must be called before using any Parameters instances.
     *
     * @param lhaFile Path to the input LHA file.
     * @param config Configuration settings to use.
     */
    void init(const std::string &lhaFile, Config config);

    /**
     * @brief Switches to a different LHA file, reloading associated parameters and spectrum.
     * @param lhaFile Path to the new LHA file.
     * @param config New configuration to apply.
     */
    void switch_lha(const std::string& lhaFile, Config config);

    /**
     * @brief Reloads user-specific input files.
     *
     * Useful for reapplying user customizations without reloading the full LHA file.
     * @param config New configuration to apply.
     */
    void reload_user_input(Config config);

    /**
     * @brief Reloads user-specific input files with a different LHA file.
     * @param lhaFile Path to the LHA file.
     * @param config New configuration to apply.
     */
    void reload_user_input(const std::string &lhaFile, Config config);

    /**
     * @brief Switches the model used for spectrum and parameters.
     *
     * Must be compatible with the loaded LHA file.
     *
     * @param model The new model to switch to (default is SM).
     * @param use_marty Whether to use external MARTY spectrum calculation.
     */
    void switch_model(Model model = Model::SM, bool use_marty = false);

    /**
     * @brief Retrieves the current memory cache.
     * 
     * @return A const reference to the internal MemoryCache.
     */
    inline const MemoryCache& getMemoryCache() { check_if_ready(); return cache; }

    /**
     * @brief Extracts specific blocks from the cached input.
     *
     * @param block_names Set of block names to extract.
     * @return Shared pointer to a new BlockAccessor containing the requested blocks.
     */
    std::shared_ptr<BlockAccessor> extract_blocks(std::unordered_set<std::string> block_names);

    MemoryManager(const MemoryManager&) = delete;
    MemoryManager& operator=(const MemoryManager&) = delete;
    MemoryManager(MemoryManager&&) noexcept = default;
    MemoryManager& operator=(MemoryManager&&) noexcept = default;

    friend class Parameters;
};

#endif // HYPERISO_MEMORY_MANAGER_H