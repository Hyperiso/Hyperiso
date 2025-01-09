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

struct MemoryCache {
    std::shared_ptr<LhaReader> reader;
    std::filesystem::path lha_path;
    std::filesystem::path obs_cov_path;
    std::filesystem::path param_cov_path;
    std::vector<ParameterType> parameter_types;
    std::thread::id thread_id;
    Model model;
    bool is_spectrum;
    bool has_wilsons;
    bool has_obs;
    bool use_marty;
    bool is_ready;
    bool param_cache_okay;
};

class MemoryManager {
private:
    MemoryCache cache;
    static MemoryManager* instance;
    MemoryManager();

    void check_if_ready();

public:

    static std::string findNearestHyperisoDirectory();
    static MemoryManager* GetInstance();

    // template<typename T>
    // std::shared_ptr<T, void(*)(void*)> makeUniquePtr(T* ptr);

    // template<typename T>
    // std::shared_ptr<T, void(*)(void*)> allocate();

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

    void init(const std::string& lhaFile, Model model = Model::SM, bool use_marty = false, bool is_spectrum=false, bool has_wilsons=false, bool has_obs=false);

    void switch_lha(const std::string& lhaFile, Model model = Model::SM, bool use_marty = false, bool is_spectrum = false, bool has_wilson = false, bool has_obs = false);
    void switch_model(Model model = Model::SM, bool use_marty = false);

    void set_observable_covariance_input_file(const std::string& path);
    void set_parameter_covariance_input_file(const std::string& path);

    std::vector<std::string> get_blocks_list(ParameterType param_type = ParameterType::SM);

    std::map<int, double> get_block_infos(const std::string& block, ParameterType param_type = ParameterType::SM);

    MemoryManager(const MemoryManager&) = delete;
    MemoryManager& operator=(const MemoryManager&) = delete;
    MemoryManager(MemoryManager&&) noexcept = default;
    MemoryManager& operator=(MemoryManager&&) noexcept = default;
};



#endif // HYPERISO_MEMORY_MANAGER_H