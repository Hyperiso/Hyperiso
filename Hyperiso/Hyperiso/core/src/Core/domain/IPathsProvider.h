#ifndef HYPERISO_IPATHSPROVIDER_H
#define HYPERISO_IPATHSPROVIDER_H

#include <filesystem>

namespace fs = std::filesystem;

/**
 * @enum APIPath
 * @brief Enumerates the filesystem paths exposed through the public API.
 *
 * LHA_PATH is a runtime path stored by MemoryManager after init(). All other
 * entries are provided by IPathsProvider and may be overridden before init()
 * through HyperisoMaster::pre_init_set_paths().
 */
enum class APIPath {
    LHA_PATH,                 ///< Path to the active LHA file.

    ASSETS_ROOT,              ///< Root directory for HyperISO read-only assets.

    DEFAULT_PARAM_VALUES,     ///< JSON file containing default parameter values.
    DEFAULT_OBS_VALUES,       ///< JSON file containing default observable values.
    DEFAULT_PARAM_CORR,       ///< JSON file containing default parameter correlations.
    DEFAULT_OBS_CORR,         ///< JSON file containing default observable correlations.
    DEFAULT_NUISANCES,        ///< JSON file containing default nuisance definitions.

    USER_SM_PARAMS,           ///< YAML/YML file containing user SM parameter overrides.
    USER_FLAVOR_PARAMS,       ///< YAML/YML file containing user flavor parameter overrides.
    USER_DECAY_PARAMS,        ///< YAML/YML file containing user decay parameter overrides.
    USER_OBS_VALUES,          ///< YAML/YML file containing user observable overrides.
    USER_PARAM_CORR,          ///< YAML/YML file containing user parameter correlation overrides.
    USER_OBS_CORR,            ///< YAML/YML file containing user observable correlation overrides.
    USER_NUISANCES,           ///< YAML/YML file containing user nuisance overrides.

    PARAM_MAPPING_DIR,        ///< Read-only directory containing MARTY/Hyperiso parameter mappings.
    TEMPLATE_DIR,             ///< Read-only directory containing generated-code templates.

    SPECTRUM_DIR,             ///< Writable cache directory used for generated spectrum files.
    MARTY_TEMP_DIR            ///< Writable cache directory used for generated MARTY files.
};

/**
 * @brief Interface providing all filesystem paths used by HyperISO.
 *
 * Implementations define where read-only package assets and writable runtime
 * cache directories are located on disk.
 */
struct IPathsProvider {
    virtual ~IPathsProvider() = default;

    /// Root directory for all read-only assets.
    virtual fs::path assets_root() const = 0;

    // read-only defaults
    virtual fs::path default_param_values() const = 0;
    virtual fs::path default_obs_values()   const = 0;
    virtual fs::path default_param_corr()   const = 0;
    virtual fs::path default_obs_corr()     const = 0;
    virtual fs::path default_nuisances()    const = 0;

    // read-only packaged user-input templates/default overrides
    virtual fs::path user_sm_params()       const = 0;
    virtual fs::path user_flavor_params()   const = 0;
    virtual fs::path user_decay_params()    const = 0;
    virtual fs::path user_obs_values()      const = 0;
    virtual fs::path user_param_corr()      const = 0;
    virtual fs::path user_obs_corr()        const = 0;
    virtual fs::path user_nuisances()       const = 0;

    // writable output/cache directories
    virtual fs::path spectrum_dir()         const = 0;
    virtual fs::path marty_temp_dir()       const = 0;

    // read-only MARTY resources
    virtual fs::path param_mapping_dir_path() const = 0;
    virtual fs::path template_dir_path() const = 0;
};

#endif
