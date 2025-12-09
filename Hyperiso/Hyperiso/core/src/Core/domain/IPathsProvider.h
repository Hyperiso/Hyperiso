#ifndef HYPERISO_IPATHSPROVIDER_H
#define HYPERISO_IPATHSPROVIDER_H

#include <filesystem>

namespace fs = std::filesystem;

/**
 * @brief Interface providing all filesystem paths used by HyperISO.
 *
 * Implementations define where default/user inputs and outputs (spectrum, templates, mappings)
 * are located on disk.
 */
struct IPathsProvider {
    virtual ~IPathsProvider() = default;

    /// Root directory for all assets.
    virtual fs::path assets_root() const = 0;

    // defaults
    virtual fs::path default_param_values() const = 0;
    virtual fs::path default_obs_values()   const = 0;
    virtual fs::path default_param_corr()   const = 0;
    virtual fs::path default_obs_corr()     const = 0;

    // user inputs
    virtual fs::path user_sm_params()       const = 0;
    virtual fs::path user_flavor_params()   const = 0;
    virtual fs::path user_decay_params()    const = 0;
    virtual fs::path user_obs_values()      const = 0;
    virtual fs::path user_param_corr()      const = 0;
    virtual fs::path user_obs_corr()      const = 0;

    // spectrum output dir
    virtual fs::path spectrum_dir()         const = 0;
    virtual fs::path param_mapping_dir_path() const = 0;
    virtual fs::path template_dir_path() const = 0;
};

#endif