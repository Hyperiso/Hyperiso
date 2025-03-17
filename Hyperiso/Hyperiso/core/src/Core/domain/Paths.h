#ifndef __PATHS_H__
#define __PATHS_H__

#include <filesystem>
#include "config.hpp"

namespace fs = std::filesystem;

struct DirPaths {
    static inline const fs::path default_dir_path           = project_assets_root.data() + std::string("default");
    static inline const fs::path user_dir_path              = project_assets_root.data() + std::string("input_files");
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
    static inline const fs::path user_obs_values_path      = DirPaths::user_dir_path/"observables/observables.yaml";
    static inline const fs::path user_obs_corr_path        = DirPaths::user_dir_path/"observables/correlations.yaml";
    static inline const fs::path user_sm_params_path       = DirPaths::user_dir_path/"parameters/sm.yaml";
    static inline const fs::path user_flavor_params_path   = DirPaths::user_dir_path/"parameters/flavor.yaml";
    static inline const fs::path user_decay_params_path    = DirPaths::user_dir_path/"parameters/decay.yaml";
    static inline const fs::path user_param_corr_path      = DirPaths::user_dir_path/"parameters/correlations.yaml";
};

#endif // __PATHS_H__
