#ifndef HYPERISO_PATHS_H
#define HYPERISO_PATHS_H

#include <filesystem>

namespace fs = std::filesystem;

/**
 * @brief Relative asset paths only.
 *
 * Runtime roots must be provided by IPathsProvider. Do not place generated
 * absolute roots here; they are not relocatable once packaged in a Python wheel.
 */
struct AssetRelativePaths {
    static inline const fs::path default_dir       = "default";
    static inline const fs::path user_dir          = "input_files";
    static inline const fs::path template_dir      = "template";
    static inline const fs::path marty_template_dir= fs::path("template") / "MARTY";
    static inline const fs::path marty_mapping_dir = fs::path("input_files") / "marty_mapping";
};

struct DefaultInputRelativePaths {
    static inline const fs::path obs_values   = fs::path("default") / "observables.json";
    static inline const fs::path obs_corr     = fs::path("default") / "observables_corr.json";
    static inline const fs::path param_values = fs::path("default") / "parameters.json";
    static inline const fs::path param_corr   = fs::path("default") / "parameters_corr.json";
    static inline const fs::path nuisances    = fs::path("default") / "nuisances.json";
};

struct UserInputRelativePaths {
    static inline const fs::path obs_values    = fs::path("input_files") / "observables" / "observables.yaml";
    static inline const fs::path obs_corr      = fs::path("input_files") / "observables" / "correlations.yaml";
    static inline const fs::path sm_params     = fs::path("input_files") / "parameters" / "sm.yaml";
    static inline const fs::path flavor_params = fs::path("input_files") / "parameters" / "flavor.yaml";
    static inline const fs::path decay_params  = fs::path("input_files") / "parameters" / "decay.yaml";
    static inline const fs::path param_corr    = fs::path("input_files") / "parameters" / "correlations.yaml";
    static inline const fs::path nuisances     = fs::path("input_files") / "parameters" / "nuisances.yaml";
};

struct CacheRelativePaths {
    static inline const fs::path spectrum_dir  = "Spectrum";
    static inline const fs::path marty_temp_dir= "MartyTemp";
};

#endif // HYPERISO_PATHS_H
