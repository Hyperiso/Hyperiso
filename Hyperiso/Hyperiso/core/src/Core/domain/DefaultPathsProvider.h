#ifndef HYPERISO_DEFAULTPATHSPROVIDER_H
#define HYPERISO_DEFAULTPATHSPROVIDER_H

#include <map>
#include <utility>

#include "IPathsProvider.h"
#include "Paths.h"

/**
 * @brief Default implementation of IPathsProvider backed by generated project paths.
 */
struct DefaultPathsProvider : public IPathsProvider {
    fs::path assets_root() const override { return fs::path(project_assets_root.data()); }

    fs::path default_param_values() const override { return FilePaths::default_param_values_path; }
    fs::path default_obs_values()   const override { return FilePaths::default_obs_values_path; }
    fs::path default_param_corr()   const override { return FilePaths::default_param_corr_path; }
    fs::path default_obs_corr()     const override { return FilePaths::default_obs_corr_path; }

    fs::path user_sm_params()       const override { return FilePaths::user_sm_params_path; }
    fs::path user_flavor_params()   const override { return FilePaths::user_flavor_params_path; }
    fs::path user_decay_params()    const override { return FilePaths::user_decay_params_path; }
    fs::path user_obs_values()      const override { return FilePaths::user_obs_values_path; }
    fs::path user_param_corr()      const override { return FilePaths::user_param_corr_path; }
    fs::path user_obs_corr()        const override { return FilePaths::user_obs_corr_path; }

    fs::path param_mapping_dir_path() const override { return DirPaths::param_mapping_dir_path; }
    fs::path template_dir_path()      const override { return DirPaths::template_dir_path; }
    fs::path spectrum_dir()          const override { return DirPaths::spectrum_dir_path; }
};

/**
 * @brief Default path provider with validated user overrides.
 *
 * This provider keeps the generated defaults for every path that is not present
 * in the override map. It is intended to be installed through
 * HyperisoMaster::pre_init_set_paths() before MemoryManager::init().
 */
struct OverridePathsProvider : public DefaultPathsProvider {
    explicit OverridePathsProvider(std::map<APIPath, fs::path> path_overrides)
        : overrides(std::move(path_overrides)) {}

    fs::path assets_root() const override {
        return get_or_default(APIPath::ASSETS_ROOT, DefaultPathsProvider::assets_root());
    }

    fs::path default_param_values() const override {
        return get_or_default(APIPath::DEFAULT_PARAM_VALUES, DefaultPathsProvider::default_param_values());
    }

    fs::path default_obs_values() const override {
        return get_or_default(APIPath::DEFAULT_OBS_VALUES, DefaultPathsProvider::default_obs_values());
    }

    fs::path default_param_corr() const override {
        return get_or_default(APIPath::DEFAULT_PARAM_CORR, DefaultPathsProvider::default_param_corr());
    }

    fs::path default_obs_corr() const override {
        return get_or_default(APIPath::DEFAULT_OBS_CORR, DefaultPathsProvider::default_obs_corr());
    }

    fs::path user_sm_params() const override {
        return get_or_default(APIPath::USER_SM_PARAMS, DefaultPathsProvider::user_sm_params());
    }

    fs::path user_flavor_params() const override {
        return get_or_default(APIPath::USER_FLAVOR_PARAMS, DefaultPathsProvider::user_flavor_params());
    }

    fs::path user_decay_params() const override {
        return get_or_default(APIPath::USER_DECAY_PARAMS, DefaultPathsProvider::user_decay_params());
    }

    fs::path user_obs_values() const override {
        return get_or_default(APIPath::USER_OBS_VALUES, DefaultPathsProvider::user_obs_values());
    }

    fs::path user_param_corr() const override {
        return get_or_default(APIPath::USER_PARAM_CORR, DefaultPathsProvider::user_param_corr());
    }

    fs::path user_obs_corr() const override {
        return get_or_default(APIPath::USER_OBS_CORR, DefaultPathsProvider::user_obs_corr());
    }

    fs::path param_mapping_dir_path() const override {
        return get_or_default(APIPath::PARAM_MAPPING_DIR, DefaultPathsProvider::param_mapping_dir_path());
    }

    fs::path template_dir_path() const override {
        return get_or_default(APIPath::TEMPLATE_DIR, DefaultPathsProvider::template_dir_path());
    }

    fs::path spectrum_dir() const override {
        return get_or_default(APIPath::SPECTRUM_DIR, DefaultPathsProvider::spectrum_dir());
    }

private:
    std::map<APIPath, fs::path> overrides;

    fs::path get_or_default(APIPath path_name, const fs::path& default_path) const {
        const auto it = overrides.find(path_name);
        return it == overrides.end() ? default_path : it->second;
    }
};

#endif
