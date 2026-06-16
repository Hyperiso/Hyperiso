#ifndef HYPERISO_DEFAULTPATHSPROVIDER_H
#define HYPERISO_DEFAULTPATHSPROVIDER_H

#include <cstdlib>
#include <map>
#include <utility>

#include "config.hpp"
#include "IPathsProvider.h"
#include "Paths.h"

/**
 * @brief Default implementation of IPathsProvider.
 *
 * Read-only assets are rooted at project_assets_root by default, but may be
 * overridden before init() through OverridePathsProvider. Writable runtime data
 * goes to the user cache directory by default, not under the assets tree.
 */
struct DefaultPathsProvider : public IPathsProvider {
    fs::path assets_root() const override { return fs::path(project_assets_root.data()); }

    fs::path default_param_values() const override { return assets_root() / DefaultInputRelativePaths::param_values; }
    fs::path default_obs_values()   const override { return assets_root() / DefaultInputRelativePaths::obs_values; }
    fs::path default_param_corr()   const override { return assets_root() / DefaultInputRelativePaths::param_corr; }
    fs::path default_obs_corr()     const override { return assets_root() / DefaultInputRelativePaths::obs_corr; }
    fs::path default_nuisances()    const override { return assets_root() / DefaultInputRelativePaths::nuisances; }

    fs::path user_sm_params()       const override { return assets_root() / UserInputRelativePaths::sm_params; }
    fs::path user_flavor_params()   const override { return assets_root() / UserInputRelativePaths::flavor_params; }
    fs::path user_decay_params()    const override { return assets_root() / UserInputRelativePaths::decay_params; }
    fs::path user_obs_values()      const override { return assets_root() / UserInputRelativePaths::obs_values; }
    fs::path user_param_corr()      const override { return assets_root() / UserInputRelativePaths::param_corr; }
    fs::path user_obs_corr()        const override { return assets_root() / UserInputRelativePaths::obs_corr; }
    fs::path user_nuisances()       const override { return assets_root() / UserInputRelativePaths::nuisances; }

    fs::path param_mapping_dir_path() const override { return assets_root() / AssetRelativePaths::marty_mapping_dir; }
    fs::path template_dir_path()      const override { return assets_root() / AssetRelativePaths::template_dir; }

    fs::path spectrum_dir()          const override { return default_cache_root() / CacheRelativePaths::spectrum_dir; }
    fs::path marty_temp_dir()        const override { return default_cache_root() / CacheRelativePaths::marty_temp_dir; }

protected:
    static fs::path default_cache_root() {
        if (const char* env = std::getenv("HYPERISO_CACHE_ROOT")) {
            if (*env) return fs::path(env);
        }

        if (const char* xdg = std::getenv("XDG_CACHE_HOME")) {
            if (*xdg) return fs::path(xdg) / "pyhyperiso";
        }

        if (const char* home = std::getenv("HOME")) {
            if (*home) return fs::path(home) / ".cache" / "pyhyperiso";
        }

        return fs::temp_directory_path() / "pyhyperiso";
    }
};

/**
 * @brief Default path provider with validated user overrides.
 *
 * This provider keeps the generated/detected defaults for every path that is
 * not present in the override map. It is intended to be installed through
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

    fs::path default_nuisances() const override {
        return get_or_default(APIPath::DEFAULT_NUISANCES, DefaultPathsProvider::default_nuisances());
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

    fs::path user_nuisances() const override {
        return get_or_default(APIPath::USER_NUISANCES, DefaultPathsProvider::user_nuisances());
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

    fs::path marty_temp_dir() const override {
        return get_or_default(APIPath::MARTY_TEMP_DIR, DefaultPathsProvider::marty_temp_dir());
    }

private:
    std::map<APIPath, fs::path> overrides;

    fs::path get_or_default(APIPath path_name, const fs::path& default_path) const {
        const auto it = overrides.find(path_name);
        return it == overrides.end() ? default_path : it->second;
    }
};

#endif
