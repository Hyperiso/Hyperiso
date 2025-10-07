#ifndef HYPERISO_DEFAULTPATHSPROVIDER_H
#define HYPERISO_DEFAULTPATHSPROVIDER_H

#include "IPathsProvider.h"
#include "Paths.h"

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

    fs::path param_mapping_dir_path()         const override { return DirPaths::param_mapping_dir_path; }
    fs::path template_dir_path()         const override { return DirPaths::template_dir_path; }
    fs::path spectrum_dir()         const override { return DirPaths::spectrum_dir_path; }
};

#endif