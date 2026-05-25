#ifndef DEFAULTNUISANCEPATHSPROVIDER_H
#define DEFAULTNUISANCEPATHSPROVIDER_H

#include "INuisancePathsProvider.h"
#include "Paths.h"

/**
 * @file DefaultNuisancePathsProvider.h
 * @brief Default nuisance-configuration path provider.
 *
 * @see INuisancePathsProvider
 * @see NuisanceReader
 */

/**
 * @struct DefaultNuisancePathsProvider
 * @brief Provides the standard application paths for nuisance configuration.
 *
 * The default registry is expected in the application default directory, while
 * user overrides are expected under the user parameter directory.
 */
struct DefaultNuisancePathsProvider : public INuisancePathsProvider {
    /**
     * @copydoc INuisancePathsProvider::default_nuisances_path
     *
     * @return `DirPaths::default_dir_path / "nuisances.json"`.
     */
    fs::path default_nuisances_path() const override {
        return DirPaths::default_dir_path / "nuisances.json";
    }

    /**
     * @copydoc INuisancePathsProvider::user_nuisances_path
     *
     * @return `DirPaths::user_dir_path / "parameters" / "nuisances.yaml"`.
     */
    fs::path user_nuisances_path() const override {
        return DirPaths::user_dir_path / "parameters" / "nuisances.yaml";
    }
};

#endif // DEFAULTNUISANCEPATHSPROVIDER_H