#ifndef DEFAULTNUISANCEPATHSPROVIDER_H
#define DEFAULTNUISANCEPATHSPROVIDER_H

#include <memory>
#include <utility>

#include "INuisancePathsProvider.h"
#include "IPathsProvider.h"
#include "DefaultPathsProvider.h"

/**
 * @file DefaultNuisancePathsProvider.h
 * @brief Default nuisance-configuration path provider.
 *
 * Nuisance paths are now resolved through IPathsProvider so they follow the
 * same read-only assets root and pre-init override system as the rest of
 * Hyperiso paths.
 *
 * @see INuisancePathsProvider
 * @see NuisanceReader
 */
struct DefaultNuisancePathsProvider : public INuisancePathsProvider {
    explicit DefaultNuisancePathsProvider(
        std::shared_ptr<IPathsProvider> paths_provider = std::make_shared<DefaultPathsProvider>()
    )
        : paths_provider_(std::move(paths_provider)) {}

    /**
     * @copydoc INuisancePathsProvider::default_nuisances_path
     *
     * Default: ASSETS_ROOT / "default" / "nuisances.json".
     * Override through APIPath::DEFAULT_NUISANCES.
     */
    fs::path default_nuisances_path() const override {
        return paths_provider_->default_nuisances();
    }

    /**
     * @copydoc INuisancePathsProvider::user_nuisances_path
     *
     * Default: ASSETS_ROOT / "input_files" / "parameters" / "nuisances.yaml".
     * Override through APIPath::USER_NUISANCES.
     */
    fs::path user_nuisances_path() const override {
        return paths_provider_->user_nuisances();
    }

private:
    std::shared_ptr<IPathsProvider> paths_provider_;
};

#endif // DEFAULTNUISANCEPATHSPROVIDER_H
