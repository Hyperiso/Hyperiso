#ifndef RUNTIME_NUISANCE_PATHS_ADAPTER_H
#define RUNTIME_NUISANCE_PATHS_ADAPTER_H

#include <memory>

#include "INuisancePathsProvider.h"
#include "IPathProvider.h"
#include "IPathsProvider.h"

/**
 * @file RuntimeNuisancePathsAdapter.h
 * @brief Bridges active Core runtime paths to Statistic nuisance-file lookup.
 */

/**
 * @class RuntimeNuisancePathsAdapter
 * @brief Implements the Statistic nuisance-path port using a Core path port.
 *
 * The production constructor uses APIAdapter, which delegates path lookup to
 * the IPathsProvider currently installed in MemoryManager. This preserves
 * packaged-wheel asset roots and pre-initialization path overrides without
 * exposing MemoryManager to the Statistic domain.
 */
class RuntimeNuisancePathsAdapter final : public INuisancePathsProvider {
public:
    /** @brief Construct an adapter backed by the active Core runtime. */
    RuntimeNuisancePathsAdapter();

    /**
     * @brief Construct an adapter with an injected Core path provider.
     * @param path_provider Provider implementing the Core APIPath port.
     * @throws std::invalid_argument when @p path_provider is null.
     */
    explicit RuntimeNuisancePathsAdapter(
        std::shared_ptr<IPathProvider<APIPath>> path_provider
    );

    /** @copydoc INuisancePathsProvider::default_nuisances_path */
    fs::path default_nuisances_path() const override;

    /** @copydoc INuisancePathsProvider::user_nuisances_path */
    fs::path user_nuisances_path() const override;

private:
    std::shared_ptr<IPathProvider<APIPath>> path_provider_;
};

#endif // RUNTIME_NUISANCE_PATHS_ADAPTER_H
