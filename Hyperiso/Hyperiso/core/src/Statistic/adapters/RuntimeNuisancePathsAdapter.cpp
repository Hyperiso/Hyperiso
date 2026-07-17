#include "RuntimeNuisancePathsAdapter.h"

#include <stdexcept>
#include <utility>

#include "APIAdapter.h"

RuntimeNuisancePathsAdapter::RuntimeNuisancePathsAdapter()
    : RuntimeNuisancePathsAdapter(std::make_shared<APIAdapter>())
{
}

RuntimeNuisancePathsAdapter::RuntimeNuisancePathsAdapter(
    std::shared_ptr<IPathProvider<APIPath>> path_provider
)
    : path_provider_(std::move(path_provider))
{
    if (!path_provider_) {
        throw std::invalid_argument(
            "RuntimeNuisancePathsAdapter requires a non-null path provider"
        );
    }
}

fs::path RuntimeNuisancePathsAdapter::default_nuisances_path() const
{
    return path_provider_->get_path(APIPath::DEFAULT_NUISANCES);
}

fs::path RuntimeNuisancePathsAdapter::user_nuisances_path() const
{
    return path_provider_->get_path(APIPath::USER_NUISANCES);
}
