#include <cassert>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>

#include "IPathProvider.h"
#include "IPathsProvider.h"
#include "RuntimeNuisancePathsAdapter.h"

namespace fs = std::filesystem;

class FakeRuntimePathProvider final : public IPathProvider<APIPath> {
public:
    fs::path get_path(APIPath path) override
    {
        switch (path) {
        case APIPath::DEFAULT_NUISANCES:
            return "/runtime/assets/default/nuisances.json";
        case APIPath::USER_NUISANCES:
            return "/runtime/assets/input_files/parameters/nuisances.yaml";
        default:
            return {};
        }
    }
};

int main()
{
    std::cout << "== Runtime nuisance paths adapter UNIT ==\n";

    auto core_paths = std::make_shared<FakeRuntimePathProvider>();
    RuntimeNuisancePathsAdapter adapter(core_paths);

    assert(adapter.default_nuisances_path() == fs::path("/runtime/assets/default/nuisances.json"));
    assert(adapter.user_nuisances_path() == fs::path("/runtime/assets/input_files/parameters/nuisances.yaml"));

    bool rejected_null = false;
    try {
        RuntimeNuisancePathsAdapter invalid(nullptr);
    } catch (const std::invalid_argument&) {
        rejected_null = true;
    }
    assert(rejected_null);

    std::cout << " UNIT OK\n";
    return 0;
}
