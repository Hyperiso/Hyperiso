#include "DefaultNuisancePathsProvider.h"
#include "NuisanceReader.h"

int main() {
    auto paths = std::make_shared<DefaultNuisancePathsProvider>();
    NuisanceReader reader(paths);

    // optionnel : si tu veux forcer un autre YAML user
    // reader.set_user_path("/tmp/my_nuisances.yaml");

    NuisanceRegistry registry = reader.load();
    // std::unordered_set<ParamId> ids = reader.load_param_ids();

    std::cout << registry << std::endl;
    
    // for (auto elem : ids) {
    //     std::cout << elem << std::endl;
    // }
    return 0;
}