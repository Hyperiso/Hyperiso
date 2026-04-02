#ifndef INUISANCEPATHSPROVIDER_H
#define INUISANCEPATHSPROVIDER_H

#include <filesystem>

namespace fs = std::filesystem;


class INuisancePathsProvider {
public:
    virtual ~INuisancePathsProvider() = default;

    virtual fs::path default_nuisances_path() const = 0;

    virtual fs::path user_nuisances_path() const = 0;
};

#endif // INUISANCEPATHSPROVIDER_H