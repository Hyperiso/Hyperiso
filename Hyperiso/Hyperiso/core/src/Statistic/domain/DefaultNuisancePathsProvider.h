#ifndef DEFAULTNUISANCEPATHSPROVIDER_H
#define DEFAULTNUISANCEPATHSPROVIDER_H

#include "INuisancePathsProvider.h"
#include "Paths.h"


struct DefaultNuisancePathsProvider : public INuisancePathsProvider {
    fs::path default_nuisances_path() const override {
        return DirPaths::default_dir_path / "nuisances.json";
    }

    fs::path user_nuisances_path() const override {
        return DirPaths::user_dir_path / "parameters" / "nuisances.yaml";
    }
};

#endif // DEFAULTNUISANCEPATHSPROVIDER_H