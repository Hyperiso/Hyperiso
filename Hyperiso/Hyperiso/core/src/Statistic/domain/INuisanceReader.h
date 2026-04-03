#ifndef INUISANCEREADER_H
#define INUISANCEREADER_H

#include <filesystem>
#include "NuisanceSpec.h"

namespace fs = std::filesystem;

class INuisanceReader {
public:
    virtual ~INuisanceReader() = default;

    virtual NuisanceRegistry load_default() const = 0;

    virtual NuisanceRegistry load_user() const = 0;

    virtual NuisanceRegistry load_user(const fs::path& user_path) const = 0;
};

#endif // INUISANCEREADER_H