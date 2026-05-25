#ifndef INUISANCEPATHSPROVIDER_H
#define INUISANCEPATHSPROVIDER_H

#include <filesystem>

namespace fs = std::filesystem;

/**
 * @file INuisancePathsProvider.h
 * @brief Interface for locating nuisance-parameter configuration files.
 *
 * @see DefaultNuisancePathsProvider
 * @see NuisanceReader
 */

/**
 * @class INuisancePathsProvider
 * @brief Abstract provider for default and user nuisance-file paths.
 *
 * This interface allows tests, applications, and command-line frontends to
 * control where nuisance specifications are loaded from without hard-coding the
 * path logic inside @ref NuisanceReader.
 */
class INuisancePathsProvider {
public:
    /**
     * @brief Virtual destructor for polymorphic use.
     */
    virtual ~INuisancePathsProvider() = default;

    /**
     * @brief Returns the path to the built-in nuisance definition file.
     *
     * @return Filesystem path to the default nuisance registry.
     */
    virtual fs::path default_nuisances_path() const = 0;

    /**
     * @brief Returns the path to the user override nuisance definition file.
     *
     * @return Filesystem path to the user nuisance registry. The returned path
     *         may be empty when no user override file is configured.
     */
    virtual fs::path user_nuisances_path() const = 0;
};

#endif // INUISANCEPATHSPROVIDER_H