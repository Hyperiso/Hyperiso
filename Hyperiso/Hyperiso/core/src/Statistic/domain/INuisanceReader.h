#ifndef INUISANCEREADER_H
#define INUISANCEREADER_H

#include <filesystem>

#include "NuisanceSpec.h"

namespace fs = std::filesystem;

/**
 * @file INuisanceReader.h
 * @brief Interface for loading nuisance-parameter specifications.
 *
 * The reader abstraction separates the source of nuisance definitions from the
 * code that consumes them. Implementations may load built-in defaults, user
 * overrides, or explicit user-provided files.
 *
 * @see NuisanceReader
 * @see NuisanceSpec
 */

/**
 * @class INuisanceReader
 * @brief Abstract interface for nuisance specification loaders.
 *
 * A nuisance reader returns registries keyed by @ref ParamId. The default and
 * user registries can be loaded independently so that callers may decide how to
 * merge them.
 */
class INuisanceReader {
public:
    /**
     * @brief Virtual destructor for polymorphic use.
     */
    virtual ~INuisanceReader() = default;

    /**
     * @brief Loads the built-in nuisance-parameter registry.
     *
     * @return Registry containing the default nuisance specifications.
     *
     * @throws std::runtime_error if the default source is missing, malformed,
     *         or cannot be parsed by the concrete implementation.
     */
    virtual NuisanceRegistry load_default() const = 0;

    /**
     * @brief Loads the nuisance registry from the configured user source.
     *
     * @return Registry containing user-defined nuisance specifications. The
     *         registry may be empty when no user source is configured.
     *
     * @throws std::runtime_error if the configured user source exists but cannot
     *         be read or parsed.
     */
    virtual NuisanceRegistry load_user() const = 0;

    /**
     * @brief Loads the nuisance registry from an explicit user file.
     *
     * @param user_path Path to the user nuisance file. An empty path may be
     *        interpreted by implementations as "no user overrides".
     *
     * @return Registry containing nuisance specifications from @p user_path.
     *
     * @throws std::runtime_error if @p user_path is non-empty but missing,
     *         unreadable, or malformed.
     */
    virtual NuisanceRegistry load_user(const fs::path& user_path) const = 0;
};

#endif // INUISANCEREADER_H