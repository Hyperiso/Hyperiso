#ifndef IPATHPROVIDER_H
#define IPATHPROVIDER_H

#include <filesystem>

namespace fs = std::filesystem;

/**
 * @file IPathProvider.h
 * @brief Interface for providing concrete filesystem paths based on enum identifiers.
 *
 * This header defines the templated IPathProvider interface, used by:
 *  - @ref APIAdapter   (for high-level API paths, see APIPath)
 *  - @ref MartyAdapter (for MARTY-specific paths, see MartyPath)
 */

/**
 * @example examples/path_provider_example.cpp
 * @brief Example showing how to use classes that implement IPathProvider.
 *
 * @defgroup PathProvidingModule Path Providing System
 * @brief Defines standardized interfaces for providing filesystem paths associated with the framework.
 *
 * ## Related Classes
 * - @ref IPathProvider
 * - @ref APIAdapter
 * - @ref MartyAdapter
 */

/**
 * @class IPathProvider
 * @ingroup PathProvidingModule
 * @brief Interface for providing specific filesystem paths based on enumerated identifiers.
 *
 * @tparam EnumType Enum describing the set of available paths for a given provider
 *         (e.g. APIPath, MartyPath).
 *
 * Implementations map EnumType values to concrete filesystem locations,
 * typically derived from configuration or a central paths provider.
 */
template <typename EnumType>
class IPathProvider {
public:
    /// Virtual destructor for polymorphic use.
    virtual ~IPathProvider() = default;

    /**
     * @brief Retrieves a filesystem path corresponding to a specific enum identifier.
     *
     * @param enum_val Enum value indicating the desired path.
     * @return Filesystem path corresponding to the enum.
     */
    virtual fs::path get_path(EnumType) = 0;
};


#endif // IPATHPROVIDER_H
