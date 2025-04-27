#ifndef IPATHPROVIDER_H
#define IPATHPROVIDER_H

#include <filesystem>

/**
 * @example examples/path_provider_example.cpp
 * @brief Example showing how to use classes that implement IPathProvider.
 * @defgroup PathProvidingModule Path Providing System
 * @brief Defines standardized interfaces for providing filesystem paths associated with the framework.
 *
 * ## Related Classes
 * - @ref IPathProvider
 * - @ref APIAdapter
 * - @ref MartyAdapter
 */

namespace fs = std::filesystem;

/**
 * @class IPathProvider
 * @ingroup PathProvidingModule
 * @brief Interface for providing specific filesystem paths based on enumerated identifiers.
 *
 * Template interface to retrieve filesystem paths dynamically at runtime depending on an Enum type.
 */
template <typename EnumType>
class IPathProvider {
public:
    virtual ~IPathProvider() = default;

    /**
     * @brief Retrieves a filesystem path corresponding to a specific enum identifier.
     * @param enum_val Enum value indicating the desired path.
     * @return Filesystem path corresponding to the enum.
     */
    virtual fs::path get_path(EnumType) = 0;
};


#endif // __IPATHPROVIDER_H__
