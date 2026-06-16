#ifndef HYPERISO_IMARTY_PATH_PROXY_H
#define HYPERISO_IMARTY_PATH_PROXY_H

#include <filesystem>
#include <optional>

#include "MartyPathTypes.h"

namespace fs = std::filesystem;

/**
 * @brief MARTY-interface-side proxy port.
 *
 * Consumers in the MARTY interface should use this interface instead of
 * MemoryManager or DefaultPathsProvider.
 */
struct IMartyPathProxy {
    virtual ~IMartyPathProxy() = default;

    virtual fs::path get_path(MartyPath path_name) const = 0;
    virtual std::optional<fs::path> get_optional_path(MartyPath path_name) const = 0;
};

#endif // HYPERISO_IMARTY_PATH_PROXY_H
