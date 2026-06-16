#ifndef HYPERISO_IMARTY_PATH_ADAPTER_H
#define HYPERISO_IMARTY_PATH_ADAPTER_H

#include <filesystem>
#include <optional>

#include "MartyPathTypes.h"

namespace fs = std::filesystem;

/**
 * @brief Core-side port implemented by MartyAdapter.
 *
 * The MARTY interface layer talks to IMartyPathAdapter through a proxy. The
 * proxy never calls MemoryManager directly.
 */
struct IMartyPathAdapter {
    virtual ~IMartyPathAdapter() = default;

    virtual fs::path get_path(MartyPath path_name) = 0;
    virtual std::optional<fs::path> get_optional_path(MartyPath path_name) = 0;
};

#endif // HYPERISO_IMARTY_PATH_ADAPTER_H
