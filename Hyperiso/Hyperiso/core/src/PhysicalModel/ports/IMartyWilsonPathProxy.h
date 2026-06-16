#ifndef IMARTY_WILSON_PATH_PROXY_H
#define IMARTY_WILSON_PATH_PROXY_H

#include <string>
#include "Include.h"

/**
 * @file IMartyWilsonPathProxy.h
 * @brief Port used by the Wilson layer to access MARTY-generated output paths.
 *
 * The Wilson domain must not know where Hyperiso stores writable runtime files.
 * Implementations bridge to the Core path system, typically through MartyAdapter
 * and MemoryManager.
 */
class IMartyWilsonPathProxy {
public:
    virtual ~IMartyWilsonPathProxy() = default;

    /**
     * @brief Writable directory containing MARTY generated files and CSV outputs.
     */
    virtual fs::path marty_temp_dir() const = 0;

    /**
     * @brief CSV output produced by MARTY for a given model name.
     */
    virtual fs::path wilson_csv_path(const std::string& model_name) const {
        return marty_temp_dir() / (model_name + "_wilson.csv");
    }
};

#endif // IMARTY_WILSON_PATH_PROVIDER_H
