#ifndef MARTY_WILSON_PATH_ADAPTER_H
#define MARTY_WILSON_PATH_ADAPTER_H

#include "IMartyWilsonPathProxy.h"
#include "MartyAdapter.h"

/**
 * @file MartyWilsonPathProxy.h
 * @brief Core-side adapter exposing MemoryManager-managed MARTY paths to Wilson.
 */
class MartyWilsonPathProxy : public IMartyWilsonPathProxy {
public:
    fs::path marty_temp_dir() const override;
    fs::path wilson_csv_path(const std::string& model_name) const override;
};

#endif // MARTY_WILSON_PATH_ADAPTER_H
