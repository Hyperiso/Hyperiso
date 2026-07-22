#ifndef MARTYADAPTER_H
#define MARTYADAPTER_H

#include <optional>
#include <string>

#include "IMonitor.h"
#include "IPathProvider.h"
#include "IMartyPathAdapter.h"
#include "MemoryManager.h"

/**
 * @file MartyAdapter.h
 * @brief Core-side adapter giving MARTY access to configuration and managed paths.
 *
 * MartyAdapter is allowed to depend on MemoryManager. MARTY-interface code should
 * not depend on MemoryManager directly; it should use MartyPathProxy instead.
 */
class MartyAdapter : public IMonitor<InternalFlag>, public IPathProvider<MartyPath>, public IMartyPathAdapter {
public:
    /**
     * @brief Retrieves a required MARTY-related path.
     */
    fs::path get_path(MartyPath path_name) override;

    /**
     * @brief Retrieves an optional MARTY-related path.
     *
     * Currently used for the user-provided BSM mapping JSON.
     */
    std::optional<fs::path> get_optional_path(MartyPath path_name) override;

    /**
     * @brief Checks the status of a specific internal flag.
     */
    bool check_flag(InternalFlag flag) override;

    /**
     * @brief Retrieves the configured MARTY model name.
     */
    std::string get_marty_model_name() const;

    /** @brief Retrieves the configured MARTY perturbative-order policy. */
    MartyOrderPolicy get_marty_order_policy() const;
};

#endif // MARTYADAPTER_H
