#ifndef MARTYADAPTER_H
#define MARTYADAPTER_H

#include <map>
#include <optional>
#include <string>
#include <vector>

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

    /** @brief Retrieves per-coefficient MARTY TreeLevel fermion-order overrides. */
    std::map<std::string, std::vector<int>> get_marty_tree_fermion_orders() const;

    /** @brief Retrieves per-coefficient MARTY OneLoop fermion-order overrides. */
    std::map<std::string, std::vector<int>> get_marty_one_loop_fermion_orders() const;
};

#endif // MARTYADAPTER_H
