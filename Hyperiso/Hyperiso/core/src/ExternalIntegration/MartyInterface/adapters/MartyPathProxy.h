#ifndef HYPERISO_MARTY_PATH_PROXY_H
#define HYPERISO_MARTY_PATH_PROXY_H

#include <memory>
#include <optional>

#include "IMartyPathAdapter.h"
#include "IMartyPathProxy.h"

/**
 * @brief MARTY-side path proxy.
 *
 * This class deliberately does not know MemoryManager. It delegates every path
 * request to an injected core-side adapter.
 */
class MartyPathProxy : public IMartyPathProxy {
public:
    explicit MartyPathProxy(std::shared_ptr<IMartyPathAdapter> adapter);

    fs::path get_path(MartyPath path_name) const override;
    std::optional<fs::path> get_optional_path(MartyPath path_name) const override;

private:
    std::shared_ptr<IMartyPathAdapter> adapter_;
};

#endif // HYPERISO_MARTY_PATH_PROXY_H
