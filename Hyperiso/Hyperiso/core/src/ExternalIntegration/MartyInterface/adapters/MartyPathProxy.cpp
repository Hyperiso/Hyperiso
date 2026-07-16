#include "MartyPathProxy.h"

#include <stdexcept>
#include <utility>

MartyPathProxy::MartyPathProxy(std::shared_ptr<IMartyPathAdapter> adapter)
    : adapter_(std::move(adapter)) {
    if (!adapter_) {
        throw std::invalid_argument("MartyPathProxy requires a non-null IMartyPathAdapter.");
    }
}

fs::path MartyPathProxy::get_path(MartyPath path_name) const {
    return adapter_->get_path(path_name);
}

std::optional<fs::path> MartyPathProxy::get_optional_path(MartyPath path_name) const {
    return adapter_->get_optional_path(path_name);
}
