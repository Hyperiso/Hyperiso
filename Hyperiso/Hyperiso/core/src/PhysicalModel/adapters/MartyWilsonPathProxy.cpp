#include "MartyWilsonPathProxy.h"

#include <stdexcept>
#include <system_error>

fs::path MartyWilsonPathProxy::marty_temp_dir() const {
    fs::path dir = MartyAdapter().get_path(MartyPath::MARTY_TEMP_DIR);

    std::error_code ec;
    fs::create_directories(dir, ec);
    if (ec) {
        throw std::runtime_error("Cannot create MARTY temp directory: " + dir.string() + " (" + ec.message() + ")");
    }

    return dir;
}

fs::path MartyWilsonPathProxy::wilson_csv_path(const std::string& model_name) const {
    return marty_temp_dir() / (model_name + "_wilson.csv");
}
