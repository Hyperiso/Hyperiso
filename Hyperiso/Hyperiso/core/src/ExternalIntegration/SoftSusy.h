#ifndef SOFT_SUSY_H
#define SOFT_SUSY_H

#include <filesystem>
#include <optional>
#include <string>

#include "Interface.h"

namespace SoftsusyRuntimeConfig {

struct Resolution {
    bool valid = false;
    std::filesystem::path executable;
    std::string error;
};

/**
 * Register a SOFTSUSY executable or installation directory at runtime.
 *
 * The path may point directly to softpoint.x, or to a directory containing one
 * of the usual layouts: softpoint.x, bin/softpoint.x, or
 * src/SOFTSUSY/softpoint.x.
 */
Resolution set_external_path(const std::string& path);

/** Return the currently configured/resolved SOFTSUSY executable, if any. */
Resolution resolve_executable();

/** Return a human-readable diagnostic for the current SOFTSUSY setup. */
std::string availability_message();

} // namespace SoftsusyRuntimeConfig

/**
 * Concrete implementation of ICalculator for SOFTSUSY.
 *
 * The class is always compiled so Python users can provide an external
 * softpoint.x at runtime. If no runtime path is provided, it can fall back to
 * the bundled Third_party/SOFTSUSY tree only when Hyperiso was configured with
 * BUILD_WITH_SOFTSUSY=ON.
 */
class SoftsusyCalculator : public ICalculator {
public:
    void calculateSpectrum(const std::string& inputFilePath, const std::string& outputFilePath) override;
};

#endif // SOFT_SUSY_H
