#include "SoftSusy.h"

#include <cstdlib>
#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

#include "Logger.h"
#include "config.hpp"

namespace fs = std::filesystem;

namespace {

std::optional<fs::path>& external_softsusy_path()
{
    static std::optional<fs::path> path;
    return path;
}

std::string shell_quote(const std::string& value)
{
    std::string quoted = "'";
    for (char c : value) {
        if (c == '\'') {
            quoted += "'\\''";
        } else {
            quoted += c;
        }
    }
    quoted += "'";
    return quoted;
}

bool is_existing_file(const fs::path& path)
{
    std::error_code ec;
    return fs::exists(path, ec) && fs::is_regular_file(path, ec);
}

fs::path normalize_path(const fs::path& path)
{
    std::error_code ec;
    const fs::path canonical = fs::weakly_canonical(path, ec);
    return ec ? fs::absolute(path) : canonical;
}

std::vector<fs::path> softsusy_candidates_from(const fs::path& raw)
{
    std::vector<fs::path> candidates;
    candidates.push_back(raw);

    std::error_code ec;
    if (fs::exists(raw, ec) && fs::is_directory(raw, ec)) {
        candidates.push_back(raw / "softpoint.x");
        candidates.push_back(raw / "bin" / "softpoint.x");
        candidates.push_back(raw / "src" / "SOFTSUSY" / "softpoint.x");
        candidates.push_back(raw / "SOFTSUSY" / "softpoint.x");
    }

    return candidates;
}

SoftsusyRuntimeConfig::Resolution resolve_from_path(const fs::path& raw)
{
    for (const auto& candidate : softsusy_candidates_from(raw)) {
        if (is_existing_file(candidate)) {
            SoftsusyRuntimeConfig::Resolution out;
            out.valid = true;
            out.executable = normalize_path(candidate);
            return out;
        }
    }

    SoftsusyRuntimeConfig::Resolution out;
    std::ostringstream oss;
    oss << "No SOFTSUSY executable found from path '" << raw.string() << "'. "
        << "Expected a softpoint.x executable, or a directory containing "
        << "softpoint.x, bin/softpoint.x, or src/SOFTSUSY/softpoint.x.";
    out.error = oss.str();
    return out;
}

std::optional<fs::path> environment_softsusy_path()
{
    if (const char* env = std::getenv("HYPERISO_SOFTSUSY")) {
        if (std::string(env).empty() == false) {
            return fs::path(env);
        }
    }
    return std::nullopt;
}

std::optional<fs::path> path_softsusy_executable()
{
    if (const char* env_path = std::getenv("PATH")) {
        std::string paths(env_path);
        std::size_t start = 0;
        while (start <= paths.size()) {
            const std::size_t end = paths.find(':', start);
            const std::string entry = paths.substr(start, end == std::string::npos ? std::string::npos : end - start);
            if (!entry.empty()) {
                fs::path candidate = fs::path(entry) / "softpoint.x";
                if (is_existing_file(candidate)) {
                    return candidate;
                }
            }
            if (end == std::string::npos) {
                break;
            }
            start = end + 1;
        }
    }
    return std::nullopt;
}

#ifdef BUILD_WITH_SOFTSUSY
std::optional<fs::path> bundled_softsusy_executable()
{
    const fs::path root_tp_file(project_tp_root.data());
    const std::vector<fs::path> candidates = {
        root_tp_file / "SOFTSUSY" / "src" / "SOFTSUSY" / "softpoint.x",
        root_tp_file / "SOFTSUSY" / "softpoint.x",
        root_tp_file / "SOFTSUSY" / "bin" / "softpoint.x"
    };

    for (const auto& candidate : candidates) {
        if (is_existing_file(candidate)) {
            return candidate;
        }
    }
    return std::nullopt;
}
#endif

fs::path resolve_input_path(const std::string& inputFilePath)
{
    const fs::path input(inputFilePath);
    if (input.is_absolute()) {
        return input;
    }
    return fs::path(project_assets_root.data()) / input;
}

} // namespace

namespace SoftsusyRuntimeConfig {

Resolution set_external_path(const std::string& path)
{
    const auto resolved = resolve_from_path(fs::path(path));
    if (resolved.valid) {
        external_softsusy_path() = resolved.executable;
    }
    return resolved;
}

Resolution resolve_executable()
{
    if (external_softsusy_path().has_value()) {
        return resolve_from_path(*external_softsusy_path());
    }

    if (auto env_path = environment_softsusy_path()) {
        auto resolved = resolve_from_path(*env_path);
        if (resolved.valid) {
            return resolved;
        }
    }

    if (auto path_candidate = path_softsusy_executable()) {
        auto resolved = resolve_from_path(*path_candidate);
        if (resolved.valid) {
            return resolved;
        }
    }

#ifdef BUILD_WITH_SOFTSUSY
    if (auto bundled = bundled_softsusy_executable()) {
        auto resolved = resolve_from_path(*bundled);
        if (resolved.valid) {
            return resolved;
        }
    }
#endif

    Resolution out;
    std::ostringstream oss;
    oss << "Cannot compute a SUSY spectrum because SOFTSUSY softpoint.x was not found. "
        << "Provide an executable before initialization with "
        << "HyperisoMaster::pre_init_set_softsusy_path('/path/to/softpoint.x') "
        << "or set HYPERISO_SOFTSUSY=/path/to/softpoint.x.";
#ifndef BUILD_WITH_SOFTSUSY
    oss << " This Hyperiso build was not configured with BUILD_WITH_SOFTSUSY=ON, "
        << "so no bundled SOFTSUSY fallback is available.";
#endif
    oss << " If your LHA file already contains a spectrum, set IS_LHA_SPECTRUM=True.";
    out.error = oss.str();
    return out;
}

std::string availability_message()
{
    const auto resolved = resolve_executable();
    if (resolved.valid) {
        return "SOFTSUSY executable: " + resolved.executable.string();
    }
    return resolved.error;
}

} // namespace SoftsusyRuntimeConfig

void SoftsusyCalculator::calculateSpectrum(const std::string& inputFilePath, const std::string& outputFilePath)
{
    const auto resolved = SoftsusyRuntimeConfig::resolve_executable();
    if (!resolved.valid) {
        throw std::runtime_error(resolved.error);
    }

    const fs::path input = resolve_input_path(inputFilePath);
    const fs::path output(outputFilePath);

    std::error_code ec;
    if (output.has_parent_path()) {
        fs::create_directories(output.parent_path(), ec);
        if (ec) {
            throw std::runtime_error("Cannot create SOFTSUSY output directory '" + output.parent_path().string() + "': " + ec.message());
        }
    }

    const std::string command = shell_quote(resolved.executable.string()) +
                                " leshouches < " + shell_quote(input.string()) +
                                " > " + shell_quote(output.string());

    LOG_DEBUG("SOFTSUSY COMMAND : " + command);

    const int result = std::system(command.c_str());
    if (result != 0) {
        throw std::runtime_error("SOFTSUSY execution failed with code " + std::to_string(result) +
                                 ". Command was: " + command);
    }

    LOG_INFO("SOFTSUSY execution successful.");
}
