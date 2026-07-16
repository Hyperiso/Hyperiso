#ifndef HYPERISO_MARTY_RUNTIME_CONFIG_H
#define HYPERISO_MARTY_RUNTIME_CONFIG_H

#include <algorithm>
#include <filesystem>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "Logger.h"
#include "config.hpp"

/**
 * @brief Runtime resolver for the MARTY installation used by generated code.
 *
 * The resolver lets API code register an existing MARTY installation before
 * Hyperiso initialization. If no explicit path is registered, it falls back to
 * the historical Third_party/MARTY/MARTY_INSTALL location produced by
 * -DBUILD_WITH_MARTY=ON.
 *
 * A valid installation prefix is expected to expose at least:
 *   - include/marty.h
 *   - lib/libmarty.so, lib/libmarty.dylib or lib/libmarty.a
 *
 * A bin/marty executable is detected and reported when present, but is not
 * required by the current generated-code pipeline, which links libmarty
 * directly.
 */
class MartyRuntimeConfig {
public:
    struct InstallInfo {
        std::filesystem::path requested_path;
        std::filesystem::path prefix;
        std::filesystem::path include_dir;
        std::filesystem::path lib_dir;
        std::filesystem::path marty_header;
        std::filesystem::path marty_library;
        std::filesystem::path marty_executable;
        std::string source;
        std::string error;
        bool has_executable = false;
        bool valid = false;
    };

    /**
     * @brief Register a user-provided MARTY installation path.
     *
     * The input can be either the install prefix itself, a path containing a
     * MARTY_INSTALL/ or install/ child, the include directory, the lib
     * directory, the marty.h file, or the libmarty library file.
     */
    static InstallInfo set_external_install_path(const std::filesystem::path& path)
    {
        external_install_path() = path;
        return inspect_install_path(path, "user-provided MARTY path");
    }

    /**
     * @brief Clear the user-provided path and return to the default fallback.
     */
    static void clear_external_install_path()
    {
        external_install_path().reset();
    }

    /**
     * @brief Resolve the MARTY install that should be used now.
     */
    static InstallInfo resolve()
    {
        if (external_install_path().has_value()) {
            return inspect_install_path(*external_install_path(), "user-provided MARTY path");
        }

        return inspect_install_path(default_install_path(), default_source_label());
    }

    /**
     * @brief Validate and return the active MARTY installation.
     *
     * Emits a detailed LOG_ERROR when neither a valid external installation nor
     * the default BUILD_WITH_MARTY installation is available.
     */
    static InstallInfo require_available(const std::string& context)
    {
        auto info = resolve();
        if (info.valid) {
            LOG_DEBUG("MARTY runtime resolved for", context, ":", info.prefix.string());
            return info;
        }

        LOG_ERROR(
            "MartyConfigError",
            "MARTY mode requested in", context,
            "but no valid MARTY installation could be resolved.",
            info.error,
            "Call HyperisoMaster::pre_init_set_marty_path(<MARTY_INSTALL>) before init(),",
            "or build/install the bundled MARTY with -DBUILD_WITH_MARTY=ON."
        );
        return info;
    }

    /**
     * @brief Shell-quote a path for command strings passed to std::system().
     */
    static std::string shell_quote(const std::filesystem::path& path)
    {
        std::string raw = path.string();
        std::string quoted = "'";
        for (char c : raw) {
            if (c == '\'') {
                quoted += "'\\''";
            } else {
                quoted += c;
            }
        }
        quoted += "'";
        return quoted;
    }

    static std::filesystem::path default_install_path()
    {
        return std::filesystem::path(project_tp_root.data()) / "MARTY" / "MARTY_INSTALL";
    }

    static bool was_external_path_provided()
    {
        return external_install_path().has_value();
    }

private:
    static std::optional<std::filesystem::path>& external_install_path()
    {
        static std::optional<std::filesystem::path> path;
        return path;
    }

    static std::string default_source_label()
    {
#ifdef BUILD_WITH_MARTY
        return "bundled MARTY install from BUILD_WITH_MARTY";
#else
        return "default bundled MARTY install path";
#endif
    }

    static bool is_marty_library_name(const std::filesystem::path& path)
    {
        const auto filename = path.filename().string();
        return filename == "libmarty.so"
            || filename == "libmarty.dylib"
            || filename == "libmarty.a";
    }

    static std::optional<std::filesystem::path> find_marty_library(const std::filesystem::path& lib_dir)
    {
        const std::vector<std::string> names = {
            "libmarty.so",
            "libmarty.dylib",
            "libmarty.a"
        };

        for (const auto& name : names) {
            std::filesystem::path candidate = lib_dir / name;
            if (std::filesystem::exists(candidate)) {
                return candidate;
            }
        }
        return std::nullopt;
    }

    static std::filesystem::path normalize_prefix_candidate(const std::filesystem::path& requested)
    {
        if (requested.empty()) {
            return requested;
        }

        std::filesystem::path p = std::filesystem::absolute(requested).lexically_normal();

        if (std::filesystem::is_regular_file(p)) {
            if (p.filename() == "marty.h") {
                return p.parent_path().parent_path();
            }
            if (is_marty_library_name(p)) {
                return p.parent_path().parent_path();
            }
        }

        if (p.filename() == "include" && std::filesystem::exists(p / "marty.h")) {
            return p.parent_path();
        }

        if (p.filename() == "lib" && find_marty_library(p).has_value()) {
            return p.parent_path();
        }

        if (std::filesystem::exists(p / "include" / "marty.h")) {
            return p;
        }

        if (std::filesystem::exists(p / "MARTY_INSTALL" / "include" / "marty.h")) {
            return p / "MARTY_INSTALL";
        }

        if (std::filesystem::exists(p / "install" / "include" / "marty.h")) {
            return p / "install";
        }

        return p;
    }

    static InstallInfo inspect_install_path(const std::filesystem::path& requested, std::string source)
    {
        InstallInfo info;
        info.requested_path = requested;
        info.source = std::move(source);
        info.prefix = normalize_prefix_candidate(requested);
        info.include_dir = info.prefix / "include";
        info.lib_dir = info.prefix / "lib";
        info.marty_header = info.include_dir / "marty.h";
        info.marty_executable = info.prefix / "bin" / "marty";

        std::vector<std::string> errors;

        if (requested.empty()) {
            errors.emplace_back("empty path");
        }

        if (!std::filesystem::exists(info.prefix)) {
            errors.emplace_back("prefix does not exist: " + info.prefix.string());
        }

        if (!std::filesystem::exists(info.marty_header)) {
            errors.emplace_back("missing header: " + info.marty_header.string());
        }

        auto lib = find_marty_library(info.lib_dir);
        if (!lib.has_value()) {
            errors.emplace_back("missing libmarty in: " + info.lib_dir.string());
        } else {
            info.marty_library = *lib;
        }

        info.has_executable = std::filesystem::exists(info.marty_executable)
                           && !std::filesystem::is_directory(info.marty_executable);

        info.valid = errors.empty();
        if (!info.valid) {
            std::ostringstream oss;
            oss << "Checked " << info.source << " at " << requested.string() << ". ";
            for (size_t i = 0; i < errors.size(); ++i) {
                if (i > 0) {
                    oss << " ";
                }
                oss << errors[i] << ".";
            }
            info.error = oss.str();
        }

        return info;
    }
};

#endif // HYPERISO_MARTY_RUNTIME_CONFIG_H
