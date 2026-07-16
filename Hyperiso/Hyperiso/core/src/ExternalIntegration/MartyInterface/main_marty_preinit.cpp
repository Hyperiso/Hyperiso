#include "HyperisoMaster.h"
#include "MartyInterface.h"
#include "MartyRuntimeConfig.h"

#include <filesystem>
#include <iostream>
#include <optional>
#include <string>
#include <unordered_map>

namespace fs = std::filesystem;

namespace {

struct CliOptions {
    std::optional<fs::path> marty_path;
    fs::path lha_file = fs::path(project_root.data()) / "Test" / "InputFiles" / "testinput_thdm.lha";
    std::string wilson = "C7";
    std::string model = "SM";
    double q_match = 81.0;
    std::optional<fs::path> model_header;
    bool validate_only = false;
};

void print_usage(const char* argv0)
{
    std::cout
        << "Usage: " << argv0 << " [options]\n\n"
        << "Options:\n"
        << "  --marty <path>          Existing MARTY install prefix, include/, lib/, marty.h, or libmarty.*\n"
        << "  --lha <path>            LHA input file used to initialize Hyperiso\n"
        << "  --model-header <path>   MARTY model header used by generated code\n"
        << "  --wilson <name>         Wilson coefficient/basis to generate, default: C7\n"
        << "  --model <name>          Model name used in generated code, default: SM\n"
        << "  --q <value>             Matching scale, default: 81\n"
        << "  --validate-only         Resolve/validate MARTY and init Hyperiso, then stop\n"
        << "  -h, --help              Show this help\n\n"
        << "Examples:\n"
        << "  " << argv0 << " --marty /opt/MARTY_INSTALL --validate-only\n"
        << "  " << argv0 << " --marty /opt/MARTY_INSTALL --wilson C7 --model SM --q 81\n";
}

fs::path default_model_header_from(const MartyRuntimeConfig::InstallInfo& marty)
{
    const fs::path installed_header = marty.include_dir / "marty" / "models" / "sm.h";
    if (fs::exists(installed_header)) {
        return installed_header;
    }

    const fs::path bundled_source_header =
        fs::path(project_tp_root.data()) / "MARTY" / "src" / "MARTY" / "src" / "marty" / "models" / "sm.h";
    if (fs::exists(bundled_source_header)) {
        return bundled_source_header;
    }

    // Last-resort fallback: keep a deterministic path and let MartyInterface's
    // normal file checks report the precise problem if the user runs generation.
    return installed_header;
}

CliOptions parse_args(int argc, char** argv)
{
    CliOptions opts;

    const std::unordered_map<std::string, std::string> aliases = {
        {"--marty-path", "--marty"},
        {"--marty-install", "--marty"},
        {"--input", "--lha"},
        {"--model-path", "--model-header"},
        {"--Q", "--q"}
    };

    for (int i = 1; i < argc; ++i) {
        std::string key = argv[i];
        if (auto alias = aliases.find(key); alias != aliases.end()) {
            key = alias->second;
        }

        auto require_value = [&](const std::string& option) -> std::string {
            if (i + 1 >= argc) {
                throw std::runtime_error("Missing value after " + option);
            }
            return argv[++i];
        };

        if (key == "-h" || key == "--help") {
            print_usage(argv[0]);
            std::exit(0);
        } else if (key == "--marty") {
            opts.marty_path = fs::path(require_value(key));
        } else if (key == "--lha") {
            opts.lha_file = fs::path(require_value(key));
        } else if (key == "--model-header") {
            opts.model_header = fs::path(require_value(key));
        } else if (key == "--wilson") {
            opts.wilson = require_value(key);
        } else if (key == "--model") {
            opts.model = require_value(key);
        } else if (key == "--q") {
            opts.q_match = std::stod(require_value(key));
        } else if (key == "--validate-only") {
            opts.validate_only = true;
        } else {
            throw std::runtime_error("Unknown option: " + key);
        }
    }

    return opts;
}

void print_marty_resolution(const MartyRuntimeConfig::InstallInfo& info)
{
    std::cout << "MARTY resolution:\n"
              << "  source    : " << info.source << "\n"
              << "  requested : " << info.requested_path.string() << "\n"
              << "  prefix    : " << info.prefix.string() << "\n"
              << "  include   : " << info.include_dir.string() << "\n"
              << "  library   : " << info.marty_library.string() << "\n";

    if (info.has_executable) {
        std::cout << "  executable: " << info.marty_executable.string() << "\n";
    } else {
        std::cout << "  executable: not found, not required by current generated-code pipeline\n";
    }
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const CliOptions opts = parse_args(argc, argv);

        HyperisoMaster hyp;

        if (opts.marty_path.has_value()) {
            // This is the new pre-init API being tested. It validates and stores
            // the external MARTY installation before Hyperiso initialization.
            hyp.pre_init_set_marty_path(opts.marty_path->string());
        }

        auto marty = MartyRuntimeConfig::resolve();
        if (!marty.valid) {
            std::cerr << "Invalid MARTY runtime: " << marty.error << "\n";
            std::cerr << "Pass --marty <MARTY_INSTALL> or build with -DBUILD_WITH_MARTY=ON.\n";
            return 2;
        }
        print_marty_resolution(marty);

        HyperisoConfig config;
        config.model = Model::SM;
        config.flags[ExternalFlag::IS_LHA_SPECTRUM] = true;
        config.flags[ExternalFlag::HAS_WILSON_INPUT] = false;
        config.flags[ExternalFlag::HAS_TH_OBSERVABLE_INPUT] = false;
        config.flags[ExternalFlag::HYP_AS_SM_MARTY] = true;

        std::cout << "Initializing Hyperiso with LHA: " << opts.lha_file.string() << "\n";
        hyp.init(opts.lha_file.string(), config);

        if (opts.validate_only) {
            std::cout << "Validation finished successfully.\n";
            return 0;
        }

        const fs::path model_header = opts.model_header.value_or(default_model_header_from(marty));
        if (!fs::exists(model_header)) {
            std::cerr << "Model header not found: " << model_header.string() << "\n";
            std::cerr << "Pass --model-header <path/to/model.h>.\n";
            return 3;
        }

        std::cout << "Running MartyInterface.calculate(" << opts.wilson
                  << ", " << opts.model
                  << ", Q=" << opts.q_match
                  << ")\n";
        std::cout << "Model header: " << model_header.string() << "\n";

        MartyInterface marty_interface;
        marty_interface.calculate(
            opts.wilson,
            opts.model,
            opts.q_match,
            model_header.string()
        );

        std::cout << "MARTY pipeline finished.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << "\n";
        return 1;
    }
}
