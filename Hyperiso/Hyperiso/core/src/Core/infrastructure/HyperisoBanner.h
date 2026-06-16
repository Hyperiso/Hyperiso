#ifndef HYPERISO_BANNER_H
#define HYPERISO_BANNER_H

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "Config.h"
#include "IPathsProvider.h"
#include "MemoryManager.h"

#ifndef HYPERISO_VERSION
#define HYPERISO_VERSION "0.1.0"
#endif

#ifndef HYPERISO_AUTHORS
#define HYPERISO_AUTHORS "Theo Reymermier and the Hyperiso contributors"
#endif

/**
 * @brief Pretty startup banner for Hyperiso sessions.
 *
 * Header-only by design: HyperisoMaster can include and call it without adding
 * another translation unit. The banner reads paths from the active
 * MemoryManager/IPathsProvider, so pre-init overrides are reflected.
 */
class HyperisoBanner {
public:
    static void print_startup(const std::string& lha_file, const HyperisoConfig& config)
    {
        auto* mm = MemoryManager::GetInstance();

        const std::string assets_root = get_path(mm, APIPath::ASSETS_ROOT);
        const std::string lha_path = resolve_lha_path(lha_file, assets_root);

        std::cout << "\n";
        print_rule('=');
        std::cout << "  H Y P E R I S O   v" << HYPERISO_VERSION << "\n";
        std::cout << "  High-energy flavour observables and Wilson coefficient engine\n";
        std::cout << "  Authors: " << HYPERISO_AUTHORS << "\n";
        print_rule('-');

        print_section("Session");
        print_kv("Model", model_name(config.model));
        print_kv("LHA input", lha_path);
        print_kv("Assets root", assets_root);

        print_section("Read-only defaults");
        print_kv("parameters", get_path(mm, APIPath::DEFAULT_PARAM_VALUES));
        print_kv("observables", get_path(mm, APIPath::DEFAULT_OBS_VALUES));
        print_kv("parameter correlations", get_path(mm, APIPath::DEFAULT_PARAM_CORR));
        print_kv("observable correlations", get_path(mm, APIPath::DEFAULT_OBS_CORR));
        print_kv("nuisances", get_path(mm, APIPath::DEFAULT_NUISANCES));

        print_section("Input overrides / packaged templates");
        print_kv("SM parameters", get_path(mm, APIPath::USER_SM_PARAMS));
        print_kv("flavour parameters", get_path(mm, APIPath::USER_FLAVOR_PARAMS));
        print_kv("decay parameters", get_path(mm, APIPath::USER_DECAY_PARAMS));
        print_kv("observables", get_path(mm, APIPath::USER_OBS_VALUES));
        print_kv("parameter correlations", get_path(mm, APIPath::USER_PARAM_CORR));
        print_kv("observable correlations", get_path(mm, APIPath::USER_OBS_CORR));
        print_kv("nuisances", get_path(mm, APIPath::USER_NUISANCES));

        print_section("Writable runtime directories");
        print_kv("spectrum cache", get_path(mm, APIPath::SPECTRUM_DIR));
        print_kv("MARTY temp/cache", get_path(mm, APIPath::MARTY_TEMP_DIR));

        if (uses_marty(config)) {
            print_section("MARTY backend");
            print_kv("model name", optional_string(config.mty_model_name));
            print_kv("model file", optional_path(config.mty_model_path));
            print_kv("template dir", get_path(mm, APIPath::TEMPLATE_DIR) + "/MARTY");
            print_kv("mapping dir", get_path(mm, APIPath::PARAM_MAPPING_DIR));
            print_kv("SM mapping", path_join(get_path(mm, APIPath::PARAM_MAPPING_DIR), "sm.json"));
            print_kv("BSM mapping", optional_path(config.mty_bsm_mapping_path));
            print_kv("Wilson CSV", marty_wilson_csv_path(mm, config));
        }

        print_rule('=');
        std::cout << std::endl;
    }

private:
    static constexpr int key_width = 24;
    static constexpr int banner_width = 92;

    static void print_rule(char c)
    {
        std::cout << std::string(banner_width, c) << "\n";
    }

    static void print_section(const std::string& name)
    {
        std::cout << "\n[" << name << "]\n";
    }

    static void print_kv(const std::string& key, const std::string& value)
    {
        std::cout << "  " << std::left << std::setw(key_width) << key << " : " << value << "\n";
    }

    static std::string get_path(MemoryManager* mm, APIPath path)
    {
        try {
            return mm->get_path(path).string();
        } catch (const std::exception& e) {
            return std::string("<unavailable: ") + e.what() + ">";
        } catch (...) {
            return "<unavailable>";
        }
    }

    static std::string path_join(const std::string& lhs, const std::string& rhs)
    {
        if (lhs.empty()) {
            return rhs;
        }
        return (std::filesystem::path(lhs) / rhs).string();
    }

    static std::string resolve_lha_path(const std::string& lha_file, const std::string& assets_root)
    {
        const std::filesystem::path input(lha_file);
        if (input.is_absolute()) {
            return input.string();
        }
        return (std::filesystem::path(assets_root) / input).string();
    }

    static std::string optional_string(const std::optional<std::string>& value)
    {
        return value.has_value() ? value.value() : std::string("<not set>");
    }

    static std::string optional_path(const std::optional<std::filesystem::path>& value)
    {
        return value.has_value() ? value.value().string() : std::string("<not set>");
    }

    static bool uses_marty(const HyperisoConfig& config)
    {
        return config.model == Model::MARTY
            || config.mty_model_name.has_value()
            || config.mty_model_path.has_value()
            || config.mty_bsm_mapping_path.has_value();
    }

    static std::string model_name(Model model)
    {
        if (model == Model::SM) return "SM";
        if (model == Model::THDM) return "THDM";
        if (model == Model::SUSY) return "SUSY";
        if (model == Model::MARTY) return "MARTY";

        std::ostringstream oss;
        oss << "Model(" << static_cast<int>(model) << ")";
        return oss.str();
    }

    static std::string marty_wilson_csv_path(MemoryManager* mm, const HyperisoConfig& config)
    {
        const std::string model = config.mty_model_name.value_or("SM");
        return (std::filesystem::path(get_path(mm, APIPath::MARTY_TEMP_DIR)) / (model + "_wilson.csv")).string();
    }
};

#endif // HYPERISO_BANNER_H
