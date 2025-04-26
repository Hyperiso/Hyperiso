#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <map>
#include <optional>
#include <filesystem>

#include "General.h"

enum class ExternalFlag { IS_LHA_SPECTRUM, HAS_WILSON_INPUT, HAS_TH_OBSERVABLE_INPUT, USE_MARTY };

struct Config {
    std::map<ExternalFlag, bool> flags {
        {ExternalFlag::IS_LHA_SPECTRUM, false},
        {ExternalFlag::HAS_WILSON_INPUT, false},
        {ExternalFlag::HAS_TH_OBSERVABLE_INPUT, false},
        {ExternalFlag::USE_MARTY, false},
    };
    Model model {Model::SM};                    ///< Model type (current model)
    std::optional<std::string> mty_model_name;  ///< MARTY model name (name of the class in MARTY) if needed
    std::optional<fs::path> mty_model_path;     ///< Path to the MARTY model file (mty_model_name.h) if needed
};

#endif // __CONFIG_H__
