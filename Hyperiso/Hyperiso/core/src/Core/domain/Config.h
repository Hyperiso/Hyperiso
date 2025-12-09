#ifndef CONFIG_H
#define CONFIG_H

#include <map>
#include <optional>
#include <filesystem>

#include "Include.h"

/**
 * @brief Flags describing external input characteristics.
 */
enum class ExternalFlag { 
    IS_LHA_SPECTRUM,            ///< Input LHA file already contains a spectrum.
    HAS_WILSON_INPUT,           ///< User provided Wilson coefficient input.
    HAS_TH_OBSERVABLE_INPUT     ///< User provided theoretical observable input
};

/**
 * @struct Config
 * @brief Configuration object controlling model and input flags.
 */
struct Config {
    /// External flags describing the nature of the inputs.
    std::map<ExternalFlag, bool> flags {
        {ExternalFlag::IS_LHA_SPECTRUM, false},
        {ExternalFlag::HAS_WILSON_INPUT, false},
        {ExternalFlag::HAS_TH_OBSERVABLE_INPUT, false}
    };
    Model model {Model::SM};                    ///< Model type (current model)
    std::optional<std::string> mty_model_name;  ///< MARTY model name (name of the class in MARTY) if needed
    std::optional<fs::path> mty_model_path;     ///< Path to the MARTY model file (mty_model_name.h) if needed
};

#endif // CONFIG_H
