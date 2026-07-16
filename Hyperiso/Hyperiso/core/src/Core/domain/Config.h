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
    HAS_TH_OBSERVABLE_INPUT,    ///< User provided theoretical observable input
    HYP_AS_SM_MARTY             ///< If true, use Hyperiso as SM values for Wilson coefficients (up to NNLO).
};

/**
 * @struct HyperisoConfig
 * @brief Configuration object controlling model, input flags and optional MARTY resources.
 */
struct HyperisoConfig {
    /// External flags describing the nature of the inputs.
    std::map<ExternalFlag, bool> flags {
        {ExternalFlag::IS_LHA_SPECTRUM, false},
        {ExternalFlag::HAS_WILSON_INPUT, false},
        {ExternalFlag::HAS_TH_OBSERVABLE_INPUT, false},
        {ExternalFlag::HYP_AS_SM_MARTY, false}
    };

    Model model {Model::SM};                    ///< Current model.
    std::optional<std::string> mty_model_name;  ///< MARTY model class name, if needed.
    std::optional<fs::path> mty_model_path;     ///< Path to the MARTY model file, if needed.

    /**
     * @brief Optional user-provided BSM mapping JSON.
     *
     * The read-only SM mapping is resolved through the path provider
     * (MartyPath::SM_MAPPING_FILE). This optional path is only for the BSM
     * mapping supplied by the user, e.g. a ZPrime mapping.
     *
     * The mapping values are intentionally strings, because identifiers may be
     * composite values such as "1_2".
     */
    std::optional<fs::path> mty_bsm_mapping_path;
};

#endif // CONFIG_H