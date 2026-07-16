#include <iostream>

#include "HyperisoMaster.h"
#include "Include.h"
#include "Logger.h"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    // Configuration for HyperisoMaster, with global flags and model's informations.
    HyperisoConfig config;
    config.flags[ExternalFlag::IS_LHA_SPECTRUM] = false;        // If true, all spectrum parameters are read from the LHA. Useful mostly with THDM or SUSY.
    config.flags[ExternalFlag::HAS_WILSON_INPUT] = false;       // If true, Hyperiso searches Wilson coefficients in the FLHA input.
    config.flags[ExternalFlag::HAS_TH_OBSERVABLE_INPUT] = false;// If true, Hyperiso searches theoretical observable values in the FLHA input.
    config.flags[ExternalFlag::HYP_AS_SM_MARTY] = false;         // If true, Hyperiso computes SM Wilsons and lets MARTY provide only BSM contributions.

    config.model = Model::SM;                                   // Model: SM, THDM, SUSY or MARTY.
    config.mty_model_name = "marty_model_name";                 // Only used in Model::MARTY mode.
    config.mty_model_path = "/my/custom/marty/path";            // Only used in Model::MARTY mode.

    // Creation of the HyperisoMaster, base of all other interfaces and singleton.
    HyperisoMaster hyp;

    // Path to the LHA file for inputs. If relative, it starts from the Assets/ directory.
    const std::string lha_file_path = "lha/si_input.flha";

    // Initialization of Hyperiso, with configuration and path to an LHA for inputs.
    hyp.init(lha_file_path, config);

    // Some flag checks.
    std::cout << "Current model enum value: " << static_cast<int>(hyp.get_model()) << "\n";
    std::cout << "Flag IS_LHA_SPECTRUM: " << hyp.check_flag(ExternalFlag::IS_LHA_SPECTRUM) << "\n";

    // If you want to switch LHA you can use this API, this will reload all parameters from the LHA.
    // Be careful if you have already done calculation, Wilson coefficients will be removed from the database and need to be recalculated.
    hyp.switch_lha("lha/testinput_thdm.lha", config);

    return 0;
}
