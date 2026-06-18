#include <filesystem>
#include <iostream>
#include <map>

#include "HyperisoMaster.h"
#include "Include.h"
#include "Logger.h"
#include "ParameterProvider.h"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    // All informations for the creation of the HyperisoMaster are provided in the base_core_example.cpp example.
    HyperisoConfig config;
    HyperisoMaster hyp;

    const auto parameters_values = std::filesystem::path(__FILE__).parent_path() / "inputs" / "sm_param_test.yaml";

    // This API allows to set special paths for Hyperiso's inputs.
    // One can change parameters through this API, Observables, Observables correlations, etc.
    // We recommend not to change DEFAULT values, even if this API allows it.
    // User inputs take the form of .yml/.yaml files, where you can specify the value and uncertainties of your parameters.
    hyp.pre_init_set_paths({
        {APIPath::USER_SM_PARAMS, parameters_values.string()}
    });

    // One can also add specials blocks through a pre-init. This is useful for BSM models with new groups not in the FLHA convention.
    // If it is a matrix block, with multiple IDs, use the options available within the API.
    hyp.pre_init_add_block("SPECIAL_BLOCK");

    // If you already have MARTY installed on your machine, please specify the path through the following API.
    // hyp.pre_init_set_marty_path("Path/To/Marty");

    hyp.init("lha/si_input.flha", config);

    std::cout << "SM::MASS[18] after user override: "
              << ParameterProvider()(ParamId(ParameterType::SM, "MASS", LhaID(18))) << "\n";

    return 0;
}
