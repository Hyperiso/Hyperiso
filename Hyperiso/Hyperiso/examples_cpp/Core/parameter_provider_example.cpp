#include <iostream>

#include "APIAdapter.h"
#include "CorrelationProvider.h"
#include "HyperisoMaster.h"
#include "Include.h"
#include "Logger.h"
#include "ParameterProvider.h"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    // All informations for the creation of the HyperisoMaster are provided in the base_core_example.cpp example.
    HyperisoConfig config;
    HyperisoMaster hyp;
    hyp.init("lha/si_input.flha", config);

    // This creates a ParameterProvider for SM parameters. This avoids specifying the ParameterType each time.
    ParameterProvider sm_param_provider(ParameterType::SM);

    // If you want parameters from differents sections, use a no-typed ParameterProvider.
    ParameterProvider no_pt_param_provider;

    // Using the typed provider, one can test if a parameter exists using the following API.
    std::cout << "SMINPUTS[2] exists: " << sm_param_provider.exists("SMINPUTS", LhaID(2)) << "\n";

    // Using the no-typed provider, one can test if a parameter exists using a ParamId.
    std::cout << "SM::SMINPUTS[2] exists: "
              << no_pt_param_provider.exists(ParamId(ParameterType::SM, "SMINPUTS", LhaID(2))) << "\n";

    // One can get the value of a Parameter through the following API.
    // If one wants to get the uncertainty, one can change the DataType.
    std::cout << "SM::MASS[25] value: " << sm_param_provider("MASS", LhaID(25), DataType::VALUE) << "\n";

    // Same with a no-typed provider.
    std::cout << "SM::MASS[25] value with ParamId: "
              << no_pt_param_provider(ParamId(ParameterType::SM, "MASS", LhaID(25)), DataType::VALUE) << "\n";

    // One can get a copy of a Parameter class with all informations: value, id, uncertainties, mode, etc.
    auto param = no_pt_param_provider.get_parameter(ParamId(ParameterType::WILSON, "B_SCALE", LhaID(1)));
    std::cout << "WILSON::B_SCALE[1] value: " << param->get_val() << "\n";

    // One might also want to get correlation between parameters or observables.
    // This is useful for the Statistic part.
    CorrelationProvider cp;

    std::cout << "Correlation between B_K[1_1_0] and B_K[1_1_1]: "
              << cp(ParamId(ParameterType::DECAY, "B_K", LhaID(1, 1, 0)),
                    ParamId(ParameterType::DECAY, "B_K", LhaID(1, 1, 1)),
                    CorrelationProvider::CorrelationType::STAT) << "\n";

    // For observables, one can use the experiment name as first argument.
    std::cout << "Correlation CMS F_L/P_1: "
              << cp("CMS",
                    Observables::F_L_B0__KSTAR0_MU_MU,
                    Observables::P_1_B0__KSTAR0_MU_MU,
                    CorrelationProvider::CorrelationType::STAT) << "\n";

    // If one doesn't know in which section a block is in, one can use the APIAdapter proxy.
    APIAdapter api;

    std::cout << "Possible ParameterTypes for SMINPUTS:\n";
    for (auto type : api.get_type_of_block("SMINPUTS")) {
        std::cout << "  " << ParameterTypeMapper::str(type) << "\n";
    }

    // This proxy is also useful to get blocks as list of parameters.
    std::cout << "First entries of SM::MASS:\n";
    std::size_t printed = 0;
    for (const auto& [id, value] : api.get_block_infos("MASS", ParameterType::SM)) {
        if (printed++ >= 5) break;
        std::cout << "  MASS[" << id << "] = " << value << "\n";
    }

    // It can also be used to get informations on paths.
    std::cout << "Current LHA path: " << api.get_path(APIPath::LHA_PATH) << "\n";

    return 0;
}
