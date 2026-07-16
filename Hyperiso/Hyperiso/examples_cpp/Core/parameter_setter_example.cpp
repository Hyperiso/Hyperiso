#include <iostream>

#include "HyperisoMaster.h"
#include "Include.h"
#include "Logger.h"
#include "ParameterProvider.h"
#include "ParameterSetter.h"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    // All informations for the creation of the HyperisoMaster are provided in the base_core_example.cpp example.
    HyperisoConfig config;
    HyperisoMaster hyp;
    hyp.init("lha/si_input.flha", config);

    // All informations for the usage of the ParameterProvider are provided in the parameter_provider_example.cpp example.
    ParameterProvider param_provider;

    // Creation of the ParameterSetter proxy, which is meant for changing value of parameters within Hyperiso.
    // Careful with this API: some parameters are DependentParameters, or inside DependentBlocks, which cannot be modified directly.
    // If you want to modify these kind of parameters, please see the dependency_pruner_example.cpp example.
    ParameterSetter param_setter;

    const ParamId ew_scale(ParameterType::WILSON, "EW_SCALE", LhaID(1));

    std::cout << "Value before the change of the Electroweak scale within Hyperiso: "
              << param_provider(ew_scale) << "\n";

    // The following API can change the value of a Parameter. One can use a double or a scalar_t.
    param_setter.mutate(ew_scale, 160.0);

    std::cout << "Value after the change of the Electroweak scale within Hyperiso: "
              << param_provider(ew_scale) << "\n";

    return 0;
}
