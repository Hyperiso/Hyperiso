#include <iostream>

#include "DependencyPruner.h"
#include "DependantBlockInfoProvider.h"
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

    // The DependencyPruner proxy allows to unlink dependent blocks or parameters from their dependencies.
    // Useful when wanting to scan over a dependent parameter or a parameter inside a DependentBlock.
    DependencyPruner dep_pruner;

    // All informations on these proxies are given in the parameter_provider_example.cpp and parameter_setter_example.cpp examples.
    ParameterProvider pp;
    ParameterSetter ps;

    // The DependantBlockInfoProvider class allows to get information on dependent blocks.
    DependantBlockInfoProvider info;

    const ParamId vckm_11(ParameterType::SM, "VCKM", LhaID(1, 1));

    std::cout << "Is VCKM a dependent block: " << info.is_dependent_block(ParameterType::SM, "VCKM") << "\n";

    std::cout << "Sources of VCKM:\n";
    for (const auto& block : info.get_source_blocks(ParameterType::SM, "VCKM")) {
        std::cout << "  " << block << "\n";
    }

    std::cout << "Current VCKM[1,1] value: " << pp(vckm_11) << "\n";

    dep_pruner.detach_block(ParameterType::SM, "VCKM");

    // VCKM is still a dependent block (it can be reattached later), but it has no source for the moment.
    std::cout << "Is VCKM a dependent block after detach: " << info.is_dependent_block(ParameterType::SM, "VCKM") << "\n";

    std::cout << "Sources of VCKM after detach:\n";
    for (const auto& block : info.get_source_blocks(ParameterType::SM, "VCKM")) {
        std::cout << "  " << block << "\n";
    }

    ps.mutate(vckm_11, 1.0);
    std::cout << "Changed VCKM[1,1] value: " << pp(vckm_11) << "\n";

    // This allows to reattach a block to its sources, and it will get back its original dependent behaviour.
    dep_pruner.reattach_block(ParameterType::SM, "VCKM");
    std::cout << "VCKM[1,1] value after reattach: " << pp(vckm_11) << "\n";

    return 0;
}
