#include <iostream>

#include "HyperisoMaster.h"
#include "Include.h"
#include "Logger.h"
#include "QCDProvider.h"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    // All informations for the creation of the HyperisoMaster are provided in the base_core_example.cpp example.
    HyperisoConfig config;
    HyperisoMaster hyp;
    hyp.init("lha/si_input.flha", config);

    // This API allows to run alpha_s and masses.
    QCDProvider qcd_runner;

    // Define the configuration, with the scale at which you want to get alpha_s, and the mass scheme.
    AlphasConfig ac(81.0, MassType::MSBAR, MassType::POLE);

    // The operator() API returns alpha_s at 81 GeV.
    std::cout << "alpha_s(81 GeV) = " << qcd_runner(ac) << "\n";

    // Configuration to run masses. Please specify the PDG id of the mass to run. Here 6 is the top quark.
    MassConfig mc(6, 81.0, MassType::MSBAR, MassType::POLE);

    // The operator() API returns the top mass at 81 GeV.
    std::cout << "m_top(81 GeV) = " << qcd_runner(mc) << "\n";

    return 0;
}
