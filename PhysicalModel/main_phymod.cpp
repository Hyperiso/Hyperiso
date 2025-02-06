#include <iostream>

#include "MemoryManager.h"
#include "Parameters.h"
#include "WilsonInterface.h"

int main() {
    auto mm = MemoryManager::GetInstance();
    mm->init("Test/InputFiles/testInput.flha", Model::SM); // Initialize program manager with LHA file

    auto wi = WilsonInterface(); // Initialize interface and build the required groups
    wi.build(
        {WGroup::B, WGroup::BPrime},                            // Coefficient groups
        2 * Parameters::Get(ParameterType::SM, "MASS", 24),     // Matching scale
        QCDHelper::mass_b_1S() / 2,                             // Hadronic scale
        QCDOrder::NNLO                                          // QCD Order
    );

    // Retrieve coefficient values
    std::ofstream fs {"out.txt"};
    fs << "C7 matching (LO + NLO + NNLO) : " << wi.getFullMatchingCoefficient(WGroup::B, WCoef::C7, QCDOrder::NNLO)     << '\n';
    fs << "C7 hadronic (LO + NLO + NNLO) : " << wi.getFullRunCoefficient(WGroup::B, WCoef::C7, QCDOrder::NNLO)          << '\n';
    fs << "CP7 matching (LO) : "             << wi.getFullMatchingCoefficient(WGroup::BPrime, WCoef::CP7, QCDOrder::LO) << '\n';
    fs << "CP7 hadronic (LO) : "             << wi.getFullRunCoefficient(WGroup::BPrime, WCoef::CP7, QCDOrder::LO)      << '\n';
    fs.close();

    return 0;
}