#include "HyperisoMaster.h"
#include "WilsonInterface.h"
#include <iostream>

int main(){
    HyperisoMaster hyp = HyperisoMaster();
    Config config;
    config.model = Model::SM;

    hyp.init("default/lha/testInput.flha", config);
    LOG_INFO("HyperisoMaster initialized");

    auto wi = WilsonInterface(); // Initialize interface and build the required groups
    LOG_INFO("WilsonInterface created");
    wi.build(
        {WGroup::B, WGroup::BPrime},                            // Coefficient groups
        81,     // Matching scale
        42,                             // Hadronic scale
        QCDOrder::NNLO                                          // QCD Order
    );

    LOG_INFO("Interface built");

    std::cout << "C7 matching (LO + NLO + NNLO) : " << wi.getFullMatchingCoefficient(WGroup::B, WCoef::C7, QCDOrder::NNLO)     << '\n';
    std::cout << "C7 hadronic (LO + NLO + NNLO) : " << wi.getFullRunCoefficient(WGroup::B, WCoef::C7, QCDOrder::NNLO)          << '\n';
    std::cout << "CP7 matching (LO) : "             << wi.getFullMatchingCoefficient(WGroup::BPrime, WCoef::CP7, QCDOrder::LO) << '\n';
    std::cout << "CP7 hadronic (LO) : "             << wi.getFullRunCoefficient(WGroup::BPrime, WCoef::CP7, QCDOrder::LO)      << '\n';

    // Retrieve coefficient values
    // std::ofstream fs {"out.txt"};
    // fs << "C7 matching (LO + NLO + NNLO) : " << wi.getFullMatchingCoefficient(WGroup::B, WCoef::C7, QCDOrder::NNLO)     << '\n';
    // fs << "C7 hadronic (LO + NLO + NNLO) : " << wi.getFullRunCoefficient(WGroup::B, WCoef::C7, QCDOrder::NNLO)          << '\n';
    // fs << "CP7 matching (LO) : "             << wi.getFullMatchingCoefficient(WGroup::BPrime, WCoef::CP7, QCDOrder::LO) << '\n';
    // fs << "CP7 hadronic (LO) : "             << wi.getFullRunCoefficient(WGroup::BPrime, WCoef::CP7, QCDOrder::LO)      << '\n';
    // fs.close();

    return 0;
}