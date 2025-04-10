#include "HyperisoMaster.h"
#include "WilsonInterface.h"
#include "ParameterSetter.h"
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
        4.7,                             // Hadronic scale
        QCDOrder::NNLO                                          // QCD Order
    );

    LOG_INFO("Interface built");
    
    ParameterSetter ps;
    ParameterProxy pa (ParameterType::WILSON);

    // LOG_INFO("Before mutation:", pa("WPARAM_MATCH_SM", LhaID(2, 1)));
    // LOG_INFO("Before mutation:", pa("B_MATCH", LhaID(3051313, 4137, 0, 0)));
    // ps.mutate({ParameterType::SM, "MASS", 24}, 100);
    // LOG_INFO("After mutation:", pa("WPARAM_MATCH_SM", LhaID(2, 1)));
    // LOG_INFO("After mutation:", pa("B_MATCH", LhaID(3051313, 4137, 0, 0)));

    // LOG_INFO("Mutating parameter");
    // ps.mutate({ParameterType::SM, "MASS", 24}, 80.3);
    
    LOG_INFO("C7(mu_h) at LO =", wi.getR(WGroup::B, WCoef::C7, QCDOrder::LO));
    LOG_INFO("C7(mu_h) at NLO =", wi.getR(WGroup::B, WCoef::C7, QCDOrder::NLO));
    LOG_INFO("C7(mu_h) at NNLO =", wi.getR(WGroup::B, WCoef::C7, QCDOrder::NNLO));

    LOG_INFO("C7(mu_h) full =", wi.getFR(WGroup::B, WCoef::C7, QCDOrder::NNLO));

    LOG_INFO("C7(mu_W) at LO =", wi.getM(WGroup::B, WCoef::C7, QCDOrder::LO));
    LOG_INFO("C7(mu_W) at NLO =", wi.getM(WGroup::B, WCoef::C7, QCDOrder::NLO));
    LOG_INFO("C7(mu_W) at NNLO =", wi.getM(WGroup::B, WCoef::C7, QCDOrder::NNLO));

    LOG_INFO("C7(mu_W) full =", wi.getFM(WGroup::B, WCoef::C7, QCDOrder::NNLO));



    // std::cout << "C7 matching (LO + NLO + NNLO) : " << wi.getFullMatchingCoefficient(WGroup::B, WCoef::C7, QCDOrder::LO)     << '\n';
    // std::cout << "C7 hadronic (LO + NLO + NNLO) : " << wi.getFullRunCoefficient(WGroup::B, WCoef::C7, QCDOrder::LO)          << '\n';
    // std::cout << "CP7 matching (LO) : "             << wi.getFullMatchingCoefficient(WGroup::BPrime, WCoef::CP7, QCDOrder::LO) << '\n';
    // std::cout << "CP7 hadronic (LO) : "             << wi.getFullRunCoefficient(WGroup::BPrime, WCoef::CP7, QCDOrder::LO)      << '\n';

    // Retrieve coefficient values
    // std::ofstream fs {"out.txt"};
    // fs << "C7 matching (LO + NLO + NNLO) : " << wi.getFullMatchingCoefficient(WGroup::B, WCoef::C7, QCDOrder::NNLO)     << '\n';
    // fs << "C7 hadronic (LO + NLO + NNLO) : " << wi.getFullRunCoefficient(WGroup::B, WCoef::C7, QCDOrder::NNLO)          << '\n';
    // fs << "CP7 matching (LO) : "             << wi.getFullMatchingCoefficient(WGroup::BPrime, WCoef::CP7, QCDOrder::LO) << '\n';
    // fs << "CP7 hadronic (LO) : "             << wi.getFullRunCoefficient(WGroup::BPrime, WCoef::CP7, QCDOrder::LO)      << '\n';
    // fs.close();

    return 0;
}
