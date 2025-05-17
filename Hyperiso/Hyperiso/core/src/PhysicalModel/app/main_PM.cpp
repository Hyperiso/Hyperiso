#include "HyperisoMaster.h"
#include "WilsonInterface.h"
#include "ParameterSetter.h"
#include "BlockProxy.h"	
#include <iostream>

int main(){
    HyperisoMaster hyp = HyperisoMaster();
    Config config;
    config.model = Model::THDM;

    hyp.init("lha/testinput_thdm.lha", config);
    LOG_INFO("HyperisoMaster initialized");

    auto wi = WilsonInterface(); // Initialize interface and build the required groups
    LOG_INFO("WilsonInterface created");

    
    wi.addWilsonGroup(WGroup::B);
    wi.init_group_matching(WGroup::B, QCDOrder::LO);
    wi.init_group_hadronic(WGroup::B, QCDOrder::LO);

    wi.addWilsonGroup(WGroup::BScalar);
    wi.init_group_matching(WGroup::BScalar, QCDOrder::LO);
    wi.init_group_hadronic(WGroup::BScalar, QCDOrder::LO);

    wi.addWilsonGroup(WGroup::BPrime);
    wi.init_group_matching(WGroup::BPrime, QCDOrder::LO);
    wi.init_group_hadronic(WGroup::BPrime, QCDOrder::LO);
    wi.post_init();
    // wi.build(
    //     WilsonBuildConfig{{WGroup::B},                            // Coefficient groups
    //     81,     // Matching scale
    //     4.7,                             // Hadronic scale
    //     QCDOrder::NNLO}                                          // QCD Order
    // );

    
    LOG_INFO("Interface built");

    // BlockProxy().log_all_blocks(ParameterType::WILSON);
    
    // LOG_INFO("C1(mu_h) at LO =", wi.getR(WGroup::B, WCoef::C1, QCDOrder::LO));
    // LOG_INFO("C2(mu_h) at LO =", wi.getR(WGroup::B, WCoef::C2, QCDOrder::LO));
    // LOG_INFO("C3(mu_h) at LO =", wi.getR(WGroup::B, WCoef::C3, QCDOrder::LO));
    // LOG_INFO("C4(mu_h) at LO =", wi.getR(WGroup::B, WCoef::C4, QCDOrder::LO));
    // LOG_INFO("C5(mu_h) at LO =", wi.getR(WGroup::B, WCoef::C5, QCDOrder::LO));
    // LOG_INFO("C6(mu_h) at LO =", wi.getR(WGroup::B, WCoef::C6, QCDOrder::LO));
    // LOG_INFO("C8(mu_h) at LO =", wi.getR(WGroup::B, WCoef::C8, QCDOrder::LO));
    // LOG_INFO("C9(mu_h) at LO =", wi.getR(WGroup::B, WCoef::C9, QCDOrder::LO));
    // LOG_INFO("C10(mu_h) at LO =", wi.getR(WGroup::B, WCoef::C10, QCDOrder::LO));

    // LOG_INFO("CQ1(mu_h) at LO =", wi.getR(WGroup::BScalar, WCoef::CQ1, QCDOrder::LO));
    // LOG_INFO("CQ2(mu_h) at LO =", wi.getR(WGroup::BScalar, WCoef::CQ2, QCDOrder::LO));
    // LOG_INFO("C10(mu_h) at LO =", wi.getR(WGroup::B, WCoef::C10, QCDOrder::LO));

    // LOG_INFO("C7(mu_h) at LO =", wi.getR(WGroup::B, WCoef::C7, QCDOrder::LO));
    // LOG_INFO("C7(mu_h) at NLO =", wi.getR(WGroup::B, WCoef::C7, QCDOrder::NLO));
    // LOG_INFO("C7(mu_h) at NNLO =", wi.getR(WGroup::B, WCoef::C7, QCDOrder::NNLO));

    // LOG_INFO("C7(mu_h) full =", wi.getFR(WGroup::B, WCoef::C7, QCDOrder::NNLO));

    // LOG_INFO("C7(mu_W) at LO =", wi.getM(WGroup::B, WCoef::C7, QCDOrder::LO));
    // LOG_INFO("C7(mu_W) at NLO =", wi.getM(WGroup::B, WCoef::C7, QCDOrder::NLO));
    // LOG_INFO("C7(mu_W) at NNLO =", wi.getM(WGroup::B, WCoef::C7, QCDOrder::NNLO));

    // LOG_INFO("C7(mu_W) full =", wi.getFM(WGroup::B, WCoef::C7, QCDOrder::NNLO));

    return 0;
}
