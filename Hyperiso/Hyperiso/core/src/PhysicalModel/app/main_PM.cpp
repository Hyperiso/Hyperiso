#include "HyperisoMaster.h"
#include "WilsonInterface.h"
#include "ParameterSetter.h"
#include "BlockProxy.h"	
#include <iostream>

int main(){
    HyperisoMaster hyp = HyperisoMaster();
    HyperisoConfig config;
    config.model = Model::SM;

    hyp.init("lha/testInput.slha", config);
    LOG_INFO("HyperisoMaster initialized");

    auto wi = WilsonInterface(); // Initialize interface and build the required groups
    LOG_INFO("WilsonInterface created");

    
    wi.build(
        WilsonBuildConfig{{WGroup::B},                            // Coefficient groups
        81,     // Matching scale
        4.7,                             // Hadronic scale
        QCDOrder::NNLO}                                          // QCD Order
    );

    BlockProxy().log_all_blocks(ParameterType::WILSON);
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
