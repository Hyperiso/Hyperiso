#include "HyperisoMaster.h"
#include "WilsonInterface.h"
#include "ParameterSetter.h"
#include "BlockProxy.h"	
#include <iostream>

int main(){
    HyperisoMaster hyp = HyperisoMaster();
    HyperisoConfig config;
    config.model = Model::SM;

    hyp.init("lha/testInput.flha", config);
    LOG_INFO("HyperisoMaster initialized");

    auto wi = WilsonInterface(); // Initialize interface and build the required groups
    LOG_INFO("WilsonInterface created");
    WilsonBuildConfig config_b{{WGroup::B}, 81, 4.18, QCDOrder::NLO};

    std::cout << "here" << std::endl;
    wi.build(config_b);

    WilsonBuildConfig config_k{{WGroup::K}, 81, 4.18, QCDOrder::NLO};

    std::cout << "here" << std::endl;
    wi.addWilsonGroup(config_k);


    
    LOG_INFO("Interface built");

    BlockProxy().log_all_blocks(ParameterType::WILSON);
    
    LOG_INFO("CK9(mu_W) at LO =", wi.getM(WGroup::K, WCoef::CK9, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CK10(mu_W) at LO =", wi.getM(WGroup::K, WCoef::CK10, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CKQ1(mu_W) at LO =", wi.getM(WGroup::K, WCoef::CKQ1, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CKQ2(mu_W) at LO =", wi.getM(WGroup::K, WCoef::CKQ2, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CPK9(mu_W) at LO =", wi.getM(WGroup::K, WCoef::CPK9, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CPK10(mu_W) at LO =", wi.getM(WGroup::K, WCoef::CPK10, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CPKQ1(mu_W) at LO =", wi.getM(WGroup::K, WCoef::CPKQ1, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CPKQ2(mu_W) at LO =", wi.getM(WGroup::K, WCoef::CPKQ2, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CK_L(mu_W) at LO =", wi.getM(WGroup::K, WCoef::CK_L, QCDOrder::LO, ContributionType::SM));

    LOG_INFO("CK9(mu_h) at LO =", wi.getR(WGroup::K, WCoef::CK9, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CK10(mu_h) at LO =", wi.getR(WGroup::K, WCoef::CK10, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CKQ1(mu_h) at LO =", wi.getR(WGroup::K, WCoef::CKQ1, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CKQ2(mu_h) at LO =", wi.getR(WGroup::K, WCoef::CKQ2, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CPK9(mu_h) at LO =", wi.getR(WGroup::K, WCoef::CPK9, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CPK10(mu_h) at LO =", wi.getR(WGroup::K, WCoef::CPK10, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CPKQ1(mu_h) at LO =", wi.getR(WGroup::K, WCoef::CPKQ1, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CPKQ2(mu_h) at LO =", wi.getR(WGroup::K, WCoef::CPKQ2, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CK_L(mu_h) at LO =", wi.getR(WGroup::K, WCoef::CK_L, QCDOrder::LO, ContributionType::SM));
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
