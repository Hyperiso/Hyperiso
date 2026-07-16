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
    WilsonBuildConfig config_b{{WGroup::CC_bc, WGroup::CC_cd, WGroup::CC_su, WGroup::CC_du}, 81, 4.18, QCDOrder::NNLO};

    std::cout << "hereeee" << std::endl;
    wi.build(config_b);
    wi.addWilsonGroup({{WGroup::CC_cs}, 81, 4.18, QCDOrder::NNLO});

    WilsonBuildConfig config_k{{WGroup::CC_bu}, 81, 4.18, QCDOrder::NNLO};

    std::cout << "here" << std::endl;
    wi.addWilsonGroup(config_k);


    
    LOG_INFO("Interface built");

    // BlockProxy().log_all_blocks(ParameterType::WILSON);
    
    LOG_INFO("C_S1_bc(mu_W) at LO =", wi.getM(WGroup::CC_bc, WCoef::C_S1_bc, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_S2_bc(mu_W) at LO =", wi.getM(WGroup::CC_bc, WCoef::C_S2_bc, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_V1_bc(mu_W) at LO =", wi.getM(WGroup::CC_bc, WCoef::C_V1_bc, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_V2_bc(mu_W) at LO =", wi.getM(WGroup::CC_bc, WCoef::C_V2_bc, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_T_bc(mu_W) at LO =", wi.getM(WGroup::CC_bc, WCoef::C_T_bc, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_S1_bu(mu_W) at LO =", wi.getM(WGroup::CC_bu, WCoef::C_S1_bu, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_S2_bu(mu_W) at LO =", wi.getM(WGroup::CC_bu, WCoef::C_S2_bu, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_V1_bu(mu_W) at LO =", wi.getM(WGroup::CC_bu, WCoef::C_V1_bu, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_V2_bu(mu_W) at LO =", wi.getM(WGroup::CC_bu, WCoef::C_V2_bu, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_T_bu(mu_W) at LO =", wi.getM(WGroup::CC_bu, WCoef::C_T_bu, QCDOrder::LO, ContributionType::SM));

    LOG_INFO("C_S1_bc(mu_W) BSM at LO =", wi.getM(WGroup::CC_bc, WCoef::C_S1_bc, QCDOrder::LO, ContributionType::BSM));
    LOG_INFO("C_S2_bc(mu_W) BSM at LO =", wi.getM(WGroup::CC_bc, WCoef::C_S2_bc, QCDOrder::LO, ContributionType::BSM));
    LOG_INFO("C_V1_bc(mu_W) BSM at LO =", wi.getM(WGroup::CC_bc, WCoef::C_V1_bc, QCDOrder::LO, ContributionType::BSM));
    LOG_INFO("C_V2_bc(mu_W) BSM at LO =", wi.getM(WGroup::CC_bc, WCoef::C_V2_bc, QCDOrder::LO, ContributionType::BSM));
    LOG_INFO("C_T_bc(mu_W) BSM at LO =", wi.getM(WGroup::CC_bc, WCoef::C_T_bc, QCDOrder::LO, ContributionType::BSM));
    LOG_INFO("C_S1_bu(mu_W) BSM at LO =", wi.getM(WGroup::CC_bu, WCoef::C_S1_bu, QCDOrder::LO, ContributionType::BSM));
    LOG_INFO("C_S2_bu(mu_W) BSM at LO =", wi.getM(WGroup::CC_bu, WCoef::C_S2_bu, QCDOrder::LO, ContributionType::BSM));
    LOG_INFO("C_V1_bu(mu_W) BSM at LO =", wi.getM(WGroup::CC_bu, WCoef::C_V1_bu, QCDOrder::LO, ContributionType::BSM));
    LOG_INFO("C_V2_bu(mu_W) BSM at LO =", wi.getM(WGroup::CC_bu, WCoef::C_V2_bu, QCDOrder::LO, ContributionType::BSM));
    LOG_INFO("C_T_bu(mu_W) BSM at LO =", wi.getM(WGroup::CC_bu, WCoef::C_T_bu, QCDOrder::LO, ContributionType::BSM));

    LOG_INFO("C_S1_bc(mu_h) at LO =", wi.getR(WGroup::CC_bc, WCoef::C_S1_bc, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_S2_bc(mu_h) at LO =", wi.getR(WGroup::CC_bc, WCoef::C_S2_bc, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_V1_bc(mu_h) at LO =", wi.getR(WGroup::CC_bc, WCoef::C_V1_bc, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_V2_bc(mu_h) at LO =", wi.getR(WGroup::CC_bc, WCoef::C_V2_bc, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_T_bc(mu_h) at LO =", wi.getR(WGroup::CC_bc, WCoef::C_T_bc, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_S1_bu(mu_h) at LO =", wi.getR(WGroup::CC_bu, WCoef::C_S1_bu, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_S2_bu(mu_h) at LO =", wi.getR(WGroup::CC_bu, WCoef::C_S2_bu, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_V1_bu(mu_h) at LO =", wi.getR(WGroup::CC_bu, WCoef::C_V1_bu, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_V2_bu(mu_h) at LO =", wi.getR(WGroup::CC_bu, WCoef::C_V2_bu, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_T_bu(mu_h) at LO =", wi.getR(WGroup::CC_bu, WCoef::C_T_bu, QCDOrder::LO, ContributionType::SM));
    // LOG_INFO("CQ1(mu_h) at LO =", wi.getR(WGroup::BScalar, WCoef::CQ1_MU, QCDOrder::LO));
    // LOG_INFO("CQ2(mu_h) at LO =", wi.getR(WGroup::BScalar, WCoef::CQ2_MU, QCDOrder::LO));
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
