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

    WilsonBuildConfig config_k{{WGroup::MESON_MIXING}, 81, 4.18, QCDOrder::LO};

    std::cout << "here" << std::endl;
    wi.build(config_k);


    
    LOG_INFO("Interface built");

    BlockProxy().log_all_blocks(ParameterType::WILSON);
    
    LOG_INFO("C_BD_1(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::C_BD_1, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CT_BD_1(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::CT_BD_1, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_BD_2(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::C_BD_2, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CT_BD_2(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::CT_BD_2, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_BD_3(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::C_BD_3, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CT_BD_3(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::CT_BD_3, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_BD_4(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::C_BD_4, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_BD_5(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::C_BD_5, QCDOrder::LO, ContributionType::SM));

    LOG_INFO("C_BS_1(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::C_BS_1, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CT_BS_1(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::CT_BS_1, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_BS_2(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::C_BS_2, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CT_BS_2(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::CT_BS_2, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_BS_3(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::C_BS_3, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CT_BS_3(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::CT_BS_3, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_BS_4(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::C_BS_4, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_BS_5(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::C_BS_5, QCDOrder::LO, ContributionType::SM));

    LOG_INFO("C_SD_1(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::C_SD_1, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CT_SD_1(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::CT_SD_1, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_SD_2(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::C_SD_2, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CT_SD_2(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::CT_SD_2, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_SD_3(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::C_SD_3, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CT_SD_3(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::CT_SD_3, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_SD_4(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::C_SD_4, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_SD_5(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::C_SD_5, QCDOrder::LO, ContributionType::SM));

    LOG_INFO("C_CU_1(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::C_CU_1, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CT_CU_1(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::CT_CU_1, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_CU_2(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::C_CU_2, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CT_CU_2(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::CT_CU_2, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_CU_3(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::C_CU_3, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CT_CU_3(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::CT_CU_3, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_CU_4(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::C_CU_4, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_CU_5(mu_W) at LO =", wi.getM(WGroup::MESON_MIXING, WCoef::C_CU_5, QCDOrder::LO, ContributionType::SM));

    
    
    LOG_INFO("C_BD_1(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::C_BD_1, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CT_BD_1(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::CT_BD_1, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_BD_2(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::C_BD_2, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CT_BD_2(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::CT_BD_2, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_BD_3(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::C_BD_3, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CT_BD_3(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::CT_BD_3, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_BD_4(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::C_BD_4, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_BD_5(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::C_BD_5, QCDOrder::LO, ContributionType::SM));

    LOG_INFO("C_BS_1(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::C_BS_1, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CT_BS_1(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::CT_BS_1, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_BS_2(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::C_BS_2, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CT_BS_2(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::CT_BS_2, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_BS_3(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::C_BS_3, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CT_BS_3(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::CT_BS_3, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_BS_4(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::C_BS_4, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_BS_5(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::C_BS_5, QCDOrder::LO, ContributionType::SM));

    LOG_INFO("C_SD_1(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::C_SD_1, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CT_SD_1(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::CT_SD_1, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_SD_2(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::C_SD_2, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CT_SD_2(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::CT_SD_2, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_SD_3(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::C_SD_3, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CT_SD_3(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::CT_SD_3, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_SD_4(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::C_SD_4, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_SD_5(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::C_SD_5, QCDOrder::LO, ContributionType::SM));

    LOG_INFO("C_CU_1(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::C_CU_1, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CT_CU_1(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::CT_CU_1, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_CU_2(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::C_CU_2, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CT_CU_2(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::CT_CU_2, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_CU_3(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::C_CU_3, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("CT_CU_3(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::CT_CU_3, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_CU_4(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::C_CU_4, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C_CU_5(mu_h) at LO =", wi.getR(WGroup::MESON_MIXING, WCoef::C_CU_5, QCDOrder::LO, ContributionType::SM));

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
