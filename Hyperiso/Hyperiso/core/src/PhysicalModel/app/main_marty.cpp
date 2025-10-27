#include "HyperisoMaster.h"
#include "WilsonInterface.h"
#include "ParameterSetter.h"
#include <iostream>
#include "BlockProxy.h"

int main() {
    HyperisoMaster hyp = HyperisoMaster();
    Config config;
    config.flags[ExternalFlag::USE_MARTY] = true;
    config.model = Model::SM;
    config.mty_model_name = "SM";
    config.mty_model_path = project_assets_root.data() + std::string("input_files/marty_model/sm.h");

    hyp.init("default/lha/testInput.flha", config);
    LOG_INFO("HyperisoMaster initialized");

    BlockProxy().log_block(ParameterType::SM, "SMINPUTS");
    auto wi = WilsonInterface();
    LOG_INFO("WilsonInterface created");
    wi.build(WilsonBuildConfig{
        {WGroup::B},
        81,
        4.7,
        QCDOrder::LO} 
    );

    LOG_INFO("Interface built");
    
    ParameterSetter ps;
    ParameterProxy pa (ParameterType::WILSON);

    LOG_INFO("Parameters created");
    LOG_INFO("C7(mu_h) at LO =", wi.getR(WGroup::B, WCoef::C7, QCDOrder::LO, ContributionType::TOTAL));
    LOG_INFO("AGAIN");
    LOG_INFO("C7(mu_h) at NLO =", wi.getR(WGroup::B, WCoef::C7, QCDOrder::NLO, ContributionType::TOTAL));
    LOG_INFO("C7(mu_h) at NNLO =", wi.getR(WGroup::B, WCoef::C7, QCDOrder::NNLO, ContributionType::TOTAL));

    LOG_INFO("C7(mu_h) full =", wi.getFR(WGroup::B, WCoef::C7, QCDOrder::NNLO, ContributionType::TOTAL));

    LOG_INFO("C7(mu_W) at LO =", wi.getM(WGroup::B, WCoef::C7, QCDOrder::LO, ContributionType::TOTAL));
    LOG_INFO("C7(mu_W) at NLO =", wi.getM(WGroup::B, WCoef::C7, QCDOrder::NLO, ContributionType::TOTAL));
    LOG_INFO("C7(mu_W) at NNLO =", wi.getM(WGroup::B, WCoef::C7, QCDOrder::NNLO, ContributionType::TOTAL));

    LOG_INFO("C7(mu_W) full =", wi.getFM(WGroup::B, WCoef::C7, QCDOrder::NNLO, ContributionType::TOTAL));

    return 0;
}