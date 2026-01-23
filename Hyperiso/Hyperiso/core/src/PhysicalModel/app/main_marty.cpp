#include "HyperisoMaster.h"
#include "WilsonInterface.h"
#include "ParameterSetter.h"
#include <iostream>
#include "BlockProxy.h"

int main() {
    HyperisoMaster hyp = HyperisoMaster();
    HyperisoConfig config;
    // config.flags[ExternalFlag::USE_MARTY] = true;
    config.model = Model::MARTY;
    config.mty_model_name = "THDM";
    config.mty_model_path = project_assets_root.data() + std::string("input_files/marty_model/thdm.h");

    // hyp.init("default/lha/testInput.flha", config);
    // hyp.init("lha/camilia.flha", config);
    hyp.init("spectrum/testinput_thdm.lha", config);
    // hyp.init("lha/testinput_thdm.lha", config);
    // hyp.init("lha/testInput.slha", config);
    LOG_INFO("HyperisoMaster initialized");

    BlockProxy().log_block(ParameterType::SM, "SMINPUTS");
    auto wi = WilsonInterface();
    LOG_INFO("WilsonInterface created");
    wi.build(WilsonBuildConfig{
        {WGroup::B},
        81,
        4.7,
        QCDOrder::NNLO} 
    );

    LOG_INFO("Interface built");
    
    ParameterSetter ps;
    ParameterProxy pa (ParameterType::WILSON);

    LOG_INFO("alpha =", QCDHelper::alpha_s(pa("EW_SCALE", 1 )));
    LOG_INFO("C9(mu_W) BSM at LO =", wi.getM(WGroup::B, WCoef::C9, QCDOrder::LO, ContributionType::BSM));
    LOG_INFO("C9(mu_W) BSM at NLO =", wi.getM(WGroup::B, WCoef::C9, QCDOrder::NLO, ContributionType::BSM));
    LOG_INFO("C9(mu_W) BSM at NNLO =", wi.getM(WGroup::B, WCoef::C9, QCDOrder::NNLO, ContributionType::BSM));
    LOG_INFO("C9(mu_W) SM at LO =", wi.getM(WGroup::B, WCoef::C9, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C9(mu_W) SM at NLO =", wi.getM(WGroup::B, WCoef::C9, QCDOrder::NLO, ContributionType::SM));
    LOG_INFO("C9(mu_W) SM at NNLO =", wi.getM(WGroup::B, WCoef::C9, QCDOrder::NNLO, ContributionType::SM));

    LOG_INFO("C7(mu_W) SM at LO =", wi.getM(WGroup::B, WCoef::C7, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C7(mu_W) BSM at LO =", wi.getM(WGroup::B, WCoef::C7, QCDOrder::LO, ContributionType::BSM));

    LOG_INFO("C7(mu_W) SM at NLO =", wi.getM(WGroup::B, WCoef::C7, QCDOrder::NLO, ContributionType::SM));
    LOG_INFO("C7(mu_W) BSM at NLO =", wi.getM(WGroup::B, WCoef::C7, QCDOrder::NLO, ContributionType::BSM));

    LOG_INFO("C7(mu_W) at LO =", wi.getM(WGroup::B, WCoef::C7, QCDOrder::LO, ContributionType::TOTAL));
    LOG_INFO("C7(mu_W) at NLO =", wi.getM(WGroup::B, WCoef::C7, QCDOrder::NLO, ContributionType::TOTAL));
    LOG_INFO("C7(mu_W) at NNLO =", wi.getM(WGroup::B, WCoef::C7, QCDOrder::NNLO, ContributionType::TOTAL));

    LOG_INFO("C7(mu_W) SM full =", wi.getFM(WGroup::B, WCoef::C7, QCDOrder::NNLO, ContributionType::SM));
    LOG_INFO("C7(mu_W) BSM full =", wi.getFM(WGroup::B, WCoef::C7, QCDOrder::NNLO, ContributionType::BSM));
    LOG_INFO("C7(mu_W) full =", wi.getFM(WGroup::B, WCoef::C7, QCDOrder::NNLO, ContributionType::TOTAL));

    return 0;
}