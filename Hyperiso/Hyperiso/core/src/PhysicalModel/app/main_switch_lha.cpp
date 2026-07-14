#include <iostream>

#include "MemoryManager.h"
#include "Parameters.h"
#include "WilsonInterface.h"
#include "BlockProxy.h"

int main(int argc, char** argv) {
    auto hyp = HyperisoMaster();
    HyperisoConfig config_hyp;
    // config_hyp.flags[ExternalFlag::USE_MARTY] = false;
    config_hyp.model = Model::THDM;
    config_hyp.mty_model_name = "ZPrime";
    config_hyp.mty_model_path = argc > 2 ? argv[2] : "Assets/input_files/marty_model/ZPrime.h";

    hyp.init(argc > 1 ? argv[1] : "Assets/lha/testinput_thdm.lha", config_hyp); // Initialize program manager with LHA file

    auto wi = WilsonInterface(); // Initialize interface and build the required groups

    WilsonBuildConfig config({WGroup::B}, 81, 42, QCDOrder::LO);

    wi.build(config);

    BlockProxy().log_block(ParameterType::WILSON, "BCoefficients_B_SCALE_STANDARD");

    LOG_INFO("Parameters created");
    LOG_INFO("C9(mu_h) at LO =", wi.getR(WGroup::B, WCoef::C9, QCDOrder::LO, ContributionType::TOTAL));
    LOG_INFO("AGAIN");
    LOG_INFO("C9(mu_h) at NLO =", wi.getR(WGroup::B, WCoef::C9, QCDOrder::NLO, ContributionType::TOTAL));
    LOG_INFO("C9(mu_h) at NNLO =", wi.getR(WGroup::B, WCoef::C9, QCDOrder::NNLO, ContributionType::TOTAL));
    
    LOG_INFO("C9(mu_h) full =", wi.getFR(WGroup::B, WCoef::C9, QCDOrder::NNLO, ContributionType::TOTAL));

    HyperisoConfig config_v2;

    BlockProxy().log_all_blocks(ParameterType::WILSON);
    // LOG_INFO("changing config");
    hyp.switch_lha("default/lha/testInput.flha", config_v2);

    BlockProxy().log_all_blocks(ParameterType::WILSON);
    
    return 0;
}