#include <iostream>

#include "MemoryManager.h"
#include "Parameters.h"
#include "WilsonInterface.h"

int main() {
    auto hyp = HyperisoMaster();
    Config config_hyp;
    config_hyp.flags[ExternalFlag::USE_MARTY] = false;
    config_hyp.model = Model::SM;
    config_hyp.mty_model_name = "ZPrime";
    config_hyp.mty_model_path = "/home/theo/hyperiso/Assets/input_files/marty_model/ZPrime.h";

    hyp.init("lha/camilia.flha", config_hyp); // Initialize program manager with LHA file

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
    return 0;
}