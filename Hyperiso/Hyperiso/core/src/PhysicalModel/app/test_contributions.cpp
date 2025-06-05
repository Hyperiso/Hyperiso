#include "HyperisoMaster.h"
#include "WilsonBuilder.h"
#include <iostream>
#include <cassert>

int main() {
    HyperisoMaster hyp = HyperisoMaster();
    Config config;
    config.model = Model::SUSY;
    // config.flags[ExternalFlag::USE_MARTY] = true; // TODO : Théo not happy
    // config.mty_model_name = "THDM";
    // config.mty_model_path = std::string(project_tp_root.data()) + "MARTY/src/MARTY/src/marty/models/thdm.h";

    hyp.init("lha/testInput.slha", config);

    WilsonBuildConfig wilson_config;
    wilson_config.groups = {WGroup::B, WGroup::BPrime, WGroup::BScalar, WGroup::BCC};
    wilson_config.matching_scale = 85.0;
    wilson_config.hadronic_scale = 4.5;
    wilson_config.order = QCDOrder::NNLO;
    WilsonBuilder builder {wilson_config};

    BlockProxy bp;

    LOG_INFO("------------- Matching blocks -------------");

    bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::MATCHING));
    bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::BPrime, ScaleType::MATCHING));
    bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::BScalar, ScaleType::MATCHING));
    bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::BCC, ScaleType::MATCHING));

    LOG_INFO("------------- Hadronic blocks -------------");

    bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC));
    bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC, WilsonBasis::B_TRADITIONAL));
    bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::BPrime, ScaleType::HADRONIC));
    bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::BScalar, ScaleType::HADRONIC));
    bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::BCC, ScaleType::HADRONIC));

    return 0;
}
