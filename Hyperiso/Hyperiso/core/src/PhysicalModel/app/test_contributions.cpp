#include "HyperisoMaster.h"
#include "WilsonBuilder.h"
#include "BlockProxy.h"
#include <iostream>
#include <cassert>

int main() {
    HyperisoMaster hyp = HyperisoMaster();
    Config config;
    config.model = Model::THDM;
    config.flags[ExternalFlag::USE_MARTY] = true; // TODO : Théo not happy
    config.mty_model_name = "THDM";
    config.mty_model_path = std::string(project_tp_root.data()) + "MARTY/src/MARTY/src/marty/models/thdm.h";

    hyp.init("lha/testinput_thdm.lha", config);

    WilsonBuildConfig wilson_config;
    wilson_config.groups = {GroupMapper::to_id(WGroup::MESON_MIXING)};
    wilson_config.matching_scale = 160.0;
    wilson_config.hadronic_scale = 3;
    wilson_config.order = QCDOrder::LO;
    WilsonBuilder builder {wilson_config};

    BlockProxy bp;

    LOG_INFO("------------- Matching blocks -------------");

    // bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::MATCHING));
    // bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::BPrime, ScaleType::MATCHING));
    // bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::BScalar, ScaleType::MATCHING));
    // bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::CC_bc, ScaleType::MATCHING));
    bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING));

    LOG_INFO("------------- Hadronic blocks -------------");

    // bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC));
    // bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC, WilsonBasis::B_TRADITIONAL));
    // bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::BPrime, ScaleType::HADRONIC));
    // bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::BScalar, ScaleType::HADRONIC));
    // bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::CC_bc, ScaleType::HADRONIC));
    bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::MESON_MIXING, ScaleType::HADRONIC));

    LOG_INFO("------------- Wilson Parameter blocks -------------");
    bp.log_block(ParameterType::SM, "QCD");
    bp.log_block(ParameterType::SM, "VCKM");
    bp.log_block(ParameterType::WILSON, "B_SCALE");
    bp.log_block(ParameterType::WILSON, "WPARAM_SI_SM");
    bp.log_block(ParameterType::WILSON, "WPARAM_MATCH_SM");
    bp.log_block(ParameterType::WILSON, "WPARAM_RUN_SM");
    bp.log_block(ParameterType::WILSON, "UM_MATRIX_5");
    bp.log_block(ParameterType::WILSON, "UM_MATRIX_4");

    return 0;
}
