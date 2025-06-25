#include "Logger.h"
#include <iostream>
#include <cassert>
#include "ObservableInterface.h"
#include "HyperisoMaster.h"
#include "config.hpp"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);
    HyperisoMaster hyp;
    Config config;
    config.model = Model::SM;
    config.flags[ExternalFlag::USE_MARTY] = false;
    config.flags[ExternalFlag::HAS_WILSON_INPUT] = true;
    config.mty_model_name = "SM";
    config.mty_model_path = project_assets_root.data() + std::string("input_files/marty_model/sm.h");
    hyp.init("lha/testInput.flha", config);
    
    LOG_INFO("HyperisoMaster initialized");

    LOG_INFO(ParameterProvider()({ParameterType::WILSON, "EW_SCALE", 1}));

    // ParameterProvider pp;
    // ParameterSetter ps;
    // WilsonInterface wi;

    // WilsonBuildConfig wilson_config;
    // wilson_config.groups = {WGroup::B, WGroup::BPrime, WGroup::BScalar, WGroup::BCC};
    // wilson_config.matching_scale = 81;
    // wilson_config.hadronic_scale = pp({ParameterType::SM, "QCD", {5, 3}}) / 2;
    // wilson_config.hadronic_scale = 42;
    // wilson_config.order = QCDOrder::LO;
    // wi.build(wilson_config);

    // BlockProvider().log_block(ParameterType::SM, "VCKM");
    // LOG_INFO(ParameterProvider()({ParameterType::WILSON, "WPARAM_MATCH_SM", 6}));
    // LOG_INFO(ParameterProvider()({ParameterType::WILSON, "WPARAM_MATCH_SM", 1}));

    // LOG_INFO("--------------- Matching --------------");
    // for (auto [g_id, _] : GroupMapper::mapping()) {
    //     for (auto c_id : WCoefMapper::get_group(g_id)) {
    //         LOG_INFO(WCoefMapper::str(c_id), "=", wi.getM(g_id, c_id, QCDOrder::LO, ContributionType::SM));
    //     }
    // }

    // LOG_INFO("--------------- Running --------------");
    // for (auto [g_id, _] : GroupMapper::mapping()) {
    //     for (auto c_id : WCoefMapper::get_group(g_id)) {
    //         LOG_INFO(WCoefMapper::str(c_id), "=", wi.getR(g_id, c_id, QCDOrder::LO, ContributionType::SM));
    //     }
    // }

    return 0;
}