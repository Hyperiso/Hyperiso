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
    config.flags[ExternalFlag::USE_MARTY] = true;
    config.mty_model_name = "SM";
    config.mty_model_path = project_assets_root.data() + std::string("input_files/marty_model/sm.h");
    hyp.init("lha/testInput.flha", config);
    
    LOG_INFO("HyperisoMaster initialized");

    ParameterProvider pp;
    ParameterSetter ps;
    WilsonInterface wi;

    WilsonBuildConfig wilson_config;
    wilson_config.groups = {WGroup::B, WGroup::BPrime};
    wilson_config.matching_scale = 2. * pp({ParameterType::SM, "MASS", 24});
    wilson_config.hadronic_scale = pp({ParameterType::SM, "QCD", {5, 3}}) / 2;
    wilson_config.order = QCDOrder::LO;
    wi.build(wilson_config);

    BlockProvider().log_block(ParameterType::SM, "VCKM");
    LOG_INFO(ParameterProvider()({ParameterType::WILSON, "WPARAM_MATCH_SM", 6}));
    LOG_INFO(ParameterProvider()({ParameterType::WILSON, "WPARAM_MATCH_SM", 1}));

    LOG_INFO("C2 =", wi.getM(WGroup::B, WCoef::C2, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C7 =", wi.getM(WGroup::B, WCoef::C7, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C8 =", wi.getM(WGroup::B, WCoef::C8, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C9 =", wi.getM(WGroup::B, WCoef::C9, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C10 =", wi.getM(WGroup::B, WCoef::C10, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C'7 =", wi.getM(WGroup::BPrime, WCoef::CP7, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C'8 =", wi.getM(WGroup::BPrime, WCoef::CP8, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C'9 =", wi.getM(WGroup::BPrime, WCoef::CP9, QCDOrder::LO, ContributionType::SM));
    LOG_INFO("C'10 =", wi.getM(WGroup::BPrime, WCoef::CP10, QCDOrder::LO, ContributionType::SM));

    return 0;
}