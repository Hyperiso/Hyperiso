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
    config.mty_model_name = "THDM";
    config.mty_model_path = project_assets_root.data() + std::string("input_files/marty_model/thdm.h");
    hyp.init("lha/testInput.flha", config);
    
    LOG_INFO("HyperisoMaster initialized");

    ParameterProvider pp;
    ParameterSetter ps;
    WilsonInterface wi;

    WilsonBuildConfig wilson_config;
    wilson_config.groups = {WGroup::B, WGroup::BPrime};
    wilson_config.matching_scale = 81;
    wilson_config.hadronic_scale = pp({ParameterType::SM, "QCD", {5, 3}}) / 2;
    // wilson_config.hadronic_scale = 42;
    wilson_config.order = QCDOrder::NNLO;
    wi.build(wilson_config);

    BlockProxy().log_block(ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC));

    ObservableInterface oi;
    oi.add_observable(Observables::BR_B__Xs_l_l__LOW_Q2, wilson_config.order, false);
    oi.add_observable(Observables::BR_B__Xs_l_l__HIGH_Q2, wilson_config.order, false);

    LOG_INFO("BR(B > X_s mu mu) (q²=[1, 6] GeV²) =", oi.compute_observable(Observables::BR_B__Xs_l_l__LOW_Q2)/*, "+-", oi.compute_uncertainty(Observables::BR_B__Xs_l_l__LOW_Q2)*/);
    LOG_INFO("BR(B > X_s mu mu) (q²>14.2 GeV²) =", oi.compute_observable(Observables::BR_B__Xs_l_l__HIGH_Q2)/*, "+-", oi.compute_uncertainty(Observables::BR_B__Xs_l_l__HIGH_Q2)*/);

    return 0;
}