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
    config.model = Model::THDM;
    config.flags[ExternalFlag::USE_MARTY] = false;
    config.mty_model_name = "THDM";
    config.mty_model_path = project_assets_root.data() + std::string("input_files/marty_model/thdm.h");
    hyp.init("lha/testinput_thdm.lha", config);
    
    LOG_INFO("HyperisoMaster initialized");

    ParameterProvider pp;
    ParameterSetter ps;
    WilsonInterface wi;

    WilsonBuildConfig wilson_config;
    wilson_config.groups = {WGroup::B};
    wilson_config.matching_scale = 81;
    wilson_config.hadronic_scale = pp({ParameterType::SM, "QCD", {5, 3}}) / 2;
    // wilson_config.hadronic_scale = 42;
    wilson_config.order = QCDOrder::LO;
    wi.build(wilson_config);

    BlockProxy bp;

    LOG_INFO("------------- Matching blocks -------------");

    bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::MATCHING));
    // bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::BPrime, ScaleType::MATCHING));
    // bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::BScalar, ScaleType::MATCHING));
    // bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::BCC, ScaleType::MATCHING));

    LOG_INFO("------------- Hadronic blocks -------------");

    bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC));
    // bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC, WilsonBasis::B_TRADITIONAL));
    // bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::BPrime, ScaleType::HADRONIC));
    // bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::BScalar, ScaleType::HADRONIC));
    // bp.log_block(ParameterType::WILSON, GroupMapper::str(WGroup::BCC, ScaleType::HADRONIC));
    ObservableInterface oi;
    oi.add_observable(Observables::BR_B_XS_GAMMA, wilson_config.order, true);

    std::ofstream ofs;
    ofs.open("B_s_gamma_SI.csv");
    ofs << "M_Hp,BR,u(BR)\n";

    double log_m_min = 2;
    double log_m_max = log10(5e3);
    size_t n = 100;
    double dl = (log_m_max - log_m_min) / n;

    double lm = log_m_min;
    for (size_t k = 0; k < n; k++) {
        double m = std::pow(10, lm);
        ps.mutate({ParameterType::BSM, "MASS", 37}, m);
        ofs << m << ","
            << oi.compute_observable(Observables::BR_B_XS_GAMMA).real() << ","
            << oi.compute_uncertainty(Observables::BR_B_XS_GAMMA).real() << "\n";
        lm += dl;
    }   

    return 0;
}