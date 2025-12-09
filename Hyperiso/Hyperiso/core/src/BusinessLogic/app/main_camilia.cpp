#include "Logger.h"
#include <iostream>
#include <cassert>
#include "ObservableInterface.h"
#include "HyperisoMaster.h"
#include "config.hpp"
#include "BlockProxy.h"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);
    HyperisoMaster hyp;
    Config config;
    config.model = Model::SM;
    hyp.init("lha/si_input.flha", config);

    QCDOrder order = QCDOrder::NNLO;
    ObservableInterface oi;

    Decays dec = Decays::B__Xs;
    // LbLllConfig dec_cfg;
    // dec_cfg.gen = LbLllConfig::Lepton::MU;
    // dec_cfg.ff_src = LbL_FF_Src::DM;
    // dec_cfg.bins = {{1.0, 4.0}};

    // oi.set_decay_config(dec, dec_cfg);
    // oi.add_observable(Observables::TEST, order);
    // oi.compute_observable(Observables::TEST);

    oi.add_observables(dec, order, false);
    for (auto o : DecayMapper::get_observables(dec)) {
        if (o == Observables::TEST) continue;

        // if (o == Observables::A_FB_B__KSTAR_L_L || o == Observables::F_L_B__KSTAR_L_L) {
            auto obs_values = oi.compute_observable(o);
        std::stringstream ss;
        ss << std::scientific << std::setprecision(3);
        if (obs_values.size() == 1) {
            ss << "= " << obs_values[0].value;
        } else {
            ss << ": ";
            for (auto ov : obs_values) {
                ss << "[" << ov.bin.value().first << ", " << ov.bin.value().second << "] = " << ov.value << ", ";
            }
        }
            
        LOG_INFO(ObservableMapper::str(o), ss.str());   
        // }
    }

    // auto Gamma = oi.compute_observable(Observables::DGAMMA_DQ2_BS__PHI_L_L)[0].value;
    // auto Gamma_bar = oi.compute_observable(Observables::DGAMMA_BAR_DQ2_BS__PHI_L_L)[0].value;

    // std::stringstream ss;
    // ss << std::scientific << std::setprecision(3);
    // ss << "= " << (Gamma + Gamma_bar) / 2;
    // LOG_INFO("BR(Bs > phi mu mu)", ss.str());   

    return 0;
}