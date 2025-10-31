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
    hyp.init("lha/si_input.flha", config);
    // LOG_INFO("HyperisoMaster initialized");

    QCDOrder order = QCDOrder::NNLO;
    ObservableInterface oi;

    for (auto& dec : DecayMapper::get_enum()) {
        if (dec == Decays::M0_Mix) continue;

        // if (dec != Decays::K__l_l) continue;

        LOG_INFO("Adding observables for decay", DecayMapper::str(dec));
        oi.add_observables(dec, order, false);
        
        LOG_INFO("Computing observables for decay", DecayMapper::str(dec));
        for (auto o : DecayMapper::get_observables(dec)) {
            auto obs_values = oi.compute_observable(o);
            std::stringstream ss;
            if (obs_values.size() == 1) {
                ss << "= " << obs_values[0].value;
            } else {
                ss << ": ";
                for (auto ov : obs_values) {
                    ss << "[" << ov.bin.value().first << ", " << ov.bin.value().second << "] = " << ov.value << ", ";
                }
            }
                
            LOG_INFO(ObservableMapper::str(o), ss.str());
        }
    }

    return 0;
}