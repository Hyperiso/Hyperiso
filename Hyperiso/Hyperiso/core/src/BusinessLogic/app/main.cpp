#include "Logger.h"
#include <iostream>
#include <cassert>
#include "ObservableInterface.h"
#include "HyperisoMaster.h"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);
    HyperisoMaster hyperiso;
    Config config;
    config.model = Model::SM;
    hyperiso.init("default/lha/testInput.flha", config);
    LOG_INFO("HyperisoMaster initialized");

    auto interface = ObservableInterface();
    interface.add_observable(Observables::BR_B_XS_GAMMA, QCDOrder::NNLO, true);
    LOG_INFO("Observable added to manager");
    
    for (auto& [k, v] : interface.compute_leading_uncertainties(Observables::BR_B_XS_GAMMA, 10)) {
        LOG_INFO(k, ":", v);
    }
    

    LOG_INFO(interface.compute_observable(Observables::BR_B_XS_GAMMA), "+-", interface.compute_uncertainty(Observables::BR_B_XS_GAMMA));
    // LOG_INFO(interface.compute_chi2());

}