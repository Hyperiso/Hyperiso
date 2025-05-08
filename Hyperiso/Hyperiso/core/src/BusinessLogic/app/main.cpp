#include "Logger.h"
#include <iostream>
#include <cassert>
#include "ObservableInterface.h"
#include "HyperisoMaster.h"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::TRACE);

    LOG_INFO("Starting HyperisoMaster test");
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);
    HyperisoMaster hyperiso;
    Config config;
    config.model = Model::SM;
    hyperiso.init("default/lha/testInput.flha", config);
    LOG_INFO("HyperisoMaster initialized");

    auto interface = ObservableInterface();

    // interface.add_observable(Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA, QCDOrder::NNLO, true);

    // auto us = interface.compute_leading_uncertainties(Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA, 5);
    // for (auto& [pid, u] : us) {
    //     LOG_INFO(pid, ":", u);
    // }

    // LOG_INFO(interface.compute_observable(Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA));

    for (Observables obs : ObservableMapper::get_enum()) {
        interface.add_observable(obs, QCDOrder::NNLO, true);
    }
    // interface.add_observable(Observables::BR_BS_MUMU_UNTAG, QCDOrder::NNLO, true);
    LOG_INFO("All observables added");
    // LOG_INFO(ObservableMapper::str(Observables::BR_BS_MUMU_UNTAG), "=", interface.compute_observable(Observables::BR_BS_MUMU_UNTAG), "+-", interface.compute_uncertainty(Observables::BR_BS_MUMU_UNTAG));
    for (Observables obs : ObservableMapper::get_enum()) {
        LOG_INFO(ObservableMapper::str(obs), "=", interface.compute_observable(obs), "+-", interface.compute_uncertainty(obs));
    }
    
    LOG_INFO(interface.compute_chi2());
}