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
    LOG_INFO("HyperisoMaster initialized");

    QCDOrder order = QCDOrder::LO;

    std::vector<Observables> obss = DecayMapper::get_observables(Decays::B__Kstar_l_l);
    LOG_INFO(obss.size());

    ObservableInterface oi;

    LOG_INFO("Interface created");

    for (auto o : obss) {
        oi.add_observable(o, order, false);
    }

    for (auto o : obss) {
        LOG_INFO(ObservableMapper::str(o), "=", oi.compute_observable(o), "+-", oi.compute_uncertainty(o));
    }
    
    // auto print_leading = [&oi] (Observables o, size_t n) {
    //     LOG_INFO("---------- Leading uncertainties for", ObservableMapper::str(o));
    //     for (const auto& [pid, u] : oi.compute_leading_uncertainties(o, n)) {
    //         LOG_INFO("\t-", pid, ":", u.real());
    //     }
    // };

    // print_leading(Observables::BR_B__Xs_mu_mu__LOW_Q2, 5);
    // print_leading(Observables::BR_B__Xs_mu_mu__HIGH_Q2, 5);

    return 0;
}