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
    hyp.init("lha/testInput.flha", config);
    LOG_INFO("HyperisoMaster initialized");

    QCDOrder order = QCDOrder::NNLO;

    ObservableInterface oi;
    oi.add_observable(Observables::BR_B__Xs_l_l__LOW_Q2, order, true);
    oi.add_observable(Observables::BR_B__Xs_l_l__HIGH_Q2, order, true);

    LOG_INFO("BR(B > X_s mu mu) (q²=[1, 6] GeV²) =", oi.compute_observable(Observables::BR_B__Xs_l_l__LOW_Q2), "+-", oi.compute_uncertainty(Observables::BR_B__Xs_l_l__LOW_Q2));
    LOG_INFO("BR(B > X_s mu mu) (q²>14.2 GeV²) =", oi.compute_observable(Observables::BR_B__Xs_l_l__HIGH_Q2), "+-", oi.compute_uncertainty(Observables::BR_B__Xs_l_l__HIGH_Q2));

    auto print_leading = [&oi] (Observables o, size_t n) {
        LOG_INFO("---------- Leading uncertainties for", ObservableMapper::str(o));
        for (const auto& [pid, u] : oi.compute_leading_uncertainties(o, n)) {
            LOG_INFO("\t-", pid, ":", u.real());
        }
    };

    print_leading(Observables::BR_B__Xs_l_l__LOW_Q2, 5);
    print_leading(Observables::BR_B__Xs_l_l__HIGH_Q2, 5);

    return 0;
}