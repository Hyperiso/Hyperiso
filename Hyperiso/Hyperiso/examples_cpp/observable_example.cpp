#include "Logger.h"
#include <iostream>
#include <cassert>
#include "ObservableInterface.h"
#include "HyperisoMaster.h"
#include "config.hpp"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);
    HyperisoMaster hyp; // Create the interface for hyperiso (reading/writing in the lha and all the options we want to use)
    Config config; // Config struct where we can put all the options we want for Hyperiso (general options)
    config.model = Model::SM; // The model we want to use, SM by default. If not THDM or SUSY, MARTY is needed.
    hyp.init("lha/si_input.flha", config); // Initialize program manager with LHA file and the config. Search in the Assets directory if relative path.
    LOG_INFO("HyperisoMaster initialized");

    QCDOrder order = QCDOrder::NNLO;

    ObservableInterface oi; // Initialize interface for observables calculation.

    oi.add_observable(Observables::BR_B__Xs_mu_mu__LOW_Q2, order, true); // add an observable to the interface, at a specific order. The last argument is used to calculate uncertainty with  all the parameters (if true).
    oi.add_observable(Observables::BR_B__Xs_mu_mu__HIGH_Q2, order, false);
    oi.add_observable(Observables::BR_B__Xs_tau_tau__HIGH_Q2, order, false);

    LOG_INFO("BR(B > X_s mu mu) (q²=[1, 6] GeV²) =", oi.compute_observable(Observables::BR_B__Xs_mu_mu__LOW_Q2), "+-", oi.compute_uncertainty(Observables::BR_B__Xs_mu_mu__LOW_Q2));
    //compute_observable calculate and retrieve the value of the observable. compute_uncertainty does the same with the uncertainty.

    
    LOG_INFO("BR(B > X_s mu mu) (q²>14.2 GeV²) =", oi.compute_observable(Observables::BR_B__Xs_mu_mu__HIGH_Q2), "+-", oi.compute_uncertainty(Observables::BR_B__Xs_mu_mu__HIGH_Q2));
    LOG_INFO("BR(B > X_s tau tau) (q²>14.2 GeV²) =", oi.compute_observable(Observables::BR_B__Xs_tau_tau__HIGH_Q2), "+-", oi.compute_uncertainty(Observables::BR_B__Xs_tau_tau__HIGH_Q2));

    auto print_leading = [&oi] (Observables o, size_t n) {
        LOG_INFO("---------- Leading uncertainties for", ObservableMapper::str(o));
        for (const auto& [pid, u] : oi.compute_leading_uncertainties(o, n)) {
            LOG_INFO("\t-", pid, ":", u.real());
        }
    };

    print_leading(Observables::BR_B__Xs_mu_mu__LOW_Q2, 5);
    print_leading(Observables::BR_B__Xs_mu_mu__HIGH_Q2, 5);

    return 0;
}