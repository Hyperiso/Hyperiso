#include "ObservableInterface.h"
#include "HyperisoMaster.h"
#include "ParameterSetter.h"

int main() {

    HyperisoMaster hyp;
    HyperisoConfig config;
    config.model = Model::SM;
    hyp.init("lha/si_input.flha", config);

    QCDOrder order = QCDOrder::NNLO;
    ObservableInterface oi;

    Decays dec = Decays::B__Xs_gamma;
    // LbLllConfig dec_cfg;
    // dec_cfg.gen = LbLllConfig::Lepton::MU;
    // dec_cfg.ff_src = LbL_FF_Src::DM;
    // dec_cfg.bins = {{1.0, 4.0}};

    // oi.set_decay_config(dec, dec_cfg);
    // oi.add_observable(Observables::TEST, order);
    // oi.compute_observable(Observables::TEST);
    ParameterSetter ps;
    oi.add_observables(dec, order, false);

    oi.add_observable(Observables::BR_BS_MUMU, QCDOrder::LO);
    oi.add_observable(Observables::BR_BS_MUMU_UNTAG, QCDOrder::LO);
    oi.add_observable(Observables::BR_BD_MUMU, QCDOrder::LO);

    auto start = std::chrono::steady_clock::now();
    int i = 0;
    for (double elem = 4.0e-3; elem < 5.0e-3; elem+=0.01e-3) {
        i++;
        ps.mutate({ParameterType::SM, "MASS", 1}, elem);
        // std::cout << "Bs_gamma : " << oi.compute_observable(Observables::BR_B_XS_GAMMA)[0].value << std::endl;
        std::cout << "BR_BS_MUMU : " << oi.compute_observable(Observables::BR_BS_MUMU)[0].value << std::endl;
        std::cout << "BR_BS_MUMU_UNTAG : " << oi.compute_observable(Observables::BR_BS_MUMU_UNTAG)[0].value << std::endl;
        std::cout << "BR_BD_MUMU : " << oi.compute_observable(Observables::BR_BD_MUMU)[0].value << std::endl;
    }
    auto stop  = std::chrono::steady_clock::now();

    auto us = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
    std::cout << "Temps : " << us*1e-6 << " s\n";
    std::cout << "num of steps : " << i << std::endl;

    return 0;
}