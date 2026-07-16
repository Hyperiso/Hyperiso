#include "DslnuDecay.h"

#include <stdexcept>

void DslnuDecay::load_params() {
    complex_t C_A = w_proxy->getFM(WGroup::CC_cs, WCoef::C_V2_cs, QCDOrder::LO) - w_proxy->getFM(WGroup::CC_cs, WCoef::C_V1_cs, QCDOrder::LO);
    complex_t C_P = w_proxy->getFM(WGroup::CC_cs, WCoef::C_S2_cs, QCDOrder::LO) - w_proxy->getFM(WGroup::CC_cs, WCoef::C_S1_cs, QCDOrder::LO);
    cache.calc_mu = PlnuCalculator(431, 13, C_A, C_P, p);
    cache.calc_tau = PlnuCalculator(431, 15, C_A, C_P, p);
}

double DslnuDecay::BR(int gen) {
    double BR = 0.0;
    switch (gen)
    {
    case 2:
        BR = cache.calc_mu.BR_0_SM() * cache.calc_mu.R_SM_BSM();
        break;
    case 3:
        BR = cache.calc_tau.BR_0_SM() * cache.calc_tau.R_SM_BSM();
        break;
    default:
        throw std::invalid_argument("DslnuDecay supports only muon (2) and tau (3) generations");
    }
    return BR;
}

std::vector<ObservableValue> DslnuDecay::compute_observable(Observables obs) {
    double value = 0.0;
    switch (obs) {
    case Observables::BR_DS__MU_NU:   
        value = BR(2);
        break;
    case Observables::BR_DS__TAU_NU:   
        value = BR(3);
        break;
    default:
        throw std::invalid_argument(
            "Observable " + ObservableMapper::str(obs) + " does not belong to " + DecayMapper::str(this->id)
        );
    }

    return {ObservableValue(ObservableMapper::to_id(obs), value)};
}


std::vector<ObservableValue> DslnuDecay::compute_observable(ObservableId obs) {
    return compute_observable(ObservableMapper::enum_of(obs).value());
}