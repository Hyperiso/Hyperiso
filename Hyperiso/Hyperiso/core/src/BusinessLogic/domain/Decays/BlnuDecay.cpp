#include "BlnuDecay.h"

void BlnuDecay::load_params() {
    complex_t C_A = w_proxy->getFM(WGroup::CC_bu, WCoef::C_V2_bu, QCDOrder::LO) - w_proxy->getFM(WGroup::CC_bu, WCoef::C_V1_bu, QCDOrder::LO);
    complex_t C_P = w_proxy->getFM(WGroup::CC_bu, WCoef::C_S2_bu, QCDOrder::LO) - w_proxy->getFM(WGroup::CC_bu, WCoef::C_S1_bu, QCDOrder::LO);
    cache.calc = PlnuCalculator(521, 15, C_A, C_P);
}

scalar_t BlnuDecay::R() {
    return cache.calc.R_SM_BSM();
}

double BlnuDecay::BR() {
    return cache.calc.BR_0_SM() * cache.calc.R_SM_BSM();
}

std::vector<ObservableValue> BlnuDecay::compute_observable(Observables obs) {
    double value;
    switch (obs) {
    case Observables::R_TAU_NU:   
        value = R();
        break;
    case Observables::BR_BU_TAU_NU:   
        value = BR();
        break;
    default:
        LOG_ERROR("IndexError", "Observable", ObservableMapper::str(obs), "doesn't belong to the decay", DecayMapper::str(this->id));
    }

    return {ObservableValue(ObservableMapper::to_id(obs), value)};
}


std::vector<ObservableValue> BlnuDecay::compute_observable(ObservableId obs) {
    return compute_observable(ObservableMapper::enum_of(obs).value());
}