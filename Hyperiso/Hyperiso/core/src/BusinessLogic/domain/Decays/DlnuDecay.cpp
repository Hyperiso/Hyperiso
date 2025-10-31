#include "DlnuDecay.h"

void DlnuDecay::load_params() {
    complex_t C_A = w_proxy->getFM(WGroup::CC_cd, WCoef::C_V2_cd, QCDOrder::LO) - w_proxy->getFM(WGroup::CC_cd, WCoef::C_V1_cd, QCDOrder::LO);
    complex_t C_P = w_proxy->getFM(WGroup::CC_cd, WCoef::C_S2_cd, QCDOrder::LO) - w_proxy->getFM(WGroup::CC_cd, WCoef::C_S1_cd, QCDOrder::LO);
    cache.calc = PlnuCalculator(411, 13, C_A, C_P);
}

double DlnuDecay::BR() {
    return cache.calc.BR_0_SM() * cache.calc.R_SM_BSM();
}

std::vector<ObservableValue> DlnuDecay::compute_observable(Observables obs) {
    double value;
    switch (obs) {
    case Observables::BR_D__MU_NU:   
        value = BR();
        break;
    default:
        LOG_ERROR("IndexError", "Observable", ObservableMapper::str(obs), "doesn't belong to the decay", DecayMapper::str(this->id));
    }

    return {ObservableValue(ObservableMapper::to_id(obs), value)};
}


std::vector<ObservableValue> DlnuDecay::compute_observable(ObservableId obs) {
    return compute_observable(ObservableMapper::enum_of(obs).value());
}