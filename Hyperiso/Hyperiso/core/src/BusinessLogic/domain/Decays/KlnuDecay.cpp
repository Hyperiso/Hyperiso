#include "KlnuDecay.h"

void KlnuDecay::load_params() {
    complex_t C_A_K = w_proxy->getFM(WGroup::CC_su, WCoef::C_V2_su, QCDOrder::LO) - w_proxy->getFM(WGroup::CC_su, WCoef::C_V1_su, QCDOrder::LO);
    complex_t C_P_K = w_proxy->getFM(WGroup::CC_su, WCoef::C_S2_su, QCDOrder::LO) - w_proxy->getFM(WGroup::CC_su, WCoef::C_S1_su, QCDOrder::LO);
    cache.calc_K = PlnuCalculator(321, 13, C_A_K, C_P_K, p);

    complex_t C_A_pi = w_proxy->getFM(WGroup::CC_du, WCoef::C_V2_du, QCDOrder::LO) - w_proxy->getFM(WGroup::CC_du, WCoef::C_V1_du, QCDOrder::LO);
    complex_t C_P_pi = w_proxy->getFM(WGroup::CC_du, WCoef::C_S2_du, QCDOrder::LO) - w_proxy->getFM(WGroup::CC_du, WCoef::C_S1_du, QCDOrder::LO);
    cache.calc_pi = PlnuCalculator(211, 13, C_A_pi, C_P_pi, p);

    double f_K = (*p)({ParameterType::FLAVOR, "FCONST", {321, 1}}, DataType::VALUE);
    double f_pi = (*p)({ParameterType::FLAVOR, "FCONST", {211, 1}}, DataType::VALUE);
    double f_K__f_pi = (*p)({ParameterType::FLAVOR, "FCONSTRATIO", {321, 211, 1, 1}}, DataType::VALUE);

    cache.f_corr = std::pow(f_K__f_pi * f_pi / f_K, 2);
    cache.delta_em = (*p)({ParameterType::DECAY, "K_lnu", 1}, DataType::VALUE);
}

double KlnuDecay::BR_K_BR_pi() {
    double BR_K = cache.calc_K.BR_0_SM() * cache.calc_K.R_SM_BSM();
    double BR_pi = cache.calc_pi.BR_0_SM() * cache.calc_pi.R_SM_BSM();
    return BR_K / BR_pi * (1 - cache.delta_em) * cache.f_corr;
}

double KlnuDecay::R_mu23() {
    return std::sqrt(cache.calc_K.R_SM_BSM() / cache.calc_pi.R_SM_BSM());
}

std::vector<ObservableValue> KlnuDecay::compute_observable(Observables obs) {
    double value;
    switch (obs) {
    case Observables::BR_K__MU_NU__BR_PI__MU_NU:   
        value = BR_K_BR_pi();
        break;
    case Observables::R_MU23:   
        value = R_mu23();
        break;
    default:
        LOG_ERROR("IndexError", "Observable", ObservableMapper::str(obs), "doesn't belong to the decay", DecayMapper::str(this->id));
    }

    return {ObservableValue(ObservableMapper::to_id(obs), value)};
}


std::vector<ObservableValue> KlnuDecay::compute_observable(ObservableId obs) {
    return compute_observable(ObservableMapper::enum_of(obs).value());
}