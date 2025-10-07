#include "BlnuDecay.h"

void BlnuDecay::load_params() {
    ObsParameterProxy p;
    cache.G_F = p(ParamId{ParameterType::SM, "SMINPUTS", 2});
    cache.m_tau = p(ParamId{ParameterType::SM, "MASS", 15});
    cache.m_b = p(ParamId{ParameterType::SM, "QCD", {5, 1}});
    cache.m_B = p(ParamId{ParameterType::FLAVOR, "FMASS", 521});
    cache.f_B = p(ParamId{ParameterType::FLAVOR, "FCONST", {521, 1}});
    cache.tau_B = p(ParamId{ParameterType::FLAVOR, "FLIFE", 521});
    cache.V_ub_2 = std::pow(std::abs(p(ParamId{ParameterType::SM, "VCKM", {0, 2}})), 2);
    cache.C_V = w_proxy->getFM(WGroup::BCC, WCoef::C_V1, QCDOrder::LO);
    cache.C_S = w_proxy->getFM(WGroup::BCC, WCoef::C_S1, QCDOrder::LO);
}

scalar_t BlnuDecay::R() {
    return std::pow(std::abs(cache.C_V + std::pow(cache.m_B, 2) * cache.C_S / (cache.m_b * cache.m_tau)), 2);
}

double BlnuDecay::BR() {
    double beta = 1 - std::pow(cache.m_tau / cache.m_B, 2);
    double pref = std::pow(cache.G_F * cache.f_B * cache.m_tau * beta, 2) * cache.tau_B * cache.m_B * cache.V_ub_2 / (8. * PI * HBAR);
    return pref * R();
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