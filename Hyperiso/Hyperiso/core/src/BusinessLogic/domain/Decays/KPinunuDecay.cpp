#include "KPinunuDecay.h"

namespace {
constexpr std::array<WCoef, 3> KNU_LEFT = {
    WCoef::CKNU_L_EE, WCoef::CKNU_L_MUMU, WCoef::CKNU_L_TAUTAU
};
constexpr std::array<WCoef, 3> KNU_RIGHT = {
    WCoef::CKNU_R_EE, WCoef::CKNU_R_MUMU, WCoef::CKNU_R_TAUTAU
};
}

void KPinunuDecay::load_params() {
    cache.alpha_s_m_Z = (*p)(ParamId{ParameterType::SM, "SMINPUTS", 3}, DataType::VALUE);
    cache.sw2 = (*p)(ParamId{ParameterType::SM, "SMINPUTS", {7, 1}}, DataType::VALUE);
    cache.m_c_m_c = (*p)(ParamId{ParameterType::SM, "MASS", 4}, DataType::VALUE);
    cache.lambda_c = (*p)(ParamId{ParameterType::SM, "VCKM", {1, 0}}, DataType::VALUE)
                   * std::conj((*p)(ParamId{ParameterType::SM, "VCKM", {1, 1}}, DataType::VALUE));
    cache.lambda_t = (*p)(ParamId{ParameterType::SM, "VCKM", {2, 0}}, DataType::VALUE)
                   * std::conj((*p)(ParamId{ParameterType::SM, "VCKM", {2, 1}}, DataType::VALUE));
    cache.lambda = (*p)(ParamId{ParameterType::SM, "VCKMIN", 1}, DataType::VALUE);

    cache.kappa_L = (*p)(ParamId{ParameterType::DECAY, "K_pi", 1}, DataType::VALUE);
    cache.kappa_p = (*p)(ParamId{ParameterType::DECAY, "K_pi", 2}, DataType::VALUE);
    cache.delta_em = (*p)(ParamId{ParameterType::DECAY, "K_pi", 3}, DataType::VALUE);

    for (std::size_t i = 0; i < KNU_LEFT.size(); ++i) {
        cache.C_L[i] = w_proxy->getFR(WGroup::KNuNu, KNU_LEFT[i], w_config.order,
                                      ContributionType::TOTAL);
        cache.C_R[i] = w_proxy->getFR(WGroup::KNuNu, KNU_RIGHT[i], w_config.order,
                                      ContributionType::TOTAL);
    }
}

double KPinunuDecay::P_c() {
    double kappa10 = 1.6624;
    double kappa01 = -2.3537;
    double kappa11 = -1.5862;
    double kappa20 = 1.5036;
    double kappa02 = -4.3477;
    double Lmc = log(cache.m_c_m_c/1.3);
    double Lalphas = log(cache.alpha_s_m_Z/0.1176);

    return pow(0.2255/cache.lambda,4.)*0.38049*pow(cache.m_c_m_c/1.3,0.5081)*pow(cache.alpha_s_m_Z/0.1176,1.0192)*(1. + (Lalphas*kappa01+Lalphas*Lalphas*kappa02+Lmc*kappa10 + Lmc*Lalphas*kappa11 + Lmc*Lmc*kappa20));
}

double KPinunuDecay::BR_L() {
    const double kappa = cache.kappa_L * 1e-10 * std::pow(cache.lambda / 0.225, 8);
    double top_sum = 0.0;
    for (std::size_t i = 0; i < cache.C_L.size(); ++i) {
        const complex_t C = cache.C_L[i] + cache.C_R[i];
        top_sum += std::pow(std::imag(cache.lambda_t * C), 2);
    }

    // kappa_L is normalized to the conventional three-flavour SM expression.
    return kappa * std::pow(cache.sw2, 2) / std::pow(cache.lambda, 10)
         * (top_sum / 3.0);
}

double KPinunuDecay::BR_p() {
    const double kappa = cache.kappa_p * 1e-10 * std::pow(cache.lambda / 0.225, 8);
    const double Xc = std::pow(cache.lambda, 4) * P_c();
    // With lambda_c = V_cd V_cs^* and lambda_t = V_td V_ts^*, the SM CKM
    // convention used by HyperIso gives negative real parts for both terms.
    // The physical charged-mode amplitude is therefore top + charm here.
    const complex_t charm = cache.lambda_c * Xc / cache.sw2;

    double flavour_sum = 0.0;
    for (std::size_t i = 0; i < cache.C_L.size(); ++i) {
        const complex_t C = cache.C_L[i] + cache.C_R[i];
        const complex_t top = cache.lambda_t * C;
        flavour_sum += std::pow(std::imag(top), 2)
                     + std::pow(std::real(top + charm), 2);
    }

    return kappa * std::pow(cache.sw2, 2) * (1 + cache.delta_em)
         / std::pow(cache.lambda, 10) * (flavour_sum / 3.0);
}

std::vector<ObservableValue> KPinunuDecay::compute_observable(Observables obs) {
    double value;
    switch (obs) {
    case Observables::BR_K__PI_NU_NU:
        value = BR_p();
        break;
    case Observables::BR_KL__PI0_NU_NU:
        value = BR_L();
        break;
    default:
        LOG_ERROR("IndexError", "Observable", ObservableMapper::str(obs), "doesn't belong to the decay", DecayMapper::str(this->id));
    }

    return {ObservableValue(ObservableMapper::to_id(obs), value)};
}

std::vector<ObservableValue> KPinunuDecay::compute_observable(ObservableId obs) {
    return compute_observable(ObservableMapper::enum_of(obs).value());
}
