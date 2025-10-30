#include "KPinunuDecay.h"

void KPinunuDecay::load_params() {
    ObsParameterProxy p;
    cache.alpha_s_m_Z = p(ParamId{ParameterType::SM, "SMINPUTS", 3});
    cache.sw2 = p(ParamId{ParameterType::SM, "SMINPUTS", {7, 1}});
    cache.m_c_m_c = p(ParamId{ParameterType::SM, "MASS", 4});
    cache.lambda_c = p(ParamId{ParameterType::SM, "VCKM", {1, 0}}) * std::conj(p(ParamId{ParameterType::SM, "VCKM", {1, 1}}));
    cache.lambda_t = p(ParamId{ParameterType::SM, "VCKM", {2, 0}}) * std::conj(p(ParamId{ParameterType::SM, "VCKM", {2, 0}}));
    cache.lambda = p(ParamId{ParameterType::SM, "VCKMIN", 1});

    cache.kappa_L = p(ParamId{ParameterType::SM, "K_pi", 1});
    cache.kappa_p = p(ParamId{ParameterType::SM, "K_pi", 2});
    cache.delta_em = p(ParamId{ParameterType::SM, "K_pi", 3});

    cache.CL = w_proxy->getFR(WGroup::K, WCoef::CK_L, w_config.order);
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
    double kappa = cache.kappa_L * 1e-10 * std::pow(cache.lambda / 0.225, 8);
    return kappa * std::pow(cache.sw2, 2) / std::pow(cache.lambda, 10) * std::pow(std::imag(cache.lambda_t * cache.CL), 2);
}

double KPinunuDecay::BR_p() {
    double kappa = cache.kappa_p * 1e-10 * std::pow(cache.lambda / 0.225, 8);
    double Xc = std::pow(cache.lambda, 4) * P_c();
    return kappa * std::pow(cache.sw2, 2) * (1 + cache.delta_em) / std::pow(cache.lambda, 10) * (
        std::pow(std::imag(cache.lambda_t * cache.CL), 2)
      + std::pow(std::real(cache.lambda_t * cache.CL - cache.lambda_c * Xc / cache.sw2), 2)
    );
}
