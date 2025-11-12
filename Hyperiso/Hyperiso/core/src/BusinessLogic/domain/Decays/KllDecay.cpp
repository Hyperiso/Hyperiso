#include "KllDecay.h"

void KllDecay::load_params() {
    ObsParameterProxy p;
    cache.G_F = p(ParamId{ParameterType::SM, "SMINPUTS", 2});
    cache.alpha_em = p(ParamId{ParameterType::SM, "EW", {1, 1}});
    cache.alpha_em_0 = p(ParamId{ParameterType::SM, "EW", {1, 4}});
    cache.alpha_s_m_Z = p(ParamId{ParameterType::SM, "SMINPUTS", 3});
    cache.sw2 = p(ParamId{ParameterType::SM, "SMINPUTS", {7, 1}});
    cache.m_l = p(ParamId{ParameterType::SM, "MASS", 9 + 2 * cfg.gen});
    cache.m_c_m_c = p(ParamId{ParameterType::SM, "MASS", 4});
    cache.m_pi = p(ParamId{ParameterType::FLAVOR, "FMASS", 111});
    cache.m_rho = p(ParamId{ParameterType::FLAVOR, "FMASS", 113});
    cache.m_K = p(ParamId{ParameterType::FLAVOR, "FMASS", 311});
    cache.f_K = p(ParamId{ParameterType::FLAVOR, "FCONST", {311, 1}});
    cache.tau_L = p(ParamId{ParameterType::FLAVOR, "FLIFE", 130}) / HBAR;
    cache.tau_S = p(ParamId{ParameterType::FLAVOR, "FLIFE", 310}) / HBAR;
    cache.lambda_c = p(ParamId{ParameterType::SM, "VCKM", {1, 0}}) * std::conj(p(ParamId{ParameterType::SM, "VCKM", {1, 1}}));
    cache.lambda_t = p(ParamId{ParameterType::SM, "VCKM", {2, 0}}) * std::conj(p(ParamId{ParameterType::SM, "VCKM", {2, 0}}));
    cache.lambda = p(ParamId{ParameterType::SM, "VCKMIN", 1});
    cache.x = cache.m_l / cache.m_K;
    cache.r_chi = cache.m_K / (p(ParamId{ParameterType::SM, "MASS", 1}) + p(ParamId{ParameterType::SM, "MASS", 3}));
    cache.beta = std::sqrt(1. - 4. * std::pow(cache.x, 2));
    cache.BR_KL_gg_exp = p(ParamId{ParameterType::DECAY, "K_ll", 1});
    cache.BR_KS_gg_exp = p(ParamId{ParameterType::DECAY, "K_ll", 2});
    cache.alpha_exp = p(ParamId{ParameterType::DECAY, "K_ll", 3});
    cache.delta_lambda = p(ParamId{ParameterType::DECAY, "K_ll", 4});
    cache.mu_b = 5 + 1e-10; // TODO : promote as decay param ?

    cache.C10 = w_proxy->getFR(WGroup::K, WCoef::CK10, w_config.order) - w_proxy->getFR(WGroup::K, WCoef::CPK10, w_config.order);
    cache.CQ1 = w_proxy->getFR(WGroup::K, WCoef::CKQ1, w_config.order) - w_proxy->getFR(WGroup::K, WCoef::CPKQ1, w_config.order);
    cache.CQ2 = w_proxy->getFR(WGroup::K, WCoef::CKQ2, w_config.order) - w_proxy->getFR(WGroup::K, WCoef::CPKQ2, w_config.order);

    printf("C10m = %.5e\n", cache.C10);
	printf("CQ1m = %.5e\n", cache.CQ1);
	printf("CQ2m = %.5e\n", cache.CQ2);
}

complex_t KllDecay::C_gg(double b) {
    complex_t x = (b - 1) / (b + 1);
    return 1 / b * (CLi2(x) + PI2 / 3 + 0.25 * std::pow(std::log(x), 2));
}

complex_t KllDecay::N_L() {
    double sign = this->cfg.N_L_sign / std::abs(this->cfg.N_L_sign);
    double L_mu_rho = std::log(cache.m_l / cache.m_rho);
    double chi_gg = -3 * cache.alpha_exp + std::pow(cache.m_K / cache.m_rho, 2) * (1. / 3 + cache.alpha_exp * (1 + 12 * L_mu_rho) / 18) - cache.delta_lambda;
    complex_t chi = chi_gg - 2.5 + 3.0 * L_mu_rho + C_gg(cache.beta);
    double N0 = 4 * cache.alpha_em_0 * cache.m_l / (PI * cache.f_K * std::pow(cache.m_K, 2)) 
                    * std::sqrt(2 * PI * cache.BR_KL_gg_exp / (cache.m_K * cache.tau_L));

    return sign * N0 * chi;
}

complex_t KllDecay::N_S() {
    complex_t I_S_mu = {-2.821, 1.216}; // TODO : perform integration over 3D space (VEGAS)
    complex_t I_S_e = {2.11, -40.41}; // TODO : perform integration over 3D space (VEGAS)
    complex_t I_S = cfg.gen == 1 ? I_S_e : I_S_mu;
    double r_pi = cache.m_pi / cache.m_K;
    double N0 = 2 * cache.alpha_em_0 * cache.m_l / (PI * cache.f_K * std::pow(cache.m_K, 2) * std::abs(KP::H(0.0, r_pi)))
                    * std::sqrt(2 * PI * cache.BR_KS_gg_exp / (cache.m_K * cache.tau_S));
    return N0 * I_S;
}

double KllDecay::P_c() {
    double kappa100 = -0.4499;
	double kappa010 = -5.9221;
	double kappa001 = 0.0114;
	double kappa110 = 3.9957;
	double kappa011 = -0.0658;
	double kappa200 = -0.1767;
	double kappa020 = 16.4465;
    double K = 0.9622;

	double L_c = log(cache.m_c_m_c / 1.3);
	double L_a = log(cache.alpha_s_m_Z / 0.1187);
	double L_b = log(cache.mu_b / 5.);
	
	return pow(0.225 / cache.lambda, 4.) * (0.1198 * pow(cache.m_c_m_c / 1.3, 2.3595) * pow(cache.alpha_s_m_Z/ 0.1187, 6.6055) * (
        K
      + kappa100 * L_c + kappa010 * L_a + kappa001 * L_b 
      + kappa110 * L_c * L_a + kappa011 * L_a * L_b + kappa200 * L_c * L_c + kappa020 * L_a * L_a)
    );
}

double KllDecay::BR_L() {
    return cache.tau_L * std::pow(cache.f_K * cache.G_F * cache.alpha_em, 2) * std::pow(cache.m_K, 3) * cache.beta / (32 * PI3) * (
        std::pow(cache.beta * std::abs(
            cache.r_chi * std::imag(cache.lambda_t * cache.CQ1)
        ), 2)
      + std::pow(std::abs(
            RT2 * PI / (cache.G_F * cache.alpha_em) * N_L()
          - 2 * cache.x * std::real(-cache.lambda_c * std::pow(cache.lambda, 4) * P_c() / cache.sw2 + cache.lambda_t * cache.C10)
          - cache.r_chi * std::real(cache.lambda_t * cache.CQ2)
        ), 2)
    );
}

double KllDecay::BR_S() {
    return cache.tau_S * std::pow(cache.f_K * cache.G_F * cache.alpha_em, 2) * std::pow(cache.m_K, 3) * cache.beta / (32 * PI3) * (
        std::pow(cache.beta * std::abs(
            RT2 * PI / (cache.G_F * cache.alpha_em) * N_S()
          - cache.r_chi * std::real(cache.lambda_t * cache.CQ1)
        ), 2)
      + std::pow(std::abs(
            2 * cache.x * std::imag(-cache.lambda_c * std::pow(cache.lambda, 4) * P_c() / cache.sw2 + cache.lambda_t * cache.C10)
          + cache.r_chi * std::imag(cache.lambda_t * cache.CQ2)
        ), 2)
    );
}

std::vector<ObservableValue> KllDecay::compute_observable(Observables obs) {
    double value;
    switch (obs) {
    case Observables::BR_KL__MU_MU:   
        value = BR_L();
        break;
    case Observables::BR_KS__MU_MU:   
        value = BR_S();
        break;
    default:
        LOG_ERROR("IndexError", "Observable", ObservableMapper::str(obs), "doesn't belong to the decay", DecayMapper::str(this->id));
    }

    return {ObservableValue(ObservableMapper::to_id(obs), value)};
}

std::vector<ObservableValue> KllDecay::compute_observable(ObservableId obs) {
    return compute_observable(ObservableMapper::enum_of(obs).value());
}
