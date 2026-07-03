#include "BllDecay.h"


void BllDecay::load_params() {
    cache.G_F = (*p)(ParamId{ParameterType::SM, "SMINPUTS", 2}, DataType::VALUE);
    cache.alpha_em = (*p)(ParamId{ParameterType::SM, "EW", {1, 2}}, DataType::VALUE);
    cache.m_e = (*p)(ParamId{ParameterType::SM, "MASS", 11}, DataType::VALUE);
    cache.m_mu = (*p)(ParamId{ParameterType::SM, "MASS", 13}, DataType::VALUE);
    cache.m_tau = (*p)(ParamId{ParameterType::SM, "MASS", 15}, DataType::VALUE);
    cache.m_Bd = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", 511}, DataType::VALUE);
    cache.m_Bs = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", 531}, DataType::VALUE);
    cache.f_Bd = (*p)(ParamId{ParameterType::FLAVOR, "FCONST", {511, 1}}, DataType::VALUE);
    cache.f_Bs = (*p)(ParamId{ParameterType::FLAVOR, "FCONST", {531, 1}}, DataType::VALUE);
    cache.tau_Bd = (*p)(ParamId{ParameterType::FLAVOR, "FLIFE", 511}, DataType::VALUE) / HBAR;
    cache.tau_Bs = (*p)(ParamId{ParameterType::FLAVOR, "FLIFE", 531}, DataType::VALUE) / HBAR;
    cache.lambda_d = (*p)(ParamId{ParameterType::SM, "VCKM", {2, 2}}, DataType::VALUE) * std::conj((*p)(ParamId{ParameterType::SM, "VCKM", {2, 0}}, DataType::VALUE));
    cache.lambda_s = (*p)(ParamId{ParameterType::SM, "VCKM", {2, 2}}, DataType::VALUE) * std::conj((*p)(ParamId{ParameterType::SM, "VCKM", {2, 1}}, DataType::VALUE));
    cache.ys = (*p)(ParamId{ParameterType::DECAY, "B_ll", 1}, DataType::VALUE);
    cache.eta_BBS = (*p)(ParamId{ParameterType::DECAY, "B_ll", 2}, DataType::VALUE);
    cache.r_d = cache.m_Bd / ((*p)(ParamId{ParameterType::SM, "QCD", {5, 2}}, DataType::VALUE) + (*p)(ParamId{ParameterType::SM, "MASS", 1}, DataType::VALUE));
    cache.r_s = cache.m_Bs / ((*p)(ParamId{ParameterType::SM, "QCD", {5, 2}}, DataType::VALUE) + (*p)(ParamId{ParameterType::SM, "MASS", 3}, DataType::VALUE));

    cache.C10_SM = w_proxy->getFR(WGroup::B, WCoef::C10, w_config.order, ContributionType::SM);
    cache.C10 = w_proxy->getFR(WGroup::B, WCoef::C10, w_config.order);
    cache.CQ1 = w_proxy->getFR(WGroup::BScalar, WCoef::CQ1, w_config.order); // ASK : Why 0 in SI ?
    cache.CQ2 = w_proxy->getFR(WGroup::BScalar, WCoef::CQ2, w_config.order);
    // cache.CQ1 = 0;
    // cache.CQ2 = 0;
    cache.C10_m = cache.C10 - w_proxy->getFR(WGroup::BPrime, WCoef::CP10, w_config.order);
    cache.CQ1_m = cache.CQ1 - w_proxy->getFR(WGroup::BPrime, WCoef::CPQ1, w_config.order);
    cache.CQ2_m = cache.CQ2 - w_proxy->getFR(WGroup::BPrime, WCoef::CPQ2, w_config.order);

    printf("f_Bs = %.5e\n", cache.f_Bs);
    printf("y_s = %.5e\n", cache.ys);
    printf("eta_BBs = %.5e\n", cache.eta_BBS);
    printf("r_s = %.5e\n", cache.r_s);
    printf("C10 = %.4e + %.4e I \n", std::real(cache.C10), std::imag(cache.C10));
    printf("Cp10 = %.4e + %.4e I \n", std::real(w_proxy->getFR(WGroup::BPrime, WCoef::CP10, w_config.order)), std::imag(w_proxy->getFR(WGroup::BPrime, WCoef::CP10, w_config.order)));
	printf("CQ1 = %.4e + %.4e I \n", std::real(cache.CQ1), std::imag(cache.CQ1));
	printf("CpQ1 = %.4e + %.4e I \n", std::real(w_proxy->getFR(WGroup::BPrime, WCoef::CPQ1, w_config.order)), std::imag(w_proxy->getFR(WGroup::BPrime, WCoef::CPQ1, w_config.order)));
	printf("CQ2 = %.4e + %.4e I \n", std::real(cache.CQ2), std::imag(cache.CQ2));
	printf("CpQ2 = %.4e + %.4e I \n", std::real(w_proxy->getFR(WGroup::BPrime, WCoef::CPQ2, w_config.order)), std::imag(w_proxy->getFR(WGroup::BPrime, WCoef::CPQ2, w_config.order)));
}

double BllDecay::BR_avg_Bq_ll(int q, int gen) {
    if (q != 1 && q != 3) LOG_ERROR("ValueError", "In Bq > l l, q can only be d (1) or s (3), found", q);
    
    double m_B, f_B, tau_B, r_q;
    complex_t lambda_q;
    if (q == 1) {
        m_B = cache.m_Bd;
        f_B = cache.f_Bd;
        tau_B = cache.tau_Bd;
        lambda_q = cache.lambda_d;
        r_q = cache.r_d;
    } else {
        m_B = cache.m_Bs;
        f_B = cache.f_Bs;
        tau_B = cache.tau_Bs;
        lambda_q = cache.lambda_s;
        r_q = cache.r_s;
    }

    double m_l = gen == 1 ? cache.m_e : gen == 2 ? cache.m_mu : cache.m_tau;
    double x_l = m_l / m_B;
    double beta_l = std::sqrt(1. - 4. * std::pow(x_l, 2));

    double pref = std::pow(cache.G_F * cache.alpha_em, 2) / (64 * std::pow(M_PI, 3)) * cache.eta_BBS;
    return pref * std::pow(f_B * std::abs(lambda_q), 2) * std::pow(m_B, 3) * tau_B * beta_l * (
        std::pow(beta_l * std::abs(r_q * cache.CQ1_m), 2) 
        + std::pow(std::abs(r_q * cache.CQ2_m + 2 * x_l * cache.C10_m), 2)
    );
}

double BllDecay::BR_untag_Bs_ll(int gen) {
    double m_l = gen == 1 ? cache.m_e : gen == 2 ? cache.m_mu : cache.m_tau;
    double x_l = m_l / cache.m_Bs;
    double beta_l = std::sqrt(1. - 4. * std::pow(x_l, 2));
    double f = cache.r_s / (2. * x_l);
    complex_t S = beta_l * f * cache.CQ1_m / cache.C10_SM;
    complex_t P = (cache.C10_m + f * cache.CQ2_m) / cache.C10_SM;
    double abs_S = std::pow(std::abs(S), 2);
    double abs_P = std::pow(std::abs(P), 2);
    double A = (abs_P * std::cos(2 * std::arg(P)) - abs_S * std::cos(2 * std::arg(S))) / (abs_P + abs_S);

    double untag_factor = (1. + A * cache.ys) / (1. - std::pow(cache.ys, 2));
    return untag_factor * BR_avg_Bq_ll(3, gen);
}

std::vector<ObservableValue> BllDecay::compute_observable(Observables obs) {
    double value;
    switch (obs) {
    case Observables::BR_BD_MUMU:   
        value = BR_avg_Bq_ll(1, 2);
        break;
    case Observables::BR_BS_MUMU:   
        value = BR_avg_Bq_ll(3, 2);
        break;
    case Observables::BR_BS_MUMU_UNTAG:   
        value = BR_untag_Bs_ll(2);
        break;
    case Observables::BR_BD_EE:   
        value = BR_avg_Bq_ll(1, 1);
        break;
    case Observables::BR_BS_EE:   
        value = BR_avg_Bq_ll(3, 1);
        break;
    case Observables::BR_BS_EE_UNTAG:   
        value = BR_untag_Bs_ll(1);
        break;
    default:
        LOG_ERROR("IndexError", "Observable", ObservableMapper::str(obs), "doesn't belong to the decay", DecayMapper::str(this->id));
    }

    return {ObservableValue(ObservableMapper::to_id(obs), value)};
}

std::vector<ObservableValue> BllDecay::compute_observable(ObservableId obs) {
    return compute_observable(ObservableMapper::enum_of(obs).value());
}
