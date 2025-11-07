#include "BllDecay.h"


void BllDecay::load_params() {
    ObsParameterProxy p;
    cache.G_F = p(ParamId{ParameterType::SM, "SMINPUTS", 2});
    cache.alpha_em = p(ParamId{ParameterType::SM, "EW", {1, 2}});
    cache.m_mu = p(ParamId{ParameterType::SM, "MASS", 13});
    cache.m_Bd = p(ParamId{ParameterType::FLAVOR, "FMASS", 511});
    cache.m_Bs = p(ParamId{ParameterType::FLAVOR, "FMASS", 531});
    cache.f_Bd = p(ParamId{ParameterType::FLAVOR, "FCONST", {511, 1}});
    cache.f_Bs = p(ParamId{ParameterType::FLAVOR, "FCONST", {531, 1}});
    cache.tau_Bd = p(ParamId{ParameterType::FLAVOR, "FLIFE", 511}) / HBAR;
    cache.tau_Bs = p(ParamId{ParameterType::FLAVOR, "FLIFE", 531}) / HBAR;
    cache.lambda_d = p(ParamId{ParameterType::SM, "VCKM", {2, 2}}) * std::conj(p(ParamId{ParameterType::SM, "VCKM", {2, 0}}));
    cache.lambda_s = p(ParamId{ParameterType::SM, "VCKM", {2, 2}}) * std::conj(p(ParamId{ParameterType::SM, "VCKM", {2, 1}}));
    cache.ys = p(ParamId{ParameterType::DECAY, "B_ll", 1});
    cache.x_d = cache.m_mu / cache.m_Bd;
    cache.x_s = cache.m_mu / cache.m_Bs;
    cache.r_d = cache.m_Bd / (p(ParamId{ParameterType::SM, "QCD", {5, 2}}) + p(ParamId{ParameterType::SM, "MASS", 1}));
    cache.r_s = cache.m_Bs / (p(ParamId{ParameterType::SM, "QCD", {5, 2}}) + p(ParamId{ParameterType::SM, "MASS", 3}));
    cache.beta_d = std::sqrt(1. - 4. * std::pow(cache.x_d, 2));
    cache.beta_s = std::sqrt(1. - 4. * std::pow(cache.x_s, 2));
    cache.C10_SM = w_proxy->getFR(WGroup::B, WCoef::C10, w_config.order, ContributionType::SM);
    cache.C10 = w_proxy->getFR(WGroup::B, WCoef::C10, w_config.order);
    cache.CQ1 = w_proxy->getFR(WGroup::BScalar, WCoef::CQ1, w_config.order);
    cache.CQ2 = w_proxy->getFR(WGroup::BScalar, WCoef::CQ2, w_config.order);
    cache.C10_m = cache.C10 - w_proxy->getFR(WGroup::BPrime, WCoef::CP10, w_config.order);
    cache.CQ1_m = cache.CQ1 - w_proxy->getFR(WGroup::BPrime, WCoef::CPQ1, w_config.order);
    cache.CQ2_m = cache.CQ2 - w_proxy->getFR(WGroup::BPrime, WCoef::CPQ2, w_config.order);

    printf("1 / alpha_em = %.4e\n", 1 / cache.alpha_em);
    printf("m_mu = %.4e\n", cache.m_mu);
    printf("m_Bd = %.4e\n", cache.m_Bd);
    printf("m_Bs = %.4e\n", cache.m_Bs);
    printf("f_Bd = %.4e\n", cache.f_Bd);
    printf("f_Bs = %.4e\n", cache.f_Bs);
    printf("tau_Bd = %.4e\n", cache.tau_Bd);
    printf("tau_Bs = %.4e\n", cache.tau_Bs);
    printf("lambda_d = %.4e\n", std::abs(cache.lambda_d));
    printf("lambda_s = %.4e\n", std::abs(cache.lambda_s));
    printf("y_s = %.4e\n", cache.ys);
    printf("r_d = %.4e\n", cache.r_d);
    printf("r_s = %.4e\n", cache.r_s);
    printf("C10 = %.4e\n", cache.C10.real());
    printf("CQ1 = %.4e\n", cache.CQ1.real());
    printf("CQ2 = %.4e\n", cache.CQ2.real());
}

double BllDecay::BR_avg_Bq_mumu(int q) {
    if (q != 1 && q != 3) LOG_ERROR("ValueError", "In Bq > mu mu, q can only be d (1) or s (3), found", q);

    double pref = std::pow(cache.G_F * cache.alpha_em, 2) / (64 * std::pow(M_PI, 3));
    if (q == 1) {
        return pref * std::pow(cache.f_Bd * std::abs(cache.lambda_d), 2) * std::pow(cache.m_Bd, 3) * cache.tau_Bd * cache.beta_d * (
            std::pow(cache.beta_d * std::abs(cache.r_d * cache.CQ1_m), 2) 
          + std::pow(std::abs(cache.r_d * cache.CQ2_m + 2 * cache.x_d * cache.C10_m), 2)
        ) * 0.995;
    } else {
        return pref * std::pow(cache.f_Bs * std::abs(cache.lambda_s), 2) * std::pow(cache.m_Bs, 3) * cache.tau_Bs * cache.beta_s * (
            std::pow(cache.beta_s * std::abs(cache.r_s * cache.CQ1_m), 2) 
          + std::pow(std::abs(cache.r_s * cache.CQ2_m + 2 * cache.x_s * cache.C10_m), 2)
        ) * 0.995;
    }
}

double BllDecay::BR_untag_Bs_mumu() {
    double f = cache.r_s / (2. * cache.x_s);
    complex_t S = cache.beta_s * f * cache.CQ1_m / cache.C10_SM;
    complex_t P = (cache.C10_m + f * cache.CQ2_m) / cache.C10_SM;
    double abs_S = std::pow(std::abs(S), 2);
    double abs_P = std::pow(std::abs(P), 2);
    double A = (abs_P * std::cos(2 * std::arg(P)) - abs_S * std::cos(2 * std::arg(S))) / (abs_P + abs_S);

    double untag_factor = (1. + A * cache.ys) / (1. - std::pow(cache.ys, 2));
    return untag_factor * BR_avg_Bq_mumu(3);
}

std::vector<ObservableValue> BllDecay::compute_observable(Observables obs) {
    double value;
    switch (obs) {
    case Observables::BR_BD_MUMU:   
        value = BR_avg_Bq_mumu(1);
        break;
    case Observables::BR_BS_MUMU:   
        value = BR_avg_Bq_mumu(3);
        break;
    case Observables::BR_BS_MUMU_UNTAG:   
        value = BR_untag_Bs_mumu();
        break;
    default:
        LOG_ERROR("IndexError", "Observable", ObservableMapper::str(obs), "doesn't belong to the decay", DecayMapper::str(this->id));
    }

    return {ObservableValue(ObservableMapper::to_id(obs), value)};
}

std::vector<ObservableValue> BllDecay::compute_observable(ObservableId obs) {
    return compute_observable(ObservableMapper::enum_of(obs).value());
}
