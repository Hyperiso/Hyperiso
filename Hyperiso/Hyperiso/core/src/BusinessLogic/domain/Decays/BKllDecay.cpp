#include "BKllDecay.h"

using Charge = BKllConfig::B_Charge;

void BKllDecay::load_params() {
    fill_wilson_cache();

    cache.ff_calculator = BPFFCalculator(
        cfg.charge == Charge::B_0 ? 511 : 521,
        cfg.charge == Charge::B_0 ? 311 : 321,
        cfg.ff_src
    );

    cache.qcdf_calculator = BPQCDfCalculator(
        cfg.charge == Charge::B_0 ? 511 : 521,
        cfg.charge == Charge::B_0 ? 311 : 321,
        w_config.hadronic_scale,
        cache.C,
        std::make_shared<BPFFCalculator>(cache.ff_calculator),
        cfg.ff_type
    );

    ObsParameterProxy p;
    cache.alpha_em = p(ParamId{ParameterType::SM, "EW", {1, 2}});
    cache.G_F = p(ParamId{ParameterType::SM, "SMINPUTS", 2});
    cache.m_l = p(ParamId{ParameterType::SM, "MASS", 11 + 2 * (int)cfg.gen});
    cache.m_s = p(ParamId{ParameterType::SM, "MASS", 3});
    cache.mu_b = w_config.hadronic_scale;
    cache.alpha_s_mu_b = ObsQCDProxy()(AlphasConfig(cache.mu_b, MassType::POLE, MassType::POLE));
    // cache.m_c_mu_b = ObsQCDProxy()(MassConfig(4, cache.mu_b, MassType::MSBAR, MassType::POLE));
    cache.m_c_mu_b = p(ParamId{ParameterType::SM, "MASS", 4}); // To match Superiso
    cache.m_b_mu_b = ObsQCDProxy()(MassConfig(5, cache.mu_b, MassType::MSBAR, MassType::POLE));
    double mu_f = sqrt(cache.mu_b * p(ParamId{ParameterType::DECAY, "B_K", 14}));
    cache.m_b_PS = p(ParamId{ParameterType::SM, "QCD", {5, 2}}) - 4 * ObsQCDProxy()(AlphasConfig(p(ParamId{ParameterType::SM, "QCD", {5, 2}}), MassType::POLE, MassType::POLE)) * mu_f / (3 * PI);
    cache.L_b = std::log(cache.mu_b / cache.m_b_PS);
    cache.m_B = p(ParamId{ParameterType::FLAVOR, "FMASS", cfg.charge == Charge::B_0 ? 511 : 521});
    cache.m_K = p(ParamId{ParameterType::FLAVOR, "FMASS", cfg.charge == Charge::B_0 ? 311 : 321});
    cache.Delta_M = -6. * cache.L_b - 4. * (1 - mu_f / cache.m_b_PS);
    cache.lambda_hat_u = std::conj(p(ParamId{ParameterType::SM, "VCKM", {0, 1}})) * p(ParamId{ParameterType::SM, "VCKM", {0, 2}}) 
                            / (std::conj(p(ParamId{ParameterType::SM, "VCKM", {2, 1}})) * p(ParamId{ParameterType::SM, "VCKM", {2, 2}}));
    cache.N_0 = std::pow(std::abs(std::conj(p(ParamId{ParameterType::SM, "VCKM", {2, 1}})) * p(ParamId{ParameterType::SM, "VCKM", {2, 2}})) * cache.G_F * cache.alpha_em, 2) / (512. * std::pow(PI, 5) * std::pow(cache.m_B, 3));
    cache.q2_min = 4 * std::pow(cache.m_l, 2);
    cache.q2_max = std::pow(cache.m_B - cache.m_K, 2);
    cache.q2_low = p(ParamId{ParameterType::DECAY, "B_K", {15, 1}});
    cache.q2_high = p(ParamId{ParameterType::DECAY, "B_K", {15, 2}});

    printf("alpha_em = %.4e\n", cache.alpha_em);
    printf("m_l = %.4e\n", cache.m_l);
    printf("m_s = %.4e\n", cache.m_s);
    printf("mu_b = %.4e\n", cache.mu_b);
    printf("alpha_s(mu_b) = %.4e\n", cache.alpha_s_mu_b);
    printf("m_c(mu_b) = %.4e\n", cache.m_c_mu_b);
    printf("m_b(mu_b) = %.4e\n", cache.m_b_mu_b);
    printf("m_b_PS = %.4e\n", cache.m_b_PS);
    printf("L_b = %.4e\n", cache.L_b);
    printf("m_B = %.4e\n", cache.m_B);
    printf("m_K = %.4e\n", cache.m_K);
    printf("Delta_M = %.4e\n", cache.Delta_M);
    printf("lambda_hat_u = %.4e + %.4e i\n", cache.lambda_hat_u.real(), cache.lambda_hat_u.imag());
    printf("N_0 = %.4e\n", cache.N_0);

    printf("f_0(s = 1.0 GeV²) = %.4e\n", cache.ff_calculator.get(BP_FF::F_0, 1.0));
    printf("f_+(s = 1.0 GeV²) = %.4e\n", cache.ff_calculator.get(BP_FF::F_PLUS, 1.0));
    printf("f_T(s = 1.0 GeV²) = %.4e\n", cache.ff_calculator.get(BP_FF::F_T, 1.0));

    // auto lam_T_P = [this] (double q2, bool bar) { return cache.qcdf_calculator.T_P(q2, bar); };
    // fill_cache(lam_T_P, cache.q2_min, cache.q2_high, cache.T_P_lookup, false); 

    // compute_binned_abc();
}

void BKllDecay::fill_wilson_cache() {
    auto b_wilsons = w_proxy->getAFR(WGroup::B, this->w_config.order);
    auto bp_wilsons = w_proxy->getAFR(WGroup::BPrime, this->w_config.order);
    auto bq_wilsons = w_proxy->getAFR(WGroup::BScalar, this->w_config.order);
    WCoef bp_cached[5] {WCoef::CP7, WCoef::CP9, WCoef::CP10, WCoef::CPQ1, WCoef::CPQ2};

    for (auto p : b_wilsons) cache.C.emplace(p); 
    for (auto p : bq_wilsons) cache.C.emplace(p);
    for (auto id : bp_cached) cache.C.emplace(std::pair{id, bp_wilsons.at(id)});
}

complex_t BKllDecay::T_P_cached(double q2) {
    return lerp(q2, cache.T_P_lookup, cache.q2_min, cache.q2_high);
}

double BKllDecay::beta_l(double q2) {
    return std::sqrt(1 - std::pow(2 * cache.m_l, 2) / q2);
}

double BKllDecay::lambda(double q2) {
    double mB2 = cache.m_B * cache.m_B;
    double mK2 = cache.m_K * cache.m_K;
    return mB2 * mB2 + mK2 * mK2 + q2 * q2 - 2. * (mB2 * mK2 + (mB2 + mK2) * q2);
}

double BKllDecay::N(double q2) {
    return cache.N_0 * std::sqrt(lambda(q2)) * beta_l(q2);
}

complex_t BKllDecay::F_V_low(double q2) {
    complex_t F, F_T;
    double m_b_local;
    complex_t had_err_factor {1.0};

    if (cfg.ff_type == B_FF_Type::SOFT) {
        had_err_factor = 1.0 + cache.A_had_err_low_0[0] + cache.A_had_err_low_1[0] * q2 / 6.0;
        F_T = T_P_cached(q2);
        F = (cache.C[WCoef::C9] + cache.C[WCoef::CP9]) * cache.ff_calculator.get(BP_FF::XI_P, q2);
        m_b_local = cache.m_b_PS;
    } else {
        complex_t T_p_err_factor = 1.0 + cache.A_had_err_low_0[0] + cache.A_had_err_low_1[0] * q2 / 6.0;
        F_T = (cache.C[WCoef::C7] + cache.C[WCoef::CP7]) * cache.ff_calculator.get(BP_FF::F_T, q2) + T_P_cached(q2) * T_p_err_factor;
        F = (cache.C[WCoef::C9] + cache.C[WCoef::CP9] + cache.qcdf_calculator.Y(q2)) * cache.ff_calculator.get(BP_FF::F_PLUS, q2);
        m_b_local = cache.m_b_PS + cache.alpha_s_mu_b * cache.Delta_M / (3 * PI);
    }

    return (F + 2. * m_b_local / (cache.m_B + cache.m_K) * F_T) * had_err_factor;
}

complex_t BKllDecay::F_A_low(double q2) {
    double ff;

    if (cfg.ff_type == B_FF_Type::SOFT) {
        ff = cache.ff_calculator.get(BP_FF::XI_P, q2);
    } else {
        ff = cache.ff_calculator.get(BP_FF::F_PLUS, q2);
    }

    return (cache.C[WCoef::C10] + cache.C[WCoef::CP10]) * ff;
}

complex_t BKllDecay::F_P_low(double q2) {
    double f_p, f0_fp;

    if (cfg.ff_type == B_FF_Type::SOFT) {
        f_p = cache.ff_calculator.get(BP_FF::XI_P, q2);
        f0_fp = 2. * (std::pow(cache.m_B, 2) + std::pow(cache.m_K, 2) - q2) / (2 * std::pow(cache.m_B, 2)) * cache.qcdf_calculator.Delta_P_0(q2);
    } else {
        f_p = cache.ff_calculator.get(BP_FF::F_PLUS, q2);
        f0_fp = cache.ff_calculator.get(BP_FF::F_0, q2) / f_p;
    }

    return f_p * ((std::pow(cache.m_B, 2) - std::pow(cache.m_K, 2)) / (2 * (cache.m_b_mu_b - cache.m_s)) * f0_fp * (cache.C[WCoef::CQ2] + cache.C[WCoef::CPQ2])
                    - cache.m_l * (cache.C[WCoef::C10] + cache.C[WCoef::CP10]) * (1. - (std::pow(cache.m_B, 2) - std::pow(cache.m_K, 2)) / q2 * (f0_fp - 1.)));
}

complex_t BKllDecay::F_S_low(double q2) {
    double ff;

    if (cfg.ff_type == B_FF_Type::SOFT) {
        double f0_fp = 2. * (std::pow(cache.m_B, 2) + std::pow(cache.m_K, 2) - q2) / (2 * std::pow(cache.m_B, 2)) * cache.qcdf_calculator.Delta_P_0(q2);
        ff = cache.ff_calculator.get(BP_FF::XI_P, q2) * f0_fp;
    } else {
        ff = cache.ff_calculator.get(BP_FF::F_0, q2);
    }

    return (std::pow(cache.m_B, 2) - std::pow(cache.m_K, 2)) / (2 * (cache.m_b_mu_b - cache.m_s)) * (cache.C[WCoef::CQ1] + cache.C[WCoef::CPQ1]) * ff;
}

complex_t BKllDecay::C7_eff(double q2) {
    double s_hat = q2 / std::pow(cache.m_b_PS, 2);
    complex_t A = BV::A_Seidel(s_hat, cache.L_b);
    return cache.C[WCoef::C7] + cache.alpha_s_mu_b / (4. * PI) * ((cache.C[WCoef::C1]-6.*cache.C[WCoef::C2])*A-cache.C[WCoef::C8]*BV::f_87(q2 / (cache.m_B * cache.m_B), cache.L_b));
}

complex_t BKllDecay::C9_eff(double q2) {
    complex_t C_h0 = 4./3.*cache.C[WCoef::C1]+cache.C[WCoef::C2]+11./2.*cache.C[WCoef::C3]-2./3.*cache.C[WCoef::C4]+52.*cache.C[WCoef::C5]-32./3.*cache.C[WCoef::C6];
    complex_t C_hb = -0.5 * (7.*cache.C[WCoef::C3]+4./3.*cache.C[WCoef::C4]+76.*cache.C[WCoef::C5]+64./3.*cache.C[WCoef::C6]);
    complex_t C_0  = 4./3.*(cache.C[WCoef::C3]+16./3.*cache.C[WCoef::C5]+16./9.*cache.C[WCoef::C6]);
    complex_t C_mc = 8. * ((4./9.*cache.C[WCoef::C1]+1./3.*cache.C[WCoef::C2])*(1.+cache.lambda_hat_u)+2.*cache.C[WCoef::C3]+20.*cache.C[WCoef::C5]);

    double s_hat = q2 / std::pow(cache.m_b_PS, 2);
    complex_t A = BV::A_Seidel(s_hat, cache.L_b);
    complex_t B = BV::B_Seidel(s_hat, cache.L_b);
    complex_t C = BV::C_Seidel(q2, cache.mu_b);

    return cache.C[WCoef::C9]
         + BV::h(q2, 0., cache.mu_b) * C_h0
         + BV::h(q2, cache.m_b_PS, cache.mu_b) * C_hb
         + C_0
         + cache.alpha_s_mu_b / (4. * PI) * (cache.C[WCoef::C1]*(B + 4. * C) - 3. * cache.C[WCoef::C2] * (2. * B - C) - cache.C[WCoef::C8] * BV::f_89(q2 / (cache.m_B * cache.m_B)))
         + std::pow(cache.m_c_mu_b, 2) / q2 * C_mc;
}

complex_t BKllDecay::F_V_high(double q2) {
    return (C9_eff(q2) + cache.C[WCoef::CP9]) * cache.ff_calculator.get(BP_FF::F_PLUS, q2) 
            + 2. * cache.m_b_mu_b / (cache.m_B + cache.m_K) * (C7_eff(q2) + cache.C[WCoef::CP7]) * cache.ff_calculator.get(BP_FF::F_T, q2);
}

complex_t BKllDecay::F_A_high(double q2) {
    return (cache.C[WCoef::C10] + cache.C[WCoef::CP10]) * cache.ff_calculator.get(BP_FF::F_PLUS, q2);
}

complex_t BKllDecay::F_P_high(double q2) {
    double f_p = cache.ff_calculator.get(BP_FF::F_PLUS, q2);
    double f0_fp = cache.ff_calculator.get(BP_FF::F_0, q2) / f_p;
    return f_p * ((std::pow(cache.m_B, 2) - std::pow(cache.m_K, 2)) / (2 * (cache.m_b_mu_b - cache.m_s)) * f0_fp * (cache.C[WCoef::CQ2] + cache.C[WCoef::CPQ2])
                    - cache.m_l * (cache.C[WCoef::C10] + cache.C[WCoef::CP10]) * (1. - (std::pow(cache.m_B, 2) - std::pow(cache.m_K, 2)) / q2 * (f0_fp - 1.)));
}

complex_t BKllDecay::F_S_high(double q2) {
    double ff = cache.ff_calculator.get(BP_FF::F_0, q2);
    return (std::pow(cache.m_B, 2) - std::pow(cache.m_K, 2)) / (2 * (cache.m_b_mu_b - cache.m_s)) * (cache.C[WCoef::CQ1] + cache.C[WCoef::CPQ1]) * ff;
}

complex_t BKllDecay::interpolate(double q2, complex_t val_low, complex_t val_high) {
    if (q2 < cache.q2_low)
        return val_low;

    if (q2 > cache.q2_high)
        return val_high;

    double t = (cache.q2_high - q2) / (cache.q2_high - cache.q2_low);
    return t * val_low + (1 - t) * val_high;
}

complex_t BKllDecay::F_V(double q2) {
    return interpolate(q2, F_V_low(q2), F_V_high(q2));
}

complex_t BKllDecay::F_A(double q2) {
    return interpolate(q2, F_A_low(q2), F_A_high(q2));
}

complex_t BKllDecay::F_P(double q2) {
    return interpolate(q2, F_P_low(q2), F_P_high(q2));
}

complex_t BKllDecay::F_S(double q2) {
    return interpolate(q2, F_S_low(q2), F_S_high(q2));
}

double BKllDecay::a(double q2) {
    return N(q2) * (
        q2 * (std::pow(beta_l(q2) * std::abs(F_S(q2)), 2) + std::pow(std::abs(F_P(q2)), 2))
      + 0.25 * lambda(q2) * (std::pow(std::abs(F_A(q2)), 2) + std::pow(std::abs(F_V(q2)), 2))
      + 2 * cache.m_l * (std::pow(cache.m_B, 2) - std::pow(cache.m_K, 2) + q2) * std::real(F_P(q2) * std::conj(F_A(q2)))
      + std::pow(2 * cache.m_l * cache.m_B, 2) * std::pow(std::abs(F_A(q2)), 2)
    );
}

double BKllDecay::b(double q2) {
    return 2 * N(q2) * cache.m_l * beta_l(q2) * std::sqrt(lambda(q2)) * std::real(F_S(q2) * std::conj(F_V(q2)));
}

double BKllDecay::c(double q2) {
    return -0.25 * N(q2) * lambda(q2) * std::pow(beta_l(q2), 2) * (std::pow(std::abs(F_A(q2)), 2) + std::pow(std::abs(F_V(q2)), 2));
}

void BKllDecay::compute_binned_abc() {
    for (auto [q2_l, q2_u] : cfg.bins) {
        cache.abc_binned[0].emplace_back(integrate([&] (double q2) { return a(q2); }, q2_l, q2_u, 1e-3));    
        cache.abc_binned[1].emplace_back(integrate([&] (double q2) { return b(q2); }, q2_l, q2_u, 1e-3));
        cache.abc_binned[2].emplace_back(integrate([&] (double q2) { return c(q2); }, q2_l, q2_u, 1e-3));
    }
}

std::vector<ObservableValue> BKllDecay::dG_dq2() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::DGAMMA_DQ2_B__K_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double res = 2 * (cache.abc_binned[0][i] + cache.abc_binned[2][i] / 3.); 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKllDecay::A_FB() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::A_FB_B__K_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double res = cache.abc_binned[1][i] / (2 * (cache.abc_binned[0][i] + cache.abc_binned[2][i] / 3.)); 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKllDecay::F_H() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::F_H_B__K_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double res = (cache.abc_binned[0][i] + cache.abc_binned[2][i]) / (cache.abc_binned[0][i] + cache.abc_binned[2][i] / 3.); 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKllDecay::compute_observable(Observables obs) {
    switch (obs) {
    case Observables::TEST:
        return {};
    case Observables::DGAMMA_DQ2_B__K_L_L:   
        return dG_dq2();
    case Observables::A_FB_B__K_L_L:   
        return A_FB();
    case Observables::F_H_B__K_L_L:   
        return F_H();
    default:
        LOG_ERROR("IndexError", "Observable", ObservableMapper::str(obs), "doesn't belong to the decay", DecayMapper::str(this->id));
    }
}

std::vector<ObservableValue> BKllDecay::compute_observable(ObservableId obs) {
    return compute_observable(ObservableMapper::enum_of(obs).value());
}
