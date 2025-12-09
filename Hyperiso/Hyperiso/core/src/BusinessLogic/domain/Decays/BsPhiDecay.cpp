#include "BsPhiDecay.h"

void BsPhiDecay::load_params() {
    LOG_INFO("Loading parameters for Bs > phi ll decay");
    fill_wilson_cache();

    cache.ff_calculator = BVFFCalculator(531, 333, cfg.ff_src);

    cache.qcdf_calculator = BVQCDfCalculator(
        531, 333,
        w_config.hadronic_scale,
        cache.C,
        std::make_shared<BVFFCalculator>(cache.ff_calculator),
        cfg.ff_type
    );

    ObsParameterProxy p;
    cache.alpha_em = p(ParamId{ParameterType::SM, "EW", {1, 2}});
    cache.G_F = p(ParamId{ParameterType::SM, "SMINPUTS", 2});
    cache.m_l = p(ParamId{ParameterType::SM, "MASS", 11 + 2 * (int)cfg.gen});
    cache.m_s = p(ParamId{ParameterType::SM, "MASS", 3});
    cache.mu_b = w_config.hadronic_scale;
    cache.alpha_s_mu_b = ObsQCDProxy()(AlphasConfig(cache.mu_b, MassType::POLE, MassType::POLE));
    cache.m_c_mu_b = ObsQCDProxy()(MassConfig(4, cache.mu_b, MassType::MSBAR, MassType::POLE));
    cache.m_b_mu_b = ObsQCDProxy()(MassConfig(5, cache.mu_b, MassType::MSBAR, MassType::POLE));
    cache.m_b_PS = p(ParamId{ParameterType::SM, "QCD", {5, 2}}) - 4 * ObsQCDProxy()(AlphasConfig(p(ParamId{ParameterType::SM, "QCD", {5, 2}}), MassType::POLE, MassType::POLE)) * sqrt(cache.mu_b * p(ParamId{ParameterType::DECAY, "B_phi", 14})) / (3 * PI);
    cache.L_b = std::log(cache.mu_b / cache.m_b_PS);
    cache.m_Bs = p(ParamId{ParameterType::FLAVOR, "FMASS", 531});
    cache.m_phi = p(ParamId{ParameterType::FLAVOR, "FMASS", 333});
    cache.lambda_hat_u = std::conj(p(ParamId{ParameterType::SM, "VCKM", {0, 1}})) * p(ParamId{ParameterType::SM, "VCKM", {0, 2}}) 
                            / (std::conj(p(ParamId{ParameterType::SM, "VCKM", {2, 1}})) * p(ParamId{ParameterType::SM, "VCKM", {2, 2}}));
    cache.kappa = 1 - 2. * cache.alpha_s_mu_b / (3. * PI) * std::log(cache.mu_b / cache.m_b_mu_b);
    cache.Delta_M = -6. * cache.L_b - 4. * (1 - sqrt(cache.mu_b * p(ParamId{ParameterType::DECAY, "B_phi", 14})) / cache.m_b_PS);
    cache.N_0 = std::conj(p(ParamId{ParameterType::SM, "VCKM", {2, 1}})) * p(ParamId{ParameterType::SM, "VCKM", {2, 2}}) * cache.G_F * cache.alpha_em / (std::sqrt(3072. * std::pow(PI, 5) * std::pow(cache.m_Bs, 3)));
    cache.q2_min = 4 * std::pow(cache.m_l, 2);
    cache.q2_max = std::pow(cache.m_Bs - cache.m_phi, 2);
    cache.q2_low = p(ParamId{ParameterType::DECAY, "B_phi", {15, 1}});
    cache.q2_high = p(ParamId{ParameterType::DECAY, "B_phi", {15, 2}});
    cache.ys = p(ParamId{ParameterType::DECAY, "B_ll", 1});
    // TODO : hardcoded in SI
    double beta_s = std::arg(-std::conj(p(ParamId{ParameterType::SM, "VCKM", {2, 2}})) * p(ParamId{ParameterType::SM, "VCKM", {2, 1}}) 
                            / (std::conj(p(ParamId{ParameterType::SM, "VCKM", {1, 2}})) * p(ParamId{ParameterType::SM, "VCKM", {1, 1}})));
    cache.phi_s = 2 * beta_s;
    cache.up = std::exp(I * cache.phi_s);
    cache.um = std::exp(-I * cache.phi_s);

    auto lam_T_perp_p = [this] (double q2, bool bar) { return cache.qcdf_calculator.T_perp_p(q2, bar); };
    fill_cache(lam_T_perp_p, cache.q2_min, cache.q2_high, cache.T_perp_p_lookup, false); 
    fill_cache(lam_T_perp_p, cache.q2_min, cache.q2_high, cache.T_perp_p_bar_lookup, true); 

    auto lam_T_perp_m = [this] (double q2, bool bar) { return cache.qcdf_calculator.T_perp_m(q2, bar); };
    fill_cache(lam_T_perp_m, cache.q2_min, cache.q2_high, cache.T_perp_m_lookup, false); 
    fill_cache(lam_T_perp_m, cache.q2_min, cache.q2_high, cache.T_perp_m_bar_lookup, true);

    auto lam_T_par_m = [this] (double q2, bool bar) { return cache.qcdf_calculator.T_par_m(q2, bar); };
    fill_cache(lam_T_par_m, cache.q2_min, cache.q2_high, cache.T_par_m_lookup, false); 
    fill_cache(lam_T_par_m, cache.q2_min, cache.q2_high, cache.T_par_m_bar_lookup, true);

    double q2 = 1.0;
    double u = 0.5;
    complex_t Tperpp = cache.qcdf_calculator.T_perp_p(q2, false);
    complex_t Tperpm = cache.qcdf_calculator.T_perp_m(q2, false);
    complex_t Tparm = cache.qcdf_calculator.T_par_m(q2, false);

    // printf("T_perp_p(s = %.3f) = %.4e + %.4e i\n", q2, std::real(Tperpp), std::imag(Tperpp));
    // printf("T_perp_m(s = %.3f) = %.4e + %.4e i\n", q2, std::real(Tperpm), std::imag(Tperpm));
    // printf("T_par_m(s = %.3f) = %.4e + %.4e i\n", q2, std::real(Tparm), std::imag(Tparm));
    // exit(0);

    // complex_t ALperp = A_perp(q2, -1, false);
    // complex_t ARperp = A_perp(q2, 1, false);
    // complex_t ALpar = A_par(q2, -1, false);
    // complex_t ARpar = A_par(q2, 1, false);
    // complex_t AL0 = A_0(q2, -1, false);
    // complex_t AR0 = A_0(q2, 1, false);
    // complex_t At = A_t(q2, false);
    // complex_t AS = A_S(q2, false);

    // complex_t ALperp = A_perp(q2, -1, true);
    // complex_t ARperp = A_perp(q2, 1, true);
    // complex_t ALpar = A_par(q2, -1, true);
    // complex_t ARpar = A_par(q2, 1, true);
    // complex_t AL0 = A_0(q2, -1, true);
    // complex_t AR0 = A_0(q2, 1, true);
    // complex_t At = A_t(q2, true);
    // complex_t AS = A_S(q2, true);

    // printf("A_L_perp(s = %.3f) = %.4e + %.4e i\n", q2, std::real(ALperp), std::imag(ALperp));
	// printf("A_R_perp(s = %.3f) = %.4e + %.4e i\n", q2, std::real(ARperp), std::imag(ARperp));
	// printf("A_L_par(s = %.3f) = %.4e + %.4e i\n", q2, std::real(ALpar), std::imag(ALpar));
	// printf("A_R_par(s = %.3f) = %.4e + %.4e i\n", q2, std::real(ARpar), std::imag(ARpar));
	// printf("A_L_0(s = %.3f) = %.4e + %.4e i\n", q2, std::real(AL0), std::imag(AL0));
	// printf("A_R_0(s = %.3f) = %.4e + %.4e i\n", q2, std::real(AR0), std::imag(AR0));
	// printf("A_t(s = %.3f) = %.4e + %.4e i\n", q2, std::real(At), std::imag(At));
	// printf("A_S(s = %.3f) = %.4e + %.4e i\n", q2, std::real(AS), std::imag(AS));
    // exit(0);
    
    // printf("J_1c(s = %.3f) = %.4e\n", q2, J1c(q2, false));
    // printf("J_1s(s = %.3f) = %.4e\n", q2, J1s(q2, false));
    // printf("J_2c(s = %.3f) = %.4e\n", q2, J2c(q2, false));
    // printf("J_2s(s = %.3f) = %.4e\n", q2, J2s(q2, false));
    // printf("J_3(s = %.3f) = %.4e\n", q2, J3(q2, false));
    // printf("J_4(s = %.3f) = %.4e\n", q2, J4(q2, false));
    // printf("J_5(s = %.3f) = %.4e\n", q2, J5(q2, false));
    // printf("J_6c(s = %.3f) = %.4e\n", q2, J6c(q2, false));
    // printf("J_6s(s = %.3f) = %.4e\n", q2, J6s(q2, false));
    // printf("J_7(s = %.3f) = %.4e\n", q2, J7(q2, false));
    // printf("J_8(s = %.3f) = %.4e\n", q2, J8(q2, false));
    // printf("J_9(s = %.3f) = %.4e\n", q2, J9(q2, false));
    // exit(0);

    compute_binned_J_i();
}

void BsPhiDecay::fill_wilson_cache() {
    auto b_wilsons = w_proxy->getAFR(WGroup::B, this->w_config.order);
    auto bp_wilsons = w_proxy->getAFR(WGroup::BPrime, this->w_config.order);
    auto bq_wilsons = w_proxy->getAFR(WGroup::BScalar, this->w_config.order);
    WCoef bp_cached[5] {WCoef::CP7, WCoef::CP9, WCoef::CP10, WCoef::CPQ1, WCoef::CPQ2};

    for (auto p : b_wilsons) cache.C.emplace(p); 
    for (auto p : bq_wilsons) cache.C.emplace(p);
    for (auto id : bp_cached) cache.C.emplace(std::pair{id, bp_wilsons.at(id)});
}

complex_t BsPhiDecay::T_perp_p_cached(double q2, bool bar) {
    return lerp(q2, bar ? cache.T_perp_p_bar_lookup : cache.T_perp_p_lookup, cache.q2_min, cache.q2_high);
}

complex_t BsPhiDecay::T_perp_m_cached(double q2, bool bar) {
    return lerp(q2, bar ? cache.T_perp_m_bar_lookup : cache.T_perp_m_lookup, cache.q2_min, cache.q2_high);
}

complex_t BsPhiDecay::T_par_m_cached(double q2, bool bar) {
    return lerp(q2, bar ? cache.T_par_m_bar_lookup : cache.T_par_m_lookup, cache.q2_min, cache.q2_high);
}

double BsPhiDecay::beta_l(double q2) {
    return std::sqrt(1 - 4. * cache.m_l * cache.m_l / q2);
}

double BsPhiDecay::lambda(double q2) {
    double mB2 = cache.m_Bs * cache.m_Bs;
    double mphi2 = cache.m_phi * cache.m_phi;
    return mB2 * mB2 + mphi2 * mphi2 + q2 * q2 - 2. * (mB2 * mphi2 + (mB2 + mphi2) * q2);
}

complex_t BsPhiDecay::N(double q2, bool bar) {
    complex_t N0 = bar ? std::conj(cache.N_0) : cache.N_0; 
    return N0 * std::sqrt(q2 * beta_l(q2) * std::sqrt(lambda(q2)));
}

complex_t BsPhiDecay::delta_A_perp(double q2, double sign, bool bar) {
    size_t id = size_t (0.5 * (1 + sign));
    double guesstimate_err = 1.0 + cache.A_had_err_low_0[id] + cache.A_had_err_low_1[id] * q2 / 6.0;
    return 2 * RT2 * cache.m_b_PS * N(q2, bar) * std::sqrt(lambda(q2)) / q2 * T_perp_p_cached(q2, bar) * guesstimate_err;
}

complex_t BsPhiDecay::delta_A_par(double q2, double sign, bool bar) {
    size_t id = 2 + size_t (0.5 * (1 + sign));
    double guesstimate_err = 1.0 + cache.A_had_err_low_0[id] + cache.A_had_err_low_1[id] * q2 / 6.0;

    return -4 * RT2 * cache.m_b_PS * N(q2, bar) * (cache.m_Bs * cache.m_Bs - cache.m_phi * cache.m_phi) * cache.ff_calculator.E(q2) / (q2 * cache.m_Bs) * T_perp_m_cached(q2, bar) * guesstimate_err;
}

complex_t BsPhiDecay::A_perp_low(double q2, double sign, bool bar) {
    complex_t F, F_T;
    complex_t delta_A {0.0};
    complex_t had_err_factor {1.0};
    double m_b_local;

    if (cfg.ff_type == B_FF_Type::SOFT) {
        complex_t w = cache.C[WCoef::C9] + cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] + cache.C[WCoef::CP10]);
        if (bar) w = std::conj(w);
        F = w * cache.ff_calculator.get(BV_FF::V, q2) / (cache.m_Bs + cache.m_phi);
        F_T = T_perp_p_cached(q2, bar);
        size_t id = size_t (0.5 * (1 + sign));
        had_err_factor = 1.0 + cache.A_had_err_low_0[id] + cache.A_had_err_low_1[id] * q2 / 6.0;
        m_b_local = cache.m_b_PS;
    } else {
        F_T = (cache.C[WCoef::C7] + cache.C[WCoef::CP7]) * cache.ff_calculator.get(BV_FF::T1, q2);
        complex_t w = cache.C[WCoef::C9] + cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] + cache.C[WCoef::CP10]);
        if (bar) {
            w = std::conj(w);
            F_T = std::conj(F_T);
        } 
        w += cache.qcdf_calculator.Y(q2);
        F = w * cache.ff_calculator.get(BV_FF::V, q2) / (cache.m_Bs + cache.m_phi);
        delta_A = delta_A_perp(q2, sign, bar);
        m_b_local = cache.m_b_PS + cache.alpha_s_mu_b * cache.Delta_M / (3 * PI);
    }

    return (N(q2, bar) * std::sqrt(2 * lambda(q2)) * (F + 2. * m_b_local * F_T / q2) + delta_A) * had_err_factor;
}

complex_t BsPhiDecay::A_par_low(double q2, double sign, bool bar) {
    complex_t F, F_T;
    complex_t delta_A {0.0};
    complex_t had_err_factor {1.0};
    double m_b_local;

    if (cfg.ff_type == B_FF_Type::SOFT) {
        complex_t w = cache.C[WCoef::C9] - cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] - cache.C[WCoef::CP10]);
        if (bar) w = std::conj(w);
        F = w * cache.ff_calculator.get(BV_FF::A1, q2) / (cache.m_Bs - cache.m_phi);
        F_T = 2. * cache.ff_calculator.E(q2) * T_perp_m_cached(q2, bar) / cache.m_Bs;
        size_t id = 2 + size_t (0.5 * (1 + sign));
        had_err_factor = 1.0 + cache.A_had_err_low_0[id] + cache.A_had_err_low_1[id] * q2 / 6.0;
        m_b_local = cache.m_b_PS;
    } else {
        F_T = (cache.C[WCoef::C7] - cache.C[WCoef::CP7]) * cache.ff_calculator.get(BV_FF::T2, q2);
        complex_t w = cache.C[WCoef::C9] - cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] - cache.C[WCoef::CP10]);
        if (bar) {
            w = std::conj(w);
            F_T = std::conj(F_T);
        } 
        w += cache.qcdf_calculator.Y(q2);
        F = w * cache.ff_calculator.get(BV_FF::A1, q2) / (cache.m_Bs - cache.m_phi);
        delta_A = delta_A_par(q2, sign, bar);
        m_b_local = cache.m_b_PS + cache.alpha_s_mu_b * cache.Delta_M / (3 * PI);
    }

    return (-N(q2, bar) * std::sqrt(2.) * (cache.m_Bs * cache.m_Bs - cache.m_phi * cache.m_phi) * (F + 2. * m_b_local * F_T / q2) + delta_A) * had_err_factor;
}

complex_t BsPhiDecay::A_0_low(double q2, double sign, bool bar) {
    double mB2 = cache.m_Bs * cache.m_Bs;
    double mK2 = cache.m_phi * cache.m_phi;
    complex_t F, F_T;
    complex_t delta_A {0.0};
    complex_t had_err_factor {1.0};
    double m_b_local;

    if (cfg.ff_type == B_FF_Type::SOFT) {
        complex_t w = cache.C[WCoef::C9] - cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] - cache.C[WCoef::CP10]);
        if (bar) w = std::conj(w);
        F = w * 16. * cache.m_Bs * mK2 * cache.ff_calculator.get(BV_FF::A12, q2);
        F_T = 2. * cache.ff_calculator.E(q2) * (mB2 + 3. * mK2 - q2) * T_perp_m_cached(q2, bar) - lambda(q2) * (T_perp_m_cached(q2, bar) + T_par_m_cached(q2, bar)) / (mB2 - mK2);
        size_t id = 4 + size_t (0.5 * (1 + sign));
        had_err_factor = 1.0 + cache.A_had_err_low_0[id] + cache.A_had_err_low_1[id] * q2 / 6.0;
        m_b_local = cache.m_b_PS;
    } else {
        F_T = (cache.C[WCoef::C7] - cache.C[WCoef::CP7]) * 8. * cache.m_Bs * mK2 / (cache.m_Bs + cache.m_phi) * cache.ff_calculator.get(BV_FF::T23, q2);
        complex_t w = cache.C[WCoef::C9] - cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] - cache.C[WCoef::CP10]);
        if (bar) {
            w = std::conj(w);
            F_T = std::conj(F_T);
        } 
        w += cache.qcdf_calculator.Y(q2);
        F = w * 16. * cache.m_Bs * mK2 * cache.ff_calculator.get(BV_FF::A12, q2);
        delta_A = delta_A_par(q2, sign, bar);
        m_b_local = cache.m_b_PS + cache.alpha_s_mu_b * cache.Delta_M / (3 * PI);
    }

    return (-N(q2, bar) / (2. * cache.m_phi * std::sqrt(q2)) * (F + 2. * m_b_local * F_T) + delta_A) * had_err_factor;
}

complex_t BsPhiDecay::A_t_low(double q2, bool bar) {
    complex_t C10 = cache.C[WCoef::C10] - cache.C[WCoef::CP10];
    complex_t CQ2 = cache.C[WCoef::CQ2] - cache.C[WCoef::CPQ2];
    if (bar) {
        C10 = std::conj(C10);
        CQ2 = std::conj(CQ2);
    }

    complex_t F;
    if (cfg.ff_type == B_FF_Type::SOFT) {
        F = cache.ff_calculator.E(q2) * cache.ff_calculator.get(BV_FF::XI_PAR, q2) / (cache.m_phi * cache.qcdf_calculator.Delta_par(q2));
    } else {
        F = cache.ff_calculator.get(BV_FF::A0, q2);
    }

    return N(q2, bar) * std::sqrt(lambda(q2) / q2) * ((2. * C10 + q2 / (cache.m_l * (cache.m_b_mu_b + cache.m_s)) * CQ2)) * F;
}

complex_t BsPhiDecay::A_S_low(double q2, bool bar) {
    complex_t CQ1 = cache.C[WCoef::CQ1] - cache.C[WCoef::CPQ1];
    if (bar) CQ1 = std::conj(CQ1);

    complex_t F;
    if (cfg.ff_type == B_FF_Type::SOFT) {
        F = cache.ff_calculator.E(q2) * cache.ff_calculator.get(BV_FF::XI_PAR, q2) / (cache.m_phi * cache.qcdf_calculator.Delta_par(q2));
    } else {
        F = cache.ff_calculator.get(BV_FF::A0, q2);
    }

    return -2. * N(q2, bar) * std::sqrt(lambda(q2)) * CQ1 / (cache.m_b_mu_b + cache.m_s) * F;
}

complex_t BsPhiDecay::delta_A_0(double q2, double sign, bool bar) {
    size_t id = 4 + size_t (0.5 * (1 + sign));
    double guesstimate_err = 1.0 + cache.A_had_err_low_0[id] + cache.A_had_err_low_1[id] * q2 / 6.0;

    double mB2 = cache.m_Bs * cache.m_Bs;
    double mB3 = cache.m_Bs * mB2;
    double mK2 = cache.m_phi * cache.m_phi;
    double f = lambda(q2) / ((mB2 - mK2) * mB2);
    complex_t delta_A = -N(q2, bar) * cache.m_b_PS * mB2 / (std::sqrt(q2) * cache.m_phi) * (
        (2 * (mB2 + 3 * mK2 - q2) * cache.ff_calculator.E(q2) / mB3 - f) * T_perp_m_cached(q2, bar)
      - f * T_par_m_cached(q2, bar)
    );
    return delta_A * guesstimate_err;
}

complex_t BsPhiDecay::C7_eff(double q2, bool bar) {
    double s_hat = q2 / std::pow(cache.m_b_PS, 2);
    complex_t A = BV::A_Seidel(s_hat, cache.L_b);
    return (bar ? std::conj(cache.C[WCoef::C7]) : cache.C[WCoef::C7]) + cache.alpha_s_mu_b / (4. * PI) * ((cache.C[WCoef::C1]-6.*cache.C[WCoef::C2])*A-cache.C[WCoef::C8]*BV::f_87(q2 / (cache.m_Bs * cache.m_Bs), cache.L_b));
}

complex_t BsPhiDecay::C9_eff(double q2, bool bar) {
    complex_t C_h0 = 4./3.*cache.C[WCoef::C1]+cache.C[WCoef::C2]+11./2.*cache.C[WCoef::C3]-2./3.*cache.C[WCoef::C4]+52.*cache.C[WCoef::C5]-32./3.*cache.C[WCoef::C6];
    complex_t C_hb = -0.5 * (7.*cache.C[WCoef::C3]+4./3.*cache.C[WCoef::C4]+76.*cache.C[WCoef::C5]+64./3.*cache.C[WCoef::C6]);
    complex_t C_0  = 4./3.*(cache.C[WCoef::C3]+16./3.*cache.C[WCoef::C5]+16./9.*cache.C[WCoef::C6]);
    complex_t l_u = bar ? std::conj(cache.lambda_hat_u) : cache.lambda_hat_u;
    complex_t C_mc = 8. * ((4./9.*cache.C[WCoef::C1]+1./3.*cache.C[WCoef::C2])*(1.+l_u)+2.*cache.C[WCoef::C3]+20.*cache.C[WCoef::C5]);

    double s_hat = q2 / std::pow(cache.m_b_PS, 2);
    complex_t A = BV::A_Seidel(s_hat, cache.L_b);
    complex_t B = BV::B_Seidel(s_hat, cache.L_b);
    complex_t C = BV::C_Seidel(q2, cache.mu_b);

    return (bar ? std::conj(cache.C[WCoef::C9]) : cache.C[WCoef::C9])
         + BV::h(q2, 0., cache.mu_b) * C_h0
         + BV::h(q2, cache.m_b_PS, cache.mu_b) * C_hb
         + C_0
         + cache.alpha_s_mu_b / (4. * PI) * (cache.C[WCoef::C1]*(B + 4. * C) - 3. * cache.C[WCoef::C2] * (2. * B - C) - cache.C[WCoef::C8] * BV::f_89(q2 / (cache.m_Bs * cache.m_Bs)))
         + std::pow(cache.m_c_mu_b, 2) / q2 * C_mc;
}

complex_t BsPhiDecay::A_perp_high(double q2, double sign, bool bar) {
    complex_t C7 = C7_eff(q2, bar) + (bar ? std::conj(cache.C[WCoef::CP7]) : cache.C[WCoef::CP7]);
    complex_t C9 = C9_eff(q2, bar) + (bar ? std::conj(cache.C[WCoef::CP9]) : cache.C[WCoef::CP9]);
    complex_t C10 = cache.C[WCoef::C10] + cache.C[WCoef::CP10];
    if (bar) C10 = std::conj(C10);
    return N(q2, bar) * (C9 + sign * C10 + 2. * cache.kappa * cache.m_b_mu_b * cache.m_Bs / q2 * C7) * cache.ff_calculator.get(BV_FF::F_PERP, q2) * (1 + cache.A_had_err_high[size_t (0.5 * (1 + sign))]);
}

complex_t BsPhiDecay::A_par_high(double q2, double sign, bool bar) {
    complex_t C7 = C7_eff(q2, bar) - (bar ? std::conj(cache.C[WCoef::CP7]) : cache.C[WCoef::CP7]);
    complex_t C9 = C9_eff(q2, bar) - (bar ? std::conj(cache.C[WCoef::CP9]) : cache.C[WCoef::CP9]);
    complex_t C10 = cache.C[WCoef::C10] - cache.C[WCoef::CP10];
    if (bar) C10 = std::conj(C10);
    return -N(q2, bar) * (C9 + sign * C10 + 2. * cache.kappa * cache.m_b_mu_b * cache.m_Bs / q2 * C7) * cache.ff_calculator.get(BV_FF::F_PAR, q2) * (1 + cache.A_had_err_high[2 + size_t (0.5 * (1 + sign))]);
}

complex_t BsPhiDecay::A_0_high(double q2, double sign, bool bar) {
    complex_t C7 = C7_eff(q2, bar) - (bar ? std::conj(cache.C[WCoef::CP7]) : cache.C[WCoef::CP7]);
    complex_t C9 = C9_eff(q2, bar) - (bar ? std::conj(cache.C[WCoef::CP9]) : cache.C[WCoef::CP9]);
    complex_t C10 = cache.C[WCoef::C10] - cache.C[WCoef::CP10];
    if (bar) C10 = std::conj(C10);
    return -N(q2, bar) * (C9 + sign * C10 + 2. * cache.kappa * cache.m_b_mu_b * cache.m_Bs / q2 * C7) * cache.ff_calculator.get(BV_FF::F_0, q2)  * (1 + cache.A_had_err_high[4 + size_t (0.5 * (1 + sign))]);
}

complex_t BsPhiDecay::A_t_high(double q2, bool bar) {
    complex_t C10 = cache.C[WCoef::C10] - cache.C[WCoef::CP10];
    complex_t CQ2 = cache.C[WCoef::CQ2] - cache.C[WCoef::CPQ2];
    if (bar) {
        C10 = std::conj(C10);
        CQ2 = std::conj(CQ2);
    }
    return N(q2, bar) / sqrt(q2 * lambda(q2)) * (2. * C10 + q2 / cache.m_l * CQ2 / (cache.m_b_mu_b + cache.m_s)) * cache.ff_calculator.get(BV_FF::A0, q2) * (1 + cache.A_had_err_high[6]);
}

complex_t BsPhiDecay::A_S_high(double q2, bool bar) {
    complex_t CQ1 = cache.C[WCoef::CQ1] - cache.C[WCoef::CPQ1];
    if (bar) CQ1 = std::conj(CQ1);
    return -2. * N(q2, bar) * sqrt(lambda(q2)) * CQ1 / (cache.m_b_mu_b + cache.m_s) * cache.ff_calculator.get(BV_FF::A0, q2) * (1 + cache.A_had_err_high[7]);
}

complex_t BsPhiDecay::interpolate(double q2, complex_t val_low, complex_t val_high) {
    if (q2 < cache.q2_low)
        return val_low;

    if (q2 > cache.q2_high)
        return val_high;

    double t = (cache.q2_high - q2) / (cache.q2_high - cache.q2_low);
    return t * val_low + (1 - t) * val_high;
}

complex_t BsPhiDecay::A_perp(double q2, double sign, bool bar) {
    return interpolate(q2, A_perp_low(q2, sign, !bar), A_perp_high(q2, sign, !bar));
}

complex_t BsPhiDecay::A_par(double q2, double sign, bool bar) {
    return interpolate(q2, A_par_low(q2, sign, !bar), A_par_high(q2, sign, !bar));
}

complex_t BsPhiDecay::A_0(double q2, double sign, bool bar) {
    return interpolate(q2, A_0_low(q2, sign, !bar), A_0_high(q2, sign, !bar));
}

complex_t BsPhiDecay::A_t(double q2, bool bar) {
    return interpolate(q2, A_t_low(q2, !bar), A_t_high(q2, !bar));
}

complex_t BsPhiDecay::A_S(double q2, bool bar) {
    return interpolate(q2, A_S_low(q2, !bar), A_S_high(q2, !bar));
}

double BsPhiDecay::J1s(double q2, bool bar) {
    return (2. + std::pow(beta_l(q2), 2)) / 4. * (
        std::pow(std::abs(A_perp(q2, -1, bar)), 2) 
      + std::pow(std::abs(A_perp(q2, 1, bar)), 2)
      + std::pow(std::abs(A_par(q2, -1, bar)), 2)
      + std::pow(std::abs(A_par(q2, 1, bar)), 2)
    ) + std::pow(2. * cache.m_l, 2) / q2 * std::real(
        A_perp(q2, -1, bar) * std::conj(A_perp(q2, 1, bar))
      + A_par(q2, -1, bar) * std::conj(A_par(q2, 1, bar))
    );
}

double BsPhiDecay::J1c(double q2, bool bar) {
    return std::pow(std::abs(A_0(q2, -1, bar)), 2) + std::pow(std::abs(A_0(q2, 1, bar)), 2) 
         + std::pow(2 * cache.m_l, 2) / q2 * (
              std::pow(std::abs(A_t(q2, bar)), 2)
            + 2. * std::real(A_0(q2, -1, bar) * std::conj(A_0(q2, 1, bar)))
           )
         + std::pow(beta_l(q2) * std::abs(A_S(q2, bar)), 2);
}

double BsPhiDecay::J2s(double q2, bool bar) {
    return std::pow(beta_l(q2), 2) / 4. * (
        std::pow(std::abs(A_perp(q2, -1, bar)), 2) 
      + std::pow(std::abs(A_perp(q2, 1, bar)), 2)
      + std::pow(std::abs(A_par(q2, -1, bar)), 2)
      + std::pow(std::abs(A_par(q2, 1, bar)), 2)
    );
}

double BsPhiDecay::J2c(double q2, bool bar) {
    return -std::pow(beta_l(q2), 2) * (
        std::pow(std::abs(A_0(q2, -1, bar)), 2) 
      + std::pow(std::abs(A_0(q2, 1, bar)), 2)
    );
}

double BsPhiDecay::J3(double q2, bool bar) {
    return std::pow(beta_l(q2), 2) / 2. * (
        std::pow(std::abs(A_perp(q2, -1, bar)), 2) 
      + std::pow(std::abs(A_perp(q2, 1, bar)), 2)
      - std::pow(std::abs(A_par(q2, -1, bar)), 2)
      - std::pow(std::abs(A_par(q2, 1, bar)), 2)
    );
}

double BsPhiDecay::J4(double q2, bool bar) {
    return std::pow(beta_l(q2), 2) / std::sqrt(2.) * (
        std::real(A_0(q2, -1, bar) * std::conj(A_par(q2, -1, bar))) 
      + std::real(A_0(q2, 1, bar) * std::conj(A_par(q2, 1, bar)))
    );
}

double BsPhiDecay::J5(double q2, bool bar) {
    return beta_l(q2) * std::sqrt(2.) * (
        std::real(A_0(q2, -1, bar) * std::conj(A_perp(q2, -1, bar))) 
      - std::real(A_0(q2, 1, bar) * std::conj(A_perp(q2, 1, bar)))
      - cache.m_l / std::sqrt(q2) * std::real((A_par(q2, -1, bar) + A_par(q2, 1, bar)) * std::conj(A_S(q2, bar)))
    );
}

double BsPhiDecay::J6s(double q2, bool bar) {
    return 2. * beta_l(q2) * (
        std::real(A_par(q2, -1, bar) * std::conj(A_perp(q2, -1, bar))) 
      - std::real(A_par(q2, 1, bar) * std::conj(A_perp(q2, 1, bar))) 
    );
}

double BsPhiDecay::J6c(double q2, bool bar) {
    return 4. * beta_l(q2) * cache.m_l / std::sqrt(q2) * (std::real((A_0(q2, -1, bar) + A_0(q2, 1, bar)) * std::conj(A_S(q2, bar))));
}

double BsPhiDecay::J7(double q2, bool bar) {
    return beta_l(q2) * std::sqrt(2.) * (
        std::imag(A_0(q2, -1, bar) * std::conj(A_par(q2, -1, bar))) 
      - std::imag(A_0(q2, 1, bar) * std::conj(A_par(q2, 1, bar)))
      + cache.m_l / std::sqrt(q2) * std::imag((A_perp(q2, -1, bar) + A_perp(q2, 1, bar)) * std::conj(A_S(q2, bar)))
    );
}

double BsPhiDecay::J8(double q2, bool bar) {
    return std::pow(beta_l(q2), 2) / std::sqrt(2.) * (
        std::imag(A_0(q2, -1, bar) * std::conj(A_perp(q2, -1, bar))) 
      + std::imag(A_0(q2, 1, bar) * std::conj(A_perp(q2, 1, bar)))
    );
}

double BsPhiDecay::J9(double q2, bool bar) {
    return std::pow(beta_l(q2), 2) * (
        std::imag(A_perp(q2, -1, bar) * std::conj(A_par(q2, -1, bar))) 
      + std::imag(A_perp(q2, 1, bar) * std::conj(A_par(q2, 1, bar))) 
    );
}

double BsPhiDecay::h1s(double q2) {
    complex_t ALperp = A_perp(q2, -1, false);
    complex_t ARperp = A_perp(q2, 1, false);
    complex_t ALpar = A_par(q2, -1, false);
    complex_t ARpar = A_par(q2, 1, false);
    complex_t ALperp_tilde = -A_perp(q2, -1, true);
    complex_t ARperp_tilde = -A_perp(q2, 1, true);
    complex_t ALpar_tilde = A_par(q2, -1, true);
    complex_t ARpar_tilde = A_par(q2, 1, true);
    complex_t sh1s_A = cache.up * (ALperp_tilde * std::conj(ALperp) + ALpar_tilde * std::conj(ALpar) + ARperp_tilde * std::conj(ARperp) + ARpar_tilde * std::conj(ARpar));
	complex_t sh1s_B = cache.up * (ALperp_tilde * std::conj(ARperp) + ALpar_tilde * std::conj(ARpar));
	complex_t sh1s_C = cache.um * (ALperp * std::conj(ARperp_tilde) + ALpar * std::conj(ARpar_tilde));
    double x = std::pow(cache.m_l, 2) / q2;

    return (2. + std::pow(beta_l(q2), 2)) / 2. * std::real(sh1s_A) + 4. * x * std::real(sh1s_B + sh1s_C);
}

double BsPhiDecay::h1c(double q2) {
    complex_t AL0 = A_0(q2, -1, false);
    complex_t AR0 = A_0(q2, 1, false);
    complex_t AL0_tilde = A_0(q2, -1, true);
    complex_t AR0_tilde = A_0(q2, 1, true);
    complex_t sh1c_A = cache.up * (AL0_tilde * std::conj(AL0) + AR0_tilde * std::conj(AR0));
	complex_t sh1c_B = cache.up * A_t(q2, true) * std::conj(A_t(q2, false));
	complex_t sh1c_C = cache.up * AL0_tilde * std::conj(AR0);
	complex_t sh1c_D = cache.um * AL0 * std::conj(AR0_tilde);
	complex_t sh1c_E = cache.up * -A_S(q2, true) * std::conj(A_S(q2, false));
    double x = std::pow(cache.m_l, 2) / q2;

    return 2. * std::real(sh1c_A) + 8. * x * (std::real(sh1c_B) + std::real(sh1c_C + sh1c_D)) + 2. * std::pow(beta_l(q2), 2) * std::real(sh1c_E);
}

double BsPhiDecay::h2s(double q2) {
    complex_t sh1s_A = cache.up * (-A_perp(q2, -1, true) * std::conj(A_perp(q2, -1, false)) + A_par(q2, -1, true) * std::conj(A_par(q2, -1, false)) - A_perp(q2, 1, true) * std::conj(A_perp(q2, 1, false)) + A_par(q2, 1, true) * std::conj(A_par(q2, 1, false)));
    return std::pow(beta_l(q2), 2) / 2. * std::real(sh1s_A);
}

double BsPhiDecay::h2c(double q2) {
    complex_t sh1c_A = cache.up * (A_0(q2, -1, true) * std::conj(A_0(q2, -1, false)) + A_0(q2, 1, true) * std::conj(A_0(q2, 1, false)));
    return -2. * std::pow(beta_l(q2), 2) * std::real(sh1c_A);
}

double BsPhiDecay::h3(double q2) {
    complex_t sh3_A = -cache.up * (A_perp(q2, -1, true) * std::conj(A_perp(q2, -1, false)) + A_par(q2, -1, true) * std::conj(A_par(q2, -1, false)) + A_perp(q2, 1, true) * std::conj(A_perp(q2, 1, false)) + A_par(q2, 1, true) * std::conj(A_par(q2, 1, false)));
    return std::pow(beta_l(q2), 2) * std::real(sh3_A);
}

double BsPhiDecay::h4(double q2) {
    complex_t sh4_A  = cache.up * (A_0(q2, -1, true) * std::conj(A_par(q2, -1, false)) + A_0(q2, 1, true) * std::conj(A_par(q2, 1, false)));
	complex_t sh4_B  = cache.um * (A_0(q2, -1, false) * std::conj(A_par(q2, -1, true)) + A_0(q2, 1, false) * std::conj(A_par(q2, 1, true)));
    return INV_RT2 * std::pow(beta_l(q2), 2) * std::real(sh4_A + sh4_B);
}

double BsPhiDecay::h5(double q2) {
    complex_t sh5_A  = cache.up * (A_0(q2, -1, true) * std::conj(A_perp(q2, -1, false)) - A_0(q2, 1, true) * std::conj(A_perp(q2, 1, false)));
	complex_t sh5_B  = cache.um * (A_0(q2, -1, false) * std::conj(-A_perp(q2, -1, true)) - A_0(q2, 1, false) * std::conj(A_par(q2, 1, true)));
	complex_t sh5_C  = cache.up * ((A_par(q2, -1, true) + A_par(q2, 1, true)) * std::conj(A_S(q2, false)));
	complex_t sh5_D  = cache.um * ((A_par(q2, -1, false) + A_par(q2, 1, false)) * std::conj(-A_S(q2, true)));
    return RT2 * beta_l(q2) * (std::real(sh5_A + sh5_B)  - cache.m_l / std::sqrt(q2) * std::real(sh5_C + sh5_D));
}

double BsPhiDecay::h6s(double q2) {
    complex_t sh6s_A = cache.up * (A_par(q2, -1, true) * std::conj(A_perp(q2, -1, false)) - A_par(q2, 1, true) * std::conj(A_perp(q2, 1, false)));
	complex_t sh6s_B = cache.um * (A_par(q2, -1, false) * std::conj(-A_perp(q2, -1, true)) - A_par(q2, 1, false) * std::conj(-A_perp(q2, 1, true)));
    return 2. * beta_l(q2) * std::real(sh6s_A + sh6s_B);
}

double BsPhiDecay::h6c(double q2) {
    complex_t sh6c_A = cache.up * ((A_0(q2, -1, true) + A_0(q2, 1, true)) * std::conj(A_S(q2, false)));
	complex_t sh6c_B = cache.um * ((A_0(q2, -1, false) + A_0(q2, 1, false)) * std::conj(-A_S(q2, true)));
    return 4. * beta_l(q2) * cache.m_l / std::sqrt(q2) * std::real(sh6c_A + sh6c_B);
}

double BsPhiDecay::h7(double q2) {
    complex_t sh7_A  = cache.up * (A_0(q2, -1, true) * std::conj(A_par(q2, -1, false)) - A_0(q2, 1, true) * std::conj(A_par(q2, 1, false)));
	complex_t sh7_B  = cache.um * (A_0(q2, -1, false) * std::conj(A_par(q2, -1, true)) - A_0(q2, 1, false) * std::conj(A_par(q2, 1, true)));
	complex_t sh7_C  = cache.up * (-(A_perp(q2, -1, true) + A_perp(q2, 1, true)) * std::conj(A_S(q2, false)));
	complex_t sh7_D  = cache.um * ((A_perp(q2, -1, false) + A_perp(q2, 1, false)) * std::conj(-A_S(q2, true)));
    return sqrt(2.) * beta_l(q2) * (std::imag(sh7_A + sh7_B) + cache.m_l / std::sqrt(q2) * std::imag(sh7_C + sh7_D));
}

double BsPhiDecay::h8(double q2) {
    complex_t sh8_A  = cache.up * (A_0(q2, -1, true) * std::conj(A_perp(q2, -1, false)) + A_0(q2, 1, true) * std::conj(A_perp(q2, 1, false)));
	complex_t sh8_B  = cache.um * (A_0(q2, -1, false) * std::conj(-A_perp(q2, -1, true)) + A_0(q2, 1, false) * std::conj(-A_perp(q2, 1, true)));
    return 1. / sqrt(q2) * std::pow(beta_l(q2), 2) * std::imag(sh8_A + sh8_B);
}

double BsPhiDecay::h9(double q2) {
    complex_t sh9_A  = cache.up * (A_par(q2, -1, true) * std::conj(A_perp(q2, -1, false)) + A_par(q2, 1, true) * std::conj(A_perp(q2, 1, false)));
	complex_t sh9_B  = cache.um * (A_par(q2, -1, false) * std::conj(-A_perp(q2, -1, true)) + A_par(q2, 1, false) * std::conj(-A_perp(q2, 1, true)));
    return -std::pow(beta_l(q2), 2) * std::imag(sh9_A + sh9_B);
}

double BsPhiDecay::s8(double q2) {
    complex_t sh8_A  = cache.up * (A_0(q2, -1, true) * std::conj(A_perp(q2, -1, false)) + A_0(q2, 1, true) * std::conj(A_perp(q2, 1, false)));
	complex_t sh8_B  = cache.um * (A_0(q2, -1, false) * std::conj(-A_perp(q2, -1, true)) + A_0(q2, 1, false) * std::conj(-A_perp(q2, 1, true)));
    return -1. / sqrt(q2) * std::pow(beta_l(q2), 2) * std::real(sh8_A - sh8_B);
}

double BsPhiDecay::s9(double q2) {
    complex_t sh9_A  = cache.up * (A_par(q2, -1, true) * std::conj(A_perp(q2, -1, false)) + A_par(q2, 1, true) * std::conj(A_perp(q2, 1, false)));
	complex_t sh9_B  = cache.um * (A_par(q2, -1, false) * std::conj(-A_perp(q2, -1, true)) + A_par(q2, 1, false) * std::conj(-A_perp(q2, 1, true)));
    return std::pow(beta_l(q2), 2) * std::real(sh9_A - sh9_B);
}

void BsPhiDecay::compute_binned_J_i() {
    auto zero_if_close = [] (double x, double tol) {
        return std::abs(x) < tol ? 0.0 : x;
    };

    for (auto [q2_l, q2_u] : cfg.bins) {
        cache.f_J_i_binned[ 0].emplace_back(integrate([&] (double q2) { return 2 * (J1s(q2, false) - cache.ys * h1s(q2) / 2) + J1c(q2, false) - cache.ys * h1c(q2) / 2 - (2 * (J2s(q2, false) - cache.ys * h2s(q2) / 2) + J2c(q2, false) - cache.ys * h2c(q2) / 2) / 3.; }, q2_l, q2_u, 1e-2));    
        cache.f_J_i_binned[ 1].emplace_back(integrate([&] (double q2) { return 2 * (J1s(q2, true ) - cache.ys * h1s(q2) / 2) + J1c(q2, true ) - cache.ys * h1c(q2) / 2 - (2 * (J2s(q2, true ) - cache.ys * h2s(q2) / 2) + J2c(q2, true ) - cache.ys * h2c(q2) / 2) / 3.; }, q2_l, q2_u, 1e-2));
        cache.f_J_i_binned[ 2].emplace_back(integrate([&] (double q2) { return J2s(q2, true) + J2s(q2, false) - cache.ys * h2s(q2); },  q2_l, q2_u, 1e-2));
        cache.f_J_i_binned[ 3].emplace_back(integrate([&] (double q2) { return J2c(q2, true) + J2c(q2, false) - cache.ys * h2c(q2); },  q2_l, q2_u, 1e-2));
        cache.f_J_i_binned[ 4].emplace_back(integrate([&] (double q2) { return J3(q2, true) + J3(q2, false) - cache.ys * h3(q2); },  q2_l, q2_u, 1e-2));
        cache.f_J_i_binned[ 5].emplace_back(integrate([&] (double q2) { return J4(q2, true) + J4(q2, false) - cache.ys * h4(q2); },  q2_l, q2_u, 1e-2));
        cache.f_J_i_binned[ 6].emplace_back(integrate([&] (double q2) { return zero_if_close(J5(q2, false) - J5(q2, true) - cache.ys * h5(q2), 1e-30); },  q2_l, q2_u, 1e-2));
        cache.f_J_i_binned[ 7].emplace_back(integrate([&] (double q2) { return zero_if_close(J6s(q2, false) - J6s(q2, true) - cache.ys * h6s(q2), 1e-30); },  q2_l, q2_u, 1e-2));
        cache.f_J_i_binned[ 8].emplace_back(integrate([&] (double q2) { return beta_l(q2) * (zero_if_close(J6s(q2, false) - J6s(q2, true) - cache.ys * h6s(q2), 1e-30)); }, q2_l, q2_u, 1e-2));
        cache.f_J_i_binned[ 9].emplace_back(integrate([&] (double q2) { return zero_if_close(J6c(q2, false) - J6c(q2, true) - cache.ys * h6c(q2), 1e-30); },  q2_l, q2_u, 1e-2));
        cache.f_J_i_binned[10].emplace_back(integrate([&] (double q2) { return J7 (q2, true) + J7 (q2, false) - cache.ys * h7(q2); },  q2_l, q2_u, 1e-2));
        cache.f_J_i_binned[11].emplace_back(integrate([&] (double q2) { return zero_if_close(J8(q2, false) - J8(q2, true) - cache.ys * h8(q2), 1e-30); },  q2_l, q2_u, 1e-2));
        cache.f_J_i_binned[12].emplace_back(integrate([&] (double q2) { return zero_if_close(J9(q2, false) - J9(q2, true) - cache.ys * h9(q2), 1e-30); },  q2_l, q2_u, 1e-2));
        cache.f_J_i_binned[13].emplace_back(integrate([&] (double q2) { return s8(q2); }, q2_l, q2_u, 1e-2));
        cache.f_J_i_binned[14].emplace_back(integrate([&] (double q2) { return s9(q2); }, q2_l, q2_u, 1e-2));
    }
}

std::vector<ObservableValue> BsPhiDecay::dG_dq2_binned(bool bar) {
    std::vector<ObservableValue> out;
    size_t idx = bar ? 1 : 0;
    ObservableId id = bar ? ObservableMapper::to_id(Observables::DGAMMA_BAR_DQ2_BS__PHI_L_L) : ObservableMapper::to_id(Observables::DGAMMA_DQ2_BS__PHI_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double res = 0.75 * cache.f_J_i_binned[idx][i]; 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

double BsPhiDecay::dG_dq2_avg_bin(size_t bin) {
    return 0.75 * (cache.f_J_i_binned[0][bin] + cache.f_J_i_binned[1][bin]);
}

std::vector<ObservableValue> BsPhiDecay::F_L() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::F_L_BS_PHI_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double res = -cache.f_J_i_binned[3][i] / dG_dq2_avg_bin(i); 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::A_T_2() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::A_T_2_BS_PHI_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double res = 0.5 * cache.f_J_i_binned[4][i] / cache.f_J_i_binned[2][i]; 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::A_T_Re_CPV() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::A_T_RE_CPV_BS_PHI_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double res = 0.25 * cache.f_J_i_binned[8][i] / cache.f_J_i_binned[2][i]; 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::A_T_Im_CPV() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::A_T_IM_CPV_BS_PHI_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double res = 0.5 * cache.f_J_i_binned[12][i] / cache.f_J_i_binned[2][i]; 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::Pp_4() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::P_PRIME_4_BS_PHI_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double res = cache.f_J_i_binned[5][i] / std::sqrt(std::abs(cache.f_J_i_binned[2][i] * cache.f_J_i_binned[3][i])); 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::Pp_6() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::P_PRIME_6_BS_PHI_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double res = -0.5 * cache.f_J_i_binned[10][i] / std::sqrt(std::abs(cache.f_J_i_binned[2][i] * cache.f_J_i_binned[3][i])); 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::S_i(int i) {
    if (!(i == 2 || i == 3 || i == 4 || i == 7)) LOG_ERROR("Value Error", "S_i(Bs > phi ll) is not defined for i =", i);

    std::map<size_t, size_t> J_idx = {{2, 2}, {3, 4}, {4, 5}, {7, 10}};
    std::map<size_t, Observables> ids = {
        {2, Observables::S_2S_BS_PHI_L_L},
        {3, Observables::S_3_BS_PHI_L_L},
        {4, Observables::S_4_BS_PHI_L_L},
        {7, Observables::S_7_BS_PHI_L_L}
    };

    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(ids[i]);
    for (size_t j = 0; j < cfg.bins.size(); j++) {
        double res = cache.f_J_i_binned[J_idx[i]][j] / dG_dq2_avg_bin(j);
        out.emplace_back(id, res, cfg.bins[j]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::A_i(int i) {
    if (!(i == 5 || i == 6 || i == 8 || i == 9)) LOG_ERROR("Value Error", "A_i(Bs > phi ll) is not defined for i =", i);

    std::map<size_t, size_t> J_idx = {{5, 6}, {6, 9}, {8, 11}, {9, 12}};
    std::map<size_t, Observables> ids = {
        {5, Observables::A_5_BS_PHI_L_L},
        {6, Observables::A_6C_BS_PHI_L_L},
        {8, Observables::A_8_BS_PHI_L_L},
        {9, Observables::A_9_BS_PHI_L_L}
    };
    double sign = i == 6 ? -1 : 1;

    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(ids[i]);
    for (size_t j = 0; j < cfg.bins.size(); j++) {
        double res = sign * cache.f_J_i_binned[J_idx[i]][j] / dG_dq2_avg_bin(j);
        out.emplace_back(id, res, cfg.bins[j]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::A_FB_CPV() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::A_FB_CPV_BS__PHI_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double res = -0.375 * (2 * cache.f_J_i_binned[7][i] + cache.f_J_i_binned[9][i]) / dG_dq2_avg_bin(i); 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::P_2_CPV() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::P_2_CPV_BS_PHI_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double res = 0.125 * cache.f_J_i_binned[7][i] / cache.f_J_i_binned[2][i]; 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::P_3_CPV() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::P_3_CPV_BS_PHI_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double res = -0.25 * cache.f_J_i_binned[12][i] / cache.f_J_i_binned[2][i]; 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::Pp_5_CPV() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::P_PRIME_5_CPV_BS_PHI_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double res = 0.5 * cache.f_J_i_binned[6][i] / std::sqrt(std::abs(cache.f_J_i_binned[2][i] * cache.f_J_i_binned[3][i])); 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::Pp_8_CPV() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::P_PRIME_8_CPV_BS_PHI_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double res = -cache.f_J_i_binned[11][i] / std::sqrt(std::abs(cache.f_J_i_binned[2][i] * cache.f_J_i_binned[3][i])); 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::Q_8_m() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::Q_8M_BS_PHI_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double res = cache.f_J_i_binned[13][i] / std::sqrt(std::abs(2 * cache.f_J_i_binned[3][i] * (cache.f_J_i_binned[2][i] - cache.f_J_i_binned[4][i]))); 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::Q_8_p() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::Q_8P_BS_PHI_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double res = cache.f_J_i_binned[13][i] / std::sqrt(std::abs(2 * cache.f_J_i_binned[3][i] * (cache.f_J_i_binned[2][i] + cache.f_J_i_binned[4][i]))); 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::Q_9() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::Q_9_BS_PHI_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double res = 0.5 * cache.f_J_i_binned[14][i] / cache.f_J_i_binned[2][i]; 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::compute_observable(Observables obs) {
    switch (obs) {
    case Observables::TEST:   
        // test_binned_obs();
        return {};
    case Observables::DGAMMA_DQ2_BS__PHI_L_L:   
        return dG_dq2_binned(false);
    case Observables::DGAMMA_BAR_DQ2_BS__PHI_L_L:   
        return dG_dq2_binned(true);
    case Observables::F_L_BS_PHI_L_L:   
        return F_L();
    case Observables::A_T_2_BS_PHI_L_L:   
        return A_T_2();
    case Observables::A_T_RE_CPV_BS_PHI_L_L:   
        return A_T_Re_CPV();
    case Observables::A_T_IM_CPV_BS_PHI_L_L:   
        return A_T_Im_CPV();
    case Observables::P_PRIME_4_BS_PHI_L_L:   
        return Pp_4();
    case Observables::P_PRIME_6_BS_PHI_L_L:   
        return Pp_6();
    case Observables::S_2S_BS_PHI_L_L:   
        return S_i(2);
    case Observables::S_3_BS_PHI_L_L:   
        return S_i(3);
    case Observables::S_4_BS_PHI_L_L:   
        return S_i(4);
    case Observables::S_7_BS_PHI_L_L:   
        return S_i(7);
    case Observables::A_5_BS_PHI_L_L:   
        return A_i(5);
    case Observables::A_6C_BS_PHI_L_L:   
        return A_i(6);
    case Observables::A_8_BS_PHI_L_L:   
        return A_i(8);
    case Observables::A_9_BS_PHI_L_L:   
        return A_i(9);
    case Observables::A_FB_CPV_BS__PHI_L_L:   
        return A_FB_CPV();
    case Observables::P_2_CPV_BS_PHI_L_L:   
        return P_2_CPV();
    case Observables::P_3_CPV_BS_PHI_L_L:   
        return P_3_CPV();
    case Observables::P_PRIME_5_CPV_BS_PHI_L_L:   
        return Pp_5_CPV();
    case Observables::P_PRIME_8_CPV_BS_PHI_L_L:   
        return Pp_8_CPV();
    case Observables::Q_8M_BS_PHI_L_L:   
        return Q_8_m();
    case Observables::Q_8P_BS_PHI_L_L:   
        return Q_8_p();
    case Observables::Q_9_BS_PHI_L_L:   
        return Q_9();
    default:
        LOG_ERROR("IndexError", "Observable", ObservableMapper::str(obs), "doesn't belong to the decay", DecayMapper::str(this->id));
    }
}

std::vector<ObservableValue> BsPhiDecay::compute_observable(ObservableId obs) {
    return compute_observable(ObservableMapper::enum_of(obs).value());
}