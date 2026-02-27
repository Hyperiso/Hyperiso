#include "LbLllDecay.h"

void LbLllDecay::load_params() {
    fill_wilson_cache();

    cache.ff_calculator = LbLFFCalculator(p, cfg.ff_src);
    cache.m_b_mu_b = (*iobs_qcdp)(MassConfig(5, w_config.hadronic_scale, MassType::MSBAR, MassType::POLE));
    cache.m_Lb = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", 5122}, DataType::VALUE);
    cache.m_L = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", 3122}, DataType::VALUE);
    cache.life_L = (*p)(ParamId{ParameterType::FLAVOR, "FLIFE", 5122}, DataType::VALUE);
    cache.alpha_L = (*p)(ParamId{ParameterType::DECAY, "Lb_L", 8}, DataType::VALUE);
    cache.N_0 = std::conj((*p)(ParamId{ParameterType::SM, "VCKM", {2, 1}}, DataType::VALUE)) * (*p)(ParamId{ParameterType::SM, "VCKM", {2, 2}}, DataType::VALUE) * (*p)(ParamId{ParameterType::SM, "SMINPUTS", 2}, DataType::VALUE) * (*p)(ParamId{ParameterType::SM, "EW", {1, 2}}, DataType::VALUE) / (std::sqrt(6144. * std::pow(PI, 5) * std::pow(cache.m_Lb, 3)));
    cache.q2_max = std::pow(cache.m_Lb - cache.m_L, 2);

    load_cfg_dep_params();

    // printf("alpha_em = %.4e\n", (*p)(ParamId{ParameterType::SM, "EW", {1, 2}}, DataType::VALUE).real());
    // printf("N0 = %.4e + %.4e i\n", cache.N_0.real(), cache.N_0.imag());

    // printf("f_perp = %.4e\n", cache.ff_calculator.get(LbL_FF::F_PERP, 1.0));
    // printf("h_perp = %.4e\n", cache.ff_calculator.get(LbL_FF::H_PERP, 1.0));
    // printf("g_perp = %.4e\n", cache.ff_calculator.get(LbL_FF::G_PERP, 1.0));
    // printf("f_+ = %.4e\n", cache.ff_calculator.get(LbL_FF::F_PLUS, 1.0));
    // printf("h_+ = %.4e\n", cache.ff_calculator.get(LbL_FF::H_PLUS, 1.0));
    // printf("g_+ = %.4e\n", cache.ff_calculator.get(LbL_FF::G_PLUS, 1.0));
    // printf("h_tilde_perp = %.4e\n", cache.ff_calculator.get(LbL_FF::H_TILDE_PERP, 1.0));
    // printf("h_tilde_+ = %.4e\n", cache.ff_calculator.get(LbL_FF::H_TILDE_PLUS, 1.0));
}

void LbLllDecay::fill_wilson_cache() {
    auto b_wilsons = w_proxy->getAFR(WGroup::B, this->w_config.order);
    auto bp_wilsons = w_proxy->getAFR(WGroup::BPrime, this->w_config.order);

    cache.C.emplace(WCoef::C7, b_wilsons[WCoef::C7]);
    cache.C.emplace(WCoef::C9, b_wilsons[WCoef::C9]); 
    cache.C.emplace(WCoef::C10, b_wilsons[WCoef::C10]);
    cache.C.emplace(WCoef::CP7, b_wilsons[WCoef::CP7]);
    cache.C.emplace(WCoef::CP9, b_wilsons[WCoef::CP9]);
    cache.C.emplace(WCoef::CP10, b_wilsons[WCoef::CP10]);
}

void LbLllDecay::set_cfg_flags(LbLllConfig::Lepton gen) {
    if (cfg.gen != gen) {
        cfg.gen = gen;
        load_cfg_dep_params();
    }
}

void LbLllDecay::load_cfg_dep_params() {
    cache.m_l = (*p)(ParamId{ParameterType::SM, "MASS", 11 + 2 * (int)cfg.gen}, DataType::VALUE);
    cache.q2_min = 4 * std::pow(cache.m_l, 2);
    compute_binned_K_i();
}

double LbLllDecay::beta_l(double q2) {
    return std::sqrt(1 - 4. * cache.m_l * cache.m_l / q2);
}

double LbLllDecay::lambda(double q2) {
    double mLb2 = cache.m_Lb * cache.m_Lb;
    double mL2 = cache.m_L * cache.m_L;
    return mLb2 * mLb2 + mL2 * mL2 + q2 * q2 - 2. * (mLb2 * mL2 + (mLb2 + mL2) * q2);
}

double LbLllDecay::s_p(double q2) {
    return std::pow(cache.m_Lb + cache.m_L, 2) - q2;
}

double LbLllDecay::s_m(double q2) {
    return std::pow(cache.m_Lb - cache.m_L, 2) - q2;
}

complex_t LbLllDecay::N(double q2, bool bar) {
    complex_t N0 = bar ? std::conj(cache.N_0) : cache.N_0; 
    return N0 * std::sqrt(q2 * beta_l(q2) * std::sqrt(lambda(q2)));
}

complex_t LbLllDecay::A_perp_1(double q2, double sign, bool bar) {
    complex_t C_V = (cache.C[WCoef::C9] + cache.C[WCoef::CP9]) + sign * (cache.C[WCoef::C10] + cache.C[WCoef::CP10]);
    complex_t C_T = cache.C[WCoef::C7] + cache.C[WCoef::CP7];
    if (bar) {
        C_V = std::conj(C_V);
        C_T = std::conj(C_T);
    }
    complex_t HVplus = -cache.ff_calculator.get(LbL_FF::F_PERP, q2) * std::sqrt(2 * s_m(q2));
    complex_t HTplus = cache.ff_calculator.get(LbL_FF::H_PERP, q2) * (cache.m_Lb + cache.m_L) * std::sqrt(2 * s_m(q2));
    return RT2 * N(q2, bar) * (C_V * HVplus - 2 * cache.m_b_mu_b / q2 * C_T * HTplus);
}

complex_t LbLllDecay::A_par_1(double q2, double sign, bool bar) {
    complex_t C_A = (cache.C[WCoef::C9] - cache.C[WCoef::CP9]) + sign * (cache.C[WCoef::C10] - cache.C[WCoef::CP10]);
    complex_t C_TA = cache.C[WCoef::C7] - cache.C[WCoef::CP7];
    if (bar) {
        C_A = std::conj(C_A);
        C_TA = std::conj(C_TA);
    }
    complex_t HAplus = -cache.ff_calculator.get(LbL_FF::G_PERP, q2) * std::sqrt(2 * s_p(q2));
    complex_t HT5plus = -cache.ff_calculator.get(LbL_FF::H_TILDE_PERP, q2) * (cache.m_Lb - cache.m_L) * std::sqrt(2 * s_p(q2));
    return -RT2 * N(q2, bar) * (C_A * HAplus + 2 * cache.m_b_mu_b / q2 * C_TA * HT5plus);
}

complex_t LbLllDecay::A_perp_0(double q2, double sign, bool bar) {
    complex_t C_V = (cache.C[WCoef::C9] + cache.C[WCoef::CP9]) + sign * (cache.C[WCoef::C10] + cache.C[WCoef::CP10]);
    complex_t C_T = cache.C[WCoef::C7] + cache.C[WCoef::CP7];
    if (bar) {
        C_V = std::conj(C_V);
        C_T = std::conj(C_T);
    }
    complex_t HV0 = cache.ff_calculator.get(LbL_FF::F_PLUS, q2) * (cache.m_Lb + cache.m_L) * std::sqrt(s_m(q2) / q2);
    complex_t HT0 = -cache.ff_calculator.get(LbL_FF::H_PLUS, q2) * std::sqrt(q2 * s_m(q2));
    return RT2 * N(q2, bar) * (C_V * HV0 - 2 * cache.m_b_mu_b / q2 * C_T * HT0);
}

complex_t LbLllDecay::A_par_0(double q2, double sign, bool bar) {
    complex_t C_A = (cache.C[WCoef::C9] - cache.C[WCoef::CP9]) + sign * (cache.C[WCoef::C10] - cache.C[WCoef::CP10]);
    complex_t C_TA = cache.C[WCoef::C7] - cache.C[WCoef::CP7];
    if (bar) {
        C_A = std::conj(C_A);
        C_TA = std::conj(C_TA);
    }
    complex_t HA0 = cache.ff_calculator.get(LbL_FF::G_PLUS, q2) * (cache.m_Lb - cache.m_L) * std::sqrt(s_p(q2) / q2);
    complex_t HT50 = cache.ff_calculator.get(LbL_FF::H_TILDE_PLUS, q2) * std::sqrt(q2 * s_p(q2));
    return -RT2 * N(q2, bar) * (C_A * HA0 + 2 * cache.m_b_mu_b / q2 * C_TA * HT50);
}

double LbLllDecay::K1ss(double q2, bool bar) {
    complex_t ALpar1 = A_par_1(q2, -1, bar);
    complex_t ALperp1 = A_perp_1(q2, -1, bar);
    complex_t ARpar1 = A_par_1(q2, 1, bar);
    complex_t ARperp1 = A_perp_1(q2, 1, bar);
    complex_t ALpar0 = A_par_0(q2, -1, bar);
    complex_t ALperp0 = A_perp_0(q2, -1, bar);
    complex_t ARpar0 = A_par_0(q2, 1, bar);
    complex_t ARperp0 = A_perp_0(q2, 1, bar);

    double b_l = beta_l(q2);
    
    return std::real(0.25*(ALpar1*conj(ALpar1)+ALperp1*conj(ALperp1)+ARpar1*conj(ARpar1)+ARperp1*conj(ARperp1))
         + 0.25*(1.+b_l*b_l)*(ALpar0*conj(ALpar0)+ALperp0*conj(ALperp0)+ARpar0*conj(ARpar0)+ARperp0*conj(ARperp0))
         + 0.5*(1.-b_l*b_l)*(ARpar1*conj(ALpar1)+ARperp1*conj(ALperp1)+ARpar0*conj(ALpar0)+ARperp0*conj(ALperp0)));
}

double LbLllDecay::K1cc(double q2, bool bar) {
    complex_t ALpar1 = A_par_1(q2, -1, bar);
    complex_t ALperp1 = A_perp_1(q2, -1, bar);
    complex_t ARpar1 = A_par_1(q2, 1, bar);
    complex_t ARperp1 = A_perp_1(q2, 1, bar);
    complex_t ALpar0 = A_par_0(q2, -1, bar);
    complex_t ALperp0 = A_perp_0(q2, -1, bar);
    complex_t ARpar0 = A_par_0(q2, 1, bar);
    complex_t ARperp0 = A_perp_0(q2, 1, bar);

    double b_l = beta_l(q2);
    
    return std::real(0.25*(1.+b_l*b_l)*(ARpar1*conj(ARpar1)+ARperp1*conj(ARperp1)+ALpar1*conj(ALpar1)+ALperp1*conj(ALperp1))
					+0.25*(1.-b_l*b_l)*(ARpar0*conj(ARpar0)+ARperp0*conj(ARperp0)+ALpar0*conj(ALpar0)+ALperp0*conj(ALperp0))
					+0.5*(1.-b_l*b_l)*(ARpar1*conj(ALpar1)+ARperp1*conj(ALperp1)+ARpar0*conj(ALpar0)+ARperp0*conj(ALperp0)));
}

double LbLllDecay::K1c(double q2, bool bar) {
    complex_t ALpar1 = A_par_1(q2, -1, bar);
    complex_t ALperp1 = A_perp_1(q2, -1, bar);
    complex_t ARpar1 = A_par_1(q2, 1, bar);
    complex_t ARperp1 = A_perp_1(q2, 1, bar);
    double b_l = beta_l(q2);
    return -b_l*std::real(ARperp1*conj(ARpar1)-ALperp1*conj(ALpar1));;
}

double LbLllDecay::K2ss(double q2, bool bar) {
    complex_t ALpar1 = A_par_1(q2, -1, bar);
    complex_t ALperp1 = A_perp_1(q2, -1, bar);
    complex_t ARpar1 = A_par_1(q2, 1, bar);
    complex_t ARperp1 = A_perp_1(q2, 1, bar);
    complex_t ALpar0 = A_par_0(q2, -1, bar);
    complex_t ALperp0 = A_perp_0(q2, -1, bar);
    complex_t ARpar0 = A_par_0(q2, 1, bar);
    complex_t ARperp0 = A_perp_0(q2, 1, bar);

    double b_l = beta_l(q2);
    
    return 0.5 * cache.alpha_L * std::real(
        ARperp1*conj(ARpar1)+ALperp1*conj(ALpar1)
      + (1. + b_l*b_l)*(ARperp0*conj(ARpar0)+ALperp0*conj(ALpar0))
      + (1. - b_l*b_l)*(ARperp1*conj(ALpar1)+ARpar1*conj(ALperp1)+ARperp0*conj(ALpar0)+ARpar0*conj(ALperp0)));
}

double LbLllDecay::K2cc(double q2, bool bar) {
    complex_t ALpar1 = A_par_1(q2, -1, bar);
    complex_t ALperp1 = A_perp_1(q2, -1, bar);
    complex_t ARpar1 = A_par_1(q2, 1, bar);
    complex_t ARperp1 = A_perp_1(q2, 1, bar);
    complex_t ALpar0 = A_par_0(q2, -1, bar);
    complex_t ALperp0 = A_perp_0(q2, -1, bar);
    complex_t ARpar0 = A_par_0(q2, 1, bar);
    complex_t ARperp0 = A_perp_0(q2, 1, bar);

    double b_l = beta_l(q2);
    
    return 0.5 * cache.alpha_L * std::real(
        (1.+b_l*b_l)*(ARperp1*conj(ARpar1)+ALperp1*conj(ALpar1))
      + (1.-b_l*b_l)*(ARpar0*conj(ARperp0)+ALpar0*conj(ALperp0))
      + (1.-b_l*b_l)*(ARperp1*conj(ALpar1)+ARpar1*conj(ALperp1)+ARperp0*conj(ALpar0)+ARpar0*conj(ALperp0)));
}

double LbLllDecay::K2c(double q2, bool bar) {
    complex_t ALpar1 = A_par_1(q2, -1, bar);
    complex_t ALperp1 = A_perp_1(q2, -1, bar);
    complex_t ARpar1 = A_par_1(q2, 1, bar);
    complex_t ARperp1 = A_perp_1(q2, 1, bar);
    double b_l = beta_l(q2);
    return -0.5*cache.alpha_L*b_l*std::real(ARpar1*conj(ARpar1)+ARperp1*conj(ARperp1)-ALpar1*conj(ALpar1)-ALperp1*conj(ALperp1));
}

void LbLllDecay::compute_binned_K_i() {
    auto fill_binned = [&] (std::array<std::vector<double>, 6>& dest, bool bar) {
        for (auto [q2_l, q2_u] : this->bins.value()) {
            dest[0].emplace_back(integrate([&] (double q2) { return K1ss(q2, bar); }, q2_l, q2_u, 1e-2));
            dest[1].emplace_back(integrate([&] (double q2) { return K1cc(q2, bar); }, q2_l, q2_u, 1e-2));
            dest[2].emplace_back(integrate([&] (double q2) { return K1c(q2, bar); }, q2_l, q2_u, 1e-2));
            dest[3].emplace_back(integrate([&] (double q2) { return K2ss(q2, bar); }, q2_l, q2_u, 1e-2));
            dest[4].emplace_back(integrate([&] (double q2) { return K2cc(q2, bar); }, q2_l, q2_u, 1e-2));
            dest[5].emplace_back(integrate([&] (double q2) { return K2c(q2, bar); }, q2_l, q2_u, 1e-2));
        }
    };

    fill_binned(cache.K_i_binned, false);
    fill_binned(cache.K_i_bar_binned, true);
}

std::vector<ObservableValue> LbLllDecay::dBR_dq2_binned(Observables oid) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double K1ss = cache.K_i_binned[0][i] + cache.K_i_bar_binned[0][i];
        double K1cc = cache.K_i_binned[1][i] + cache.K_i_bar_binned[1][i];
        double res = (2 * K1ss + K1cc) / 2 * cache.life_L;
        out.emplace_back(ObservableMapper::to_id(oid), res, this->bins.value()[i]);
    }   
    return out;
}

double LbLllDecay::dG_dq2_avg_bin(size_t bin) {
    double K1ss = cache.K_i_binned[0][bin] + cache.K_i_bar_binned[0][bin];
    double K1cc = cache.K_i_binned[1][bin] + cache.K_i_bar_binned[1][bin];
    return (2 * K1ss + K1cc) / 2;
}

std::vector<ObservableValue> LbLllDecay::A_FB_l(Observables oid) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double K1c = cache.K_i_binned[2][i] + cache.K_i_bar_binned[2][i];
        double res = 0.75 * K1c / dG_dq2_avg_bin(i);
        out.emplace_back(ObservableMapper::to_id(oid), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> LbLllDecay::A_FB_h(Observables oid) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double K2ss = cache.K_i_binned[3][i] + cache.K_i_bar_binned[3][i];
        double K2cc = cache.K_i_binned[4][i] + cache.K_i_bar_binned[4][i];
        double res = 0.25 * (2 * K2ss + K2cc) / dG_dq2_avg_bin(i);
        out.emplace_back(ObservableMapper::to_id(oid), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> LbLllDecay::A_FB_lh(Observables oid) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double K2c = cache.K_i_binned[5][i] + cache.K_i_bar_binned[5][i];
        double res = 0.375 * K2c / dG_dq2_avg_bin(i);
        out.emplace_back(ObservableMapper::to_id(oid), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> LbLllDecay::F_L(Observables oid) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double K1ss = cache.K_i_binned[0][i] + cache.K_i_bar_binned[0][i];
        double K1cc = cache.K_i_binned[1][i] + cache.K_i_bar_binned[1][i];
        double res = 0.5 * (2 * K1ss - K1cc) / dG_dq2_avg_bin(i);
        out.emplace_back(ObservableMapper::to_id(oid), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> LbLllDecay::F_T(Observables oid) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double K1cc = cache.K_i_binned[1][i] + cache.K_i_bar_binned[1][i];
        double res = K1cc / dG_dq2_avg_bin(i);
        out.emplace_back(ObservableMapper::to_id(oid), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> LbLllDecay::compute_observable(Observables obs) {
    switch (obs) {
    case Observables::DBR_DQ2_LAMBDA_B__LAMBDA_E_E:   
        set_cfg_flags(LbLllConfig::Lepton::E);
        return dBR_dq2_binned(obs);
    case Observables::A_FB_L_LAMBDA_B__LAMBDA_E_E:   
        set_cfg_flags(LbLllConfig::Lepton::E);
        return A_FB_l(obs);
    case Observables::A_FB_H_LAMBDA_B__LAMBDA_E_E:
        set_cfg_flags(LbLllConfig::Lepton::E);
        return A_FB_h(obs);
    case Observables::A_FB_LH_LAMBDA_B__LAMBDA_E_E:
        set_cfg_flags(LbLllConfig::Lepton::E);
        return A_FB_lh(obs);
    case Observables::F_L_LAMBDA_B__LAMBDA_E_E:
        set_cfg_flags(LbLllConfig::Lepton::E);
        return F_L(obs);
    case Observables::F_T_LAMBDA_B__LAMBDA_E_E:
        set_cfg_flags(LbLllConfig::Lepton::E);
        return F_T(obs);
    case Observables::DBR_DQ2_LAMBDA_B__LAMBDA_MU_MU:
        set_cfg_flags(LbLllConfig::Lepton::MU);
        return dBR_dq2_binned(obs);
    case Observables::A_FB_L_LAMBDA_B__LAMBDA_MU_MU:
        set_cfg_flags(LbLllConfig::Lepton::MU);
        return A_FB_l(obs);
    case Observables::A_FB_H_LAMBDA_B__LAMBDA_MU_MU:
        set_cfg_flags(LbLllConfig::Lepton::MU);   
        return A_FB_h(obs);
    case Observables::A_FB_LH_LAMBDA_B__LAMBDA_MU_MU:
        set_cfg_flags(LbLllConfig::Lepton::MU);   
        return A_FB_lh(obs);
    case Observables::F_L_LAMBDA_B__LAMBDA_MU_MU:
        set_cfg_flags(LbLllConfig::Lepton::MU);   
        return F_L(obs);
    case Observables::F_T_LAMBDA_B__LAMBDA_MU_MU:
        set_cfg_flags(LbLllConfig::Lepton::MU);   
        return F_T(obs);
    case Observables::DBR_DQ2_LAMBDA_B__LAMBDA_TAU_TAU:
        set_cfg_flags(LbLllConfig::Lepton::TAU);   
        return dBR_dq2_binned(obs);
    case Observables::A_FB_L_LAMBDA_B__LAMBDA_TAU_TAU:
        set_cfg_flags(LbLllConfig::Lepton::TAU);   
        return A_FB_l(obs);
    case Observables::A_FB_H_LAMBDA_B__LAMBDA_TAU_TAU:
        set_cfg_flags(LbLllConfig::Lepton::TAU);   
        return A_FB_h(obs);
    case Observables::A_FB_LH_LAMBDA_B__LAMBDA_TAU_TAU:
        set_cfg_flags(LbLllConfig::Lepton::TAU);   
        return A_FB_lh(obs);
    case Observables::F_L_LAMBDA_B__LAMBDA_TAU_TAU:
        set_cfg_flags(LbLllConfig::Lepton::TAU);   
        return F_L(obs);
    case Observables::F_T_LAMBDA_B__LAMBDA_TAU_TAU:
        set_cfg_flags(LbLllConfig::Lepton::TAU);   
        return F_T(obs);
    default:
        LOG_ERROR("IndexError", "Observable", ObservableMapper::str(obs), "doesn't belong to the decay", DecayMapper::str(this->id));
    }
}

std::vector<ObservableValue> LbLllDecay::compute_observable(ObservableId obs) {
    return compute_observable(ObservableMapper::enum_of(obs).value());
}
