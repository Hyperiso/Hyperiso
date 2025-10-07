#include "M0_Mixing.h"

void M0Mixing::load_params() {
    ObsParameterProxy p;
    cache.G_F = p(ParamId{ParameterType::SM, "SMINPUTS", 2});
    cache.m_W = p(ParamId{ParameterType::SM, "MASS", 24});
    cache.mu_W = w_config.matching_scale;
    cache.x_c = pow(p(ParamId{ParameterType::SM, "MASS", 4}) / cache.m_W, 2);
    cache.x_t = std::pow(ObsQCDProxy()(MassConfig{6, cache.mu_W, MassType::POLE, MassType::POLE}) / cache.m_W, 2);
    cache.lambda_c = std::conj(p(ParamId{ParameterType::SM, "VCKM", {1, 0}})) * p(ParamId{ParameterType::SM, "VCKM", {1, 1}});
    cache.lambda_t = std::conj(p(ParamId{ParameterType::SM, "VCKM", {2, 0}})) * p(ParamId{ParameterType::SM, "VCKM", {2, 1}});

    cache.m_Bd = p(ParamId{ParameterType::FLAVOR, "FMASS", 511});
    cache.m_Bs = p(ParamId{ParameterType::FLAVOR, "FMASS", 531});
    cache.m_K = p(ParamId{ParameterType::FLAVOR, "FMASS", 311});
    cache.m_D = p(ParamId{ParameterType::FLAVOR, "FMASS", 421});
    cache.tau_D = p(ParamId{ParameterType::FLAVOR, "FLIFE", 421});
    cache.G12_s = p(ParamId{ParameterType::DECAY, "M0_Mix", {5, 1}}) * std::exp(I * p(ParamId{ParameterType::DECAY, "M0_Mix", {5, 2}}));
    cache.delta_G_s = p(ParamId{ParameterType::DECAY, "M0_Mix", 6});
    cache.kappa_e = p(ParamId{ParameterType::DECAY, "M0_Mix", 7});
    cache.eta_cc = p(ParamId{ParameterType::DECAY, "M0_Mix", 11});
    cache.eta_ct = p(ParamId{ParameterType::DECAY, "M0_Mix", 12});
    cache.alpha_s_mu_W = ObsQCDProxy()(AlphasConfig{cache.mu_W, MassType::POLE, MassType::POLE});
    cache.alpha_s_mu_b = ObsQCDProxy()(AlphasConfig{p(ParamId{ParameterType::DECAY, "M0_Mix", 8}), MassType::POLE, MassType::POLE});
    cache.alpha_s_mu_c = ObsQCDProxy()(AlphasConfig{p(ParamId{ParameterType::DECAY, "M0_Mix", 9}), MassType::POLE, MassType::POLE});

    std::array<double, 5> B {
        p(ParamId{ParameterType::FLAVOR, "FBAG", LhaID(511, 1)}),
        p(ParamId{ParameterType::FLAVOR, "FBAG", LhaID(511, 2)}),
        p(ParamId{ParameterType::FLAVOR, "FBAG", LhaID(511, 3)}),
        p(ParamId{ParameterType::FLAVOR, "FBAG", LhaID(511, 4)}),
        p(ParamId{ParameterType::FLAVOR, "FBAG", LhaID(511, 5)}),
    };
    double mf2 = pow(cache.m_Bd * p(ParamId{ParameterType::FLAVOR, "FCONST", {511, 1}}), 2);
    populate_Q_from_bag(cache.Q_Bd, B, p(ParamId{ParameterType::DECAY, "M0_Mix", 1}), mf2, true);

    B = {
        p(ParamId{ParameterType::FLAVOR, "FBAG", LhaID(311, 1)}),
        p(ParamId{ParameterType::FLAVOR, "FBAG", LhaID(311, 2)}),
        p(ParamId{ParameterType::FLAVOR, "FBAG", LhaID(311, 3)}),
        p(ParamId{ParameterType::FLAVOR, "FBAG", LhaID(311, 4)}),
        p(ParamId{ParameterType::FLAVOR, "FBAG", LhaID(311, 5)}),
    };
    mf2 = pow(cache.m_K * p(ParamId{ParameterType::FLAVOR, "FCONST", {311, 1}}), 2);
    populate_Q_from_bag(cache.Q_K, B, p(ParamId{ParameterType::DECAY, "M0_Mix", 2}), mf2, true);

    cache.Q_D = {
        p(ParamId{ParameterType::DECAY, "M0_Mix", LhaID(4, 1)}),
        p(ParamId{ParameterType::DECAY, "M0_Mix", LhaID(4, 2)}),
        p(ParamId{ParameterType::DECAY, "M0_Mix", LhaID(4, 3)}),
        p(ParamId{ParameterType::DECAY, "M0_Mix", LhaID(4, 4)}),
        p(ParamId{ParameterType::DECAY, "M0_Mix", LhaID(4, 5)}),
    };

    populate_C(cache.C_Bd, p(ParamId{ParameterType::DECAY, "M0_Mix", 8}), 0);
    populate_C(cache.C_Bs, p(ParamId{ParameterType::DECAY, "M0_Mix", 8}), 8);
    populate_C(cache.C_K, p(ParamId{ParameterType::DECAY, "M0_Mix", 9}), 16);
    populate_C(cache.C_D, p(ParamId{ParameterType::DECAY, "M0_Mix", 10}), 24);

    cache.C1_Bd_SM = std::pow(cache.G_F * cache.m_W * std::conj(p(ParamId{ParameterType::SM, "VCKM", {2, 0}})) * p(ParamId{ParameterType::SM, "VCKM", {2, 2}}), 2) * S0(cache.x_t) / (4 * PI2);
    cache.C1_Bs_SM = std::pow(cache.G_F * cache.m_W * std::conj(p(ParamId{ParameterType::SM, "VCKM", {2, 1}})) * p(ParamId{ParameterType::SM, "VCKM", {2, 2}}), 2) * S0(cache.x_t) / (4 * PI2);
}

double M0Mixing::S_18(double x) {
    return -(64.-68.*x-17.*x*x+11.*x*x*x)/(4.*pow(1.-x,2.)) + (32.-68.*x+32.*x*x-28.*x*x*x+3.*pow(x,4.))/(2.*pow(1.-x,3.))*log(x) + x*x*(4.-7.*x+7.*x*x-2.*x*x*x)/(2.*pow(1.-x,4.))*pow(log(x),2.) + 2.*x*(4. - 7.*x - 7.*x*x + x*x*x) / pow(1.-x,3.) * Li2(1.-x) + 16./x*(PI2/6. - Li2(1.-x));;
}

double M0Mixing::S_11(double x) {
    return -x*(4.-39.*x+168.*x*x+11.*x*x*x)/(4.*pow(1.-x,3.)) - 3.*x*(4.-24.*x+36.*x*x+7.*x*x*x+pow(x,4.))/(2.*pow(1.-x,4.))*log(x) + 3.*x*x*x*(13.+4.*x+x*x)/(2.*pow(1.-x,4.))*pow(log(x),2.) - 3*x*x*x*(5.+x)/pow(1.-x,3.)*Li2(1.-x);
}

double M0Mixing::S_1t(double x) {
    return (S_18(x) + 4. * S_11(x)) / 3.;
}

double M0Mixing::dS0_dx(double x) {
    return 3.*x*x/(2.*pow(x-1.,3.)) + x*(-11.+2.*x)/(4.*pow(x-1.,2.)) + (4.-11.*x+x*x)/(4.*pow(x-1.,2.)) - x*(4.-11.*x+x*x)/(2.*pow(x-1.,3.)) + 9.*x*x*log(x)/(2.*pow(x-1.,3.)) - 9.*x*x*x*log(x)/(2.*pow(x-1.,4.));
}

double M0Mixing::S0_ct(double x_c, double x_t) {
    return x_c*(std::log(x_t/x_c) - 3.*x_t/(4.*(1-x_t)) -3*x_t*x_t*std::log(x_t)/(4.*std::pow(1.-x_t,2.)));
}

double M0Mixing::F_S1(double x, double m_W, double mu_W, double mu_t) {
    return 8 * x * dS0_dx(x) * std::log(std::pow(mu_t / m_W, 2)) + 2 * S0(x) * std::log(std::pow(mu_W / m_W, 2));
}

double M0Mixing::Q_i(int i, double B_i, double r_chi, double mf2, bool is_B) {
    double Q_i = N_i.at(i) * mf2 * B_i;
    if (i > 1) {
        Q_i *= is_B ? r_chi + d_i.at(i) : r_chi;
    }
    return Q_i;
}

void M0Mixing::populate_Q_from_bag(std::array<double, 5>& Q, const std::array<double, 5>& B, double r_chi, double mf2, bool is_B) {
    for (size_t i = 0; i < 5; i++) {
        Q[i] = Q_i(i, B[i], r_chi, mf2, is_B);
    }
}

void M0Mixing::populate_C(std::array<complex_t, 8>& C, double hadronic_scale, size_t offset) {
    ObsParameterMutator().set(ParamId{ParameterType::WILSON, "B_SCALE", 1}, hadronic_scale);
    auto ids = WCoefMapper::get_group(WGroup::MESON_MIXING);
    for (size_t i = 0; i < 8; i++) {
        C[i] = this->w_proxy->getFR(WGroup::MESON_MIXING, ids[i + offset], this->w_config.order, ContributionType::BSM);
    }
}

complex_t M0Mixing::M_12_NP(const std::array<complex_t, 8>& C, const std::array<double, 5>& Q, double m_M) {
    complex_t M_12 {0.0};
    for (size_t i = 0; i < 8; i++) {
        M_12 += C[i] * Q[i % 5];
    }
    
    return M_12 / (2 * m_M);
}

complex_t M0Mixing::M_12_B_SM(int gen) {
    double eta_2B = std::pow(cache.alpha_s_mu_W, 6./23) * (1 + cache.alpha_s_mu_W / (4 * PI) * ((S_1t(cache.x_t) + F_S1(cache.x_t, cache.m_W, cache.mu_W, cache.mu_W)) / S0(cache.x_t) + B_t - J_5));
    double c_RGI_B = std::pow(cache.alpha_s_mu_b, -6./23) * (1 + cache.alpha_s_mu_b * J_5 / (4 * PI));
    return gen == 1 ? cache.C1_Bd_SM * eta_2B * c_RGI_B * cache.Q_Bd[0] / (2 * cache.m_Bd) : cache.C1_Bs_SM * eta_2B * c_RGI_B * cache.Q_Bs[0] / (2 * cache.m_Bs);
}

complex_t M0Mixing::M_12_K_SM() {
    double eta_tt = std::pow(cache.alpha_s_mu_c, 2./9.)*std::pow(cache.alpha_s_mu_b/cache.alpha_s_mu_c, 6./25.)*std::pow(cache.alpha_s_mu_W/cache.alpha_s_mu_b, 6./23.)*( 1. + cache.alpha_s_mu_c/4./PI*(J_4-J_3) + cache.alpha_s_mu_b/4./PI*(J_5-J_4) +cache.alpha_s_mu_W/4./PI*( (S_1t(cache.x_t) + F_S1(cache.x_t, cache.m_W, cache.mu_W, cache.mu_W))/S0(cache.x_t) + B_t - J_5));
    complex_t C_1_sd = std::pow(cache.lambda_c, 2) * cache.eta_cc * cache.x_c + std::pow(cache.lambda_t, 2) * eta_tt * S0(cache.x_t) + 2. * cache.lambda_t * cache.lambda_c * cache.eta_ct * S0_ct(cache.x_c, cache.x_t);
    double c_RGI_K = std::pow(cache.alpha_s_mu_c, -2./9) * (1 + cache.alpha_s_mu_c * J_3 / (4 * PI));
    return std::pow(cache.G_F * cache.m_W / (2 * PI), 2) / (2 * cache.m_K) * c_RGI_K * C_1_sd * cache.Q_K[0];
}

// B mixing observables

double M0Mixing::delta_M_B(int gen) {
    complex_t M_12 = gen == 1 ? M_12_NP(cache.C_Bd, cache.Q_Bd, cache.m_Bd) : M_12_NP(cache.C_Bs, cache.Q_Bs, cache.m_Bs);
    M_12 += M_12_B_SM(gen);
    return 2 * std::abs(M_12) * GEV_TO_INV_PS;
}

double M0Mixing::phi_q(int gen) {
    complex_t M_12 = gen == 1 ? M_12_NP(cache.C_Bd, cache.Q_Bd, cache.m_Bd) : M_12_NP(cache.C_Bs, cache.Q_Bs, cache.m_Bs);
    M_12 += M_12_B_SM(gen);
    return std::arg(M_12);
}

// B_s specific observables

double M0Mixing::a_fs() {
    complex_t M_12 = M_12_NP(cache.C_Bs, cache.Q_Bs, cache.m_Bs) + M_12_B_SM(2);
    return std::tan(std::arg(-M_12 / cache.G12_s)) * cache.delta_G_s / (2 * std::abs(M_12));
}

// K mixing observables

double M0Mixing::epsilon_K() {
    complex_t M_12 = M_12_NP(cache.C_K, cache.Q_K, cache.m_K) + M_12_K_SM();
    return cache.kappa_e / (2 * RT2) * std::abs(std::imag(M_12) / std::real(M_12)); 
}

double M0Mixing::delta_M_K() {
    complex_t M_12 = M_12_NP(cache.C_K, cache.Q_K, cache.m_K) + M_12_K_SM();
    return 2. * std::real(M_12) * GEV_TO_INV_PS;
}

// D mixing observables

double M0Mixing::x_D() {
    complex_t M_12 = M_12_NP(cache.C_D, cache.Q_D, cache.m_D);
    return 2 * cache.tau_D * std::abs(M_12) * GEV_TO_INV_S;
}

std::vector<ObservableValue> M0Mixing::compute_observable(Observables obs) {
    double value;
    switch (obs) {
    case Observables::PHI_D:   
        value = phi_q(1);
        break;
    case Observables::DELTA_M_BD:   
        value = delta_M_B(1);
        break;
    case Observables::PHI_S:   
        value = phi_q(2);
        break;
    case Observables::DELTA_M_BS:   
        value = delta_M_B(2);
        break;
    case Observables::A_FS:   
        value = a_fs();
        break;
    case Observables::DELTA_M_K:   
        value = delta_M_K();
        break;
    case Observables::ABS_EPSILON_K:   
        value = epsilon_K();
        break;
    case Observables::X_D:   
        value = x_D();
        break;
    default:
        LOG_ERROR("IndexError", "Observable", ObservableMapper::str(obs), "doesn't belong to the decay", DecayMapper::str(this->id));
    }

    return {ObservableValue(ObservableMapper::to_id(obs), value)};
}

std::vector<ObservableValue> M0Mixing::compute_observable(ObservableId obs) {
    return compute_observable(ObservableMapper::enum_of(obs).value());
}