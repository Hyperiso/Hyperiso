#include "M0_Mixing.h"

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

double M0Mixing::alpha(double mu) {
    return ObsQCDProxy()(AlphasConfig{mu, MassType::POLE, MassType::POLE});
}

double M0Mixing::x_t(double mu, double m_W) {
    return std::pow(ObsQCDProxy()(MassConfig{6, mu, MassType::POLE, MassType::POLE}) / m_W, 2);
}

double M0Mixing::c_RGI_B(double alpha_s_mu_b) {
    return std::pow(alpha_s_mu_b, -6./23) * (1 + alpha_s_mu_b * J_5 / (4 * PI));
}

double M0Mixing::eta_2B(double alpha_s_mu_W, double x_t, double m_W) {
    double mu_W = this->w_config.matching_scale;
    return std::pow(alpha_s_mu_W, 6./23) * (1 + alpha_s_mu_W / (4 * PI) * ((S_1t(x_t) + F_S1(x_t, m_W, mu_W, mu_W)) / S0(x_t) + B_t - J_5));
}

complex_t M0Mixing::lambda(complex_t V_qd, complex_t V_qs) {
    return std::conj(V_qd) * V_qs;
}

double M0Mixing::c_RGI_K(double alpha_s_mu_c) {
    return std::pow(alpha_s_mu_c, -2./9) * (1 + alpha_s_mu_c * J_3 / (4 * PI));
}

double M0Mixing::eta_tt(double alpha_s_mu_c, double alpha_s_mu_b, double alpha_s_mu_W, double x_t, double m_W) {
    double mu_W = this->w_config.matching_scale;
    return std::pow(alpha_s_mu_c,2./9.)*std::pow(alpha_s_mu_b/alpha_s_mu_c, 6./25.)*std::pow(alpha_s_mu_W/alpha_s_mu_b, 6./23.)*( 1. + alpha_s_mu_c/4./PI*(J_4-J_3) + alpha_s_mu_b/4./PI*(J_5-J_4) +alpha_s_mu_W/4./PI*( (S_1t(x_t) + F_S1(x_t, m_W, mu_W, mu_W))/S0(x_t) + B_t - J_5));
}

complex_t M0Mixing::C_sd_SM(double x_t, double x_c, complex_t lambda_c, complex_t lambda_t, double eta_cc, double eta_tt, double eta_ct) {
    return std::pow(lambda_c, 2) * eta_cc * x_c + std::pow(lambda_t, 2) * eta_tt * S0(x_t) + 2. * lambda_t * lambda_c * eta_ct * S0_ct(x_c, x_t);
}

double M0Mixing::Q_i(int i, double B_i, double r_chi, double mf2, bool is_B) {
    double Q_i = N_i.at(i) * mf2 * B_i;
    if (i > 1) {
        Q_i *= is_B ? r_chi + d_i.at(i) : r_chi;
    }
    return Q_i;
}

void M0Mixing::populate_B(double B_1, double B_2, double B_3, double B_4, double B_5) {
    std::array<double, 5> B_temp {B_1, B_2, B_3, B_4, B_5};
    this->B = B_temp;
}

void M0Mixing::populate_Q_from_bag(double r_chi, double mf2, bool is_B) {
    for (size_t i = 0; i < 5; i++) {
        this->Q[i] = Q_i(i, this->B[i], r_chi, mf2, is_B);
    }
}

void M0Mixing::populate_Q_from_mels(double Q_1, double Q_2, double Q_3, double Q_4, double Q_5) {
    std::array<double, 5> Q_temp {Q_1, Q_2, Q_3, Q_4, Q_5};
    this->Q = Q_temp;
}

void M0Mixing::populate_C(double hadronic_scale, size_t offset) {
    ObsParameterMutator().set(ParamId{ParameterType::WILSON, "B_SCALE", 1}, hadronic_scale);
    auto ids = WCoefMapper::get_group(WGroup::MESON_MIXING);
    for (size_t i = 0; i < 8; i++) {
        this->C[i] = this->w_proxy->getFR(WGroup::MESON_MIXING, ids[i + offset], this->w_config.order, ContributionType::BSM);
    }
}

complex_t M0Mixing::M_12_NP(double m_M) {
    complex_t M_12 {0.0};
    for (size_t i = 0; i < 8; i++) {
        M_12 += this->C[i] * this->Q[i % 5];
    }
    
    return M_12 / (2 * m_M);
}

complex_t M0Mixing::M_12_B_SM(double m_B, double eta_2B, double c_RGI_B, int gen) {
    complex_t C_1 = this->w_proxy->getFM(WGroup::MESON_MIXING, gen == 1 ? WCoef::C_BD_1 : WCoef::C_BS_1, this->w_config.order, ContributionType::SM);
    return C_1 * eta_2B * c_RGI_B * Q[0] / (2 * m_B);
}

complex_t M0Mixing::M_12_K_SM(double m_M, double G_F, double m_W, double c_RGI_K, complex_t C_1_sd) {
    return std::pow(G_F * m_W / (2 * PI), 2) / (2 * m_M) * c_RGI_K * C_1_sd * Q[0];
}

// B mixing observables

double M0Mixing::delta_M_B(complex_t M_12) {
    return 2 * std::abs(M_12) * GEV_TO_INV_PS;
}

double M0Mixing::phi_q(complex_t M_12) {
    return std::arg(M_12);
}

// B_s specific observables

double M0Mixing::a_fs(complex_t M_12, complex_t G_12, double delta_M, double delta_G) {
    return std::tan(std::arg(-M_12 / G_12)) * delta_G / delta_M;
}

// K mixing observables

double M0Mixing::epsilon_K(complex_t M_12, double kappa_e) {
    return kappa_e / (2 * RT2) * std::abs(std::imag(M_12) / std::real(M_12)); 
}

double M0Mixing::delta_M_K(complex_t M_12) {
    return 2. * std::real(M_12) * GEV_TO_INV_PS;
}

// D mixing observables

double M0Mixing::x_D(complex_t M_12, double tau_D0) {
    return 2 * tau_D0 * std::abs(M_12) * GEV_TO_INV_S;
}

void M0Mixing::build_op_tree() {

    // SM Parameters
    auto G_F = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "SMINPUTS", 2));
    auto m_c = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 4));
    auto m_W = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 24));
    auto V_td = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "VCKM", {2, 0}));
    auto V_ts = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "VCKM", {2, 1}));
    auto V_cd = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "VCKM", {1, 0}));
    auto V_cs = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "VCKM", {1, 1}));

    // Wilson node
    auto wilson = this->get_wilson_node();

    // Flavor Parameters
    auto m_Bs = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FMASS", 531));
    auto m_Bd = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FMASS", 511));
    auto m_K = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FMASS", 311));
    auto m_D = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FMASS", 421));

    auto tau_D = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FLIFE", 421));

    auto f_Bs = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FCONST", LhaID(531,1)));
    auto f_Bd = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FCONST", LhaID(511,1)));
    auto f_K = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FCONST", LhaID(321,1)));

    auto B1_Bs = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FBAG", LhaID(531,1)));
    auto B2_Bs = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FBAG", LhaID(531,2)));
    auto B3_Bs = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FBAG", LhaID(531,3)));
    auto B4_Bs = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FBAG", LhaID(531,4)));
    auto B5_Bs = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FBAG", LhaID(531,5)));

    auto B1_Bd = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FBAG", LhaID(511,1)));
    auto B2_Bd = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FBAG", LhaID(511,2)));
    auto B3_Bd = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FBAG", LhaID(511,3)));
    auto B4_Bd = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FBAG", LhaID(511,4)));
    auto B5_Bd = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FBAG", LhaID(511,5)));

    auto B1_K = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FBAG", LhaID(311,1)));
    auto B2_K = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FBAG", LhaID(311,2)));
    auto B3_K = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FBAG", LhaID(311,3)));
    auto B4_K = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FBAG", LhaID(311,4)));
    auto B5_K = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FBAG", LhaID(311,5)));

    // Misc experimental input
    auto r_chi_Bd = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "M0_Mix", 1));
    auto r_chi_Bs = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "M0_Mix", 2));
    auto r_chi_K = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "M0_Mix", 3));
    auto Q1_D = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "M0_Mix", {4, 1}));
    auto Q2_D = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "M0_Mix", {4, 2}));
    auto Q3_D = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "M0_Mix", {4, 3}));
    auto Q4_D = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "M0_Mix", {4, 4}));
    auto Q5_D = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "M0_Mix", {4, 5}));
    auto abs_G12_s = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "M0_Mix", {5, 1}));
    auto arg_G12_s = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "M0_Mix", {5, 2}));
    auto delta_G_s = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "M0_Mix", 6));
    auto kappa_e = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "M0_Mix", 7));
    auto eta_cc = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "M0_Mix", 11));
    auto eta_ct = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "M0_Mix", 12));

    auto mu_B = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "M0_Mix", 8));
    auto mu_K = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "M0_Mix", 9));
    auto mu_D = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "M0_Mix", 10));
    
    // Operator nodes for B_d observables
    auto alpha_s_mu_W = std::make_shared<OperatorNode>("alpha_s(mu_b)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return alpha(this->w_config.matching_scale); });
    auto alpha_s_mu_b = std::make_shared<OperatorNode>("alpha_s(mu_W)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return alpha(values[0]); });
    alpha_s_mu_b->addChild(mu_B);
    auto xt = std::make_shared<OperatorNode>("x_t", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return x_t(this->w_config.matching_scale, values[0]); });
    xt->addChild(m_W);
    auto eta = std::make_shared<OperatorNode>("eta_2B", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return eta_2B(values[0], values[1], values[2]); });
    eta->addChildren({alpha_s_mu_W, xt, m_W});
    auto crgi = std::make_shared<OperatorNode>("c_RGI_B", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return c_RGI_B(values[0]); });
    crgi->addChildren({alpha_s_mu_b});
    auto wilson_cache_Bd = std::make_shared<OperatorNode>("Wilson cache", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { populate_C(values[0]); return 0; });
    wilson_cache_Bd->addChildren({mu_B, wilson});
    auto mf2_Bd = std::make_shared<OperatorNode>("(f_Bd m_Bd)^2", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return pow(values[0] * values[1], 2); });
    mf2_Bd->addChildren({f_Bd, m_Bd});
    auto Q_Bd = std::make_shared<OperatorNode>("Q_Bd", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { 
        populate_B(values[0], values[1], values[2], values[3], values[4]);
        populate_Q_from_bag(values[5], values[6], true);
        LOG_INFO("Q_Bd_1 =", Q[0]);
        return 0;
    });
    Q_Bd->addChildren({B1_Bd, B2_Bd, B3_Bd, B4_Bd, B5_Bd, r_chi_Bd, mf2_Bd});
    auto M12_Bd_SM = std::make_shared<OperatorNode>("M12_Bd_SM", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return M_12_B_SM(values[0], values[1], values[2], 1); });
    M12_Bd_SM->addChildren({m_Bd, eta, crgi, wilson, Q_Bd});
    auto M12_Bd_NP = std::make_shared<OperatorNode>("M12_Bd_NP", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return M_12_NP(values[0]); });
    M12_Bd_NP->addChildren({m_Bd, Q_Bd, wilson_cache_Bd});
    auto M12_Bd = std::make_shared<OperatorNode>("M12_Bd", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] + values[1]; });
    M12_Bd->addChildren({M12_Bd_SM, M12_Bd_NP});
    auto phi_d = std::make_shared<OperatorNode>("phi_d", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return phi_q(values[0]); });
    phi_d->addChild(M12_Bd);
    roots.emplace(ObservableMapper::to_id(Observables::PHI_D), phi_d);
    auto delta_M_Bd = std::make_shared<OperatorNode>("Delta_M_Bd", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return delta_M_B(values[0]); });
    delta_M_Bd->addChild(M12_Bd);
    roots.emplace(ObservableMapper::to_id(Observables::DELTA_M_BD), delta_M_Bd);

    // Operator nodes for B_s observables
    auto wilson_cache_Bs = std::make_shared<OperatorNode>("Wilson cache", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { populate_C(values[0], 8); return 0; });
    wilson_cache_Bs->addChildren({mu_B, wilson});
    auto mf2_Bs = std::make_shared<OperatorNode>("(f_Bs m_Bs)^2", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return pow(values[0] * values[1], 2); });
    mf2_Bs->addChildren({f_Bs, m_Bs});
    auto Q_Bs = std::make_shared<OperatorNode>("Q_Bs", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { 
        populate_B(values[0], values[1], values[2], values[3], values[4]);
        populate_Q_from_bag(values[5], values[6], true);
        LOG_INFO("Q_Bs_1 =", Q[0]);
        return 0;
    });
    Q_Bs->addChildren({B1_Bs, B2_Bs, B3_Bs, B4_Bs, B5_Bs, r_chi_Bs, mf2_Bs});
    auto M12_Bs_SM = std::make_shared<OperatorNode>("M12_Bs_SM", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return M_12_B_SM(values[0], values[1], values[2], 2); });
    M12_Bs_SM->addChildren({m_Bs, eta, crgi, wilson, Q_Bs});
    auto M12_Bs_NP = std::make_shared<OperatorNode>("M12_Bs", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return M_12_NP(values[0]); });
    M12_Bs_NP->addChildren({m_Bs, Q_Bs, wilson_cache_Bs});
    auto M12_Bs = std::make_shared<OperatorNode>("M12_Bs", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] + values[1]; });
    M12_Bs->addChildren({M12_Bs_SM, M12_Bs_NP});
    auto phi_s = std::make_shared<OperatorNode>("phi_s", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return phi_q(values[0]); });
    phi_s->addChild(M12_Bs);
    roots.emplace(ObservableMapper::to_id(Observables::PHI_S), phi_s);
    auto delta_M_Bs = std::make_shared<OperatorNode>("Delta_M_Bs", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return delta_M_B(values[0]); });
    delta_M_Bs->addChild(M12_Bs);
    roots.emplace(ObservableMapper::to_id(Observables::DELTA_M_BS), delta_M_Bs);
    auto Gamma_12_Bs = std::make_shared<OperatorNode>("Gamma_12_Bs", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] * std::exp(I * values[1]); });
    Gamma_12_Bs->addChildren({abs_G12_s, arg_G12_s});
    auto afs = std::make_shared<OperatorNode>("afs", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return a_fs(values[0], values[1], values[2], values[3]); });
    afs->addChildren({M12_Bs, Gamma_12_Bs, delta_M_Bs, delta_G_s});
    roots.emplace(ObservableMapper::to_id(Observables::A_FS), afs);

    // Operator nodes for K observables
    auto wilson_cache_K = std::make_shared<OperatorNode>("Wilson cache", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { populate_C(values[0], 16); return 0; });
    wilson_cache_K->addChildren({mu_K, wilson});
    auto mf2_K = std::make_shared<OperatorNode>("(f_K m_K)^2", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return pow(values[0] * values[1], 2); });
    mf2_K->addChildren({f_K, m_K});
    auto Q_K = std::make_shared<OperatorNode>("Q_K", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { 
        populate_B(values[0], values[1], values[2], values[3], values[4]);
        populate_Q_from_bag(values[5], values[6], false);
        return 0;
    });
    Q_K->addChildren({B1_K, B2_K, B3_K, B4_K, B5_K, r_chi_K, mf2_K});
    auto alpha_s_mu_c = std::make_shared<OperatorNode>("alpha_s(mu_K)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return alpha(values[0]); });
    alpha_s_mu_c->addChild(mu_K);
    auto x_c = std::make_shared<OperatorNode>("x_c", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return std::pow(values[0] / values[1], 2); });
    x_c->addChildren({m_c, m_W});
    auto n_lambda_c = std::make_shared<OperatorNode>("lambda_c", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return lambda(values[0], values[1]); });
    n_lambda_c->addChildren({V_cd, V_cs});
    auto n_lambda_t = std::make_shared<OperatorNode>("lambda_t", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return lambda(values[0], values[1]); });
    n_lambda_t->addChildren({V_td, V_ts});
    auto n_eta_tt = std::make_shared<OperatorNode>("eta_tt", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return eta_tt(values[0], values[1], values[2], values[3], values[4]); });
    n_eta_tt->addChildren({alpha_s_mu_c, alpha_s_mu_b, alpha_s_mu_W, xt, m_W});
    auto n_C_sd_SM = std::make_shared<OperatorNode>("C_sd_SM", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return C_sd_SM(values[0], values[1], values[2], values[3], values[4], values[5], values[6]); });
    n_C_sd_SM->addChildren({xt, x_c, n_lambda_c, n_lambda_t, eta_cc, n_eta_tt, eta_ct});
    auto n_c_RGI_K = std::make_shared<OperatorNode>("c_RGI_K", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return c_RGI_K(values[0]); });
    n_c_RGI_K->addChildren({alpha_s_mu_c});
    auto M12_K_SM = std::make_shared<OperatorNode>("M12_K_SM", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return M_12_K_SM(values[0], values[1], values[2], values[3], values[4]); });
    M12_K_SM->addChildren({m_K, G_F, m_W, n_c_RGI_K, n_C_sd_SM, Q_K});
    auto M12_K_NP = std::make_shared<OperatorNode>("M12_K_NP", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return M_12_NP(values[0]); });
    M12_K_NP->addChildren({m_K, Q_K, wilson_cache_K});
    auto M12_K = std::make_shared<OperatorNode>("M12_K", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] + values[1]; });
    M12_K->addChildren({M12_K_SM, M12_K_NP});
    auto abs_epsilon_K = std::make_shared<OperatorNode>("|epsilon_K|", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return std::abs(epsilon_K(values[0], values[1])); });
    abs_epsilon_K->addChildren({M12_K, kappa_e});
    roots.emplace(ObservableMapper::to_id(Observables::ABS_EPSILON_K), abs_epsilon_K);
    auto delta_M_K0 = std::make_shared<OperatorNode>("Delta_M_K", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return delta_M_K(values[0]); });
    delta_M_K0->addChild(M12_K);
    roots.emplace(ObservableMapper::to_id(Observables::DELTA_M_K), delta_M_K0);

    // Operator nodes for D observables
    auto wilson_cache_D = std::make_shared<OperatorNode>("Wilson cache", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { populate_C(values[0], 24); return 0; });
    wilson_cache_D->addChildren({mu_D, wilson});
    auto Q_D = std::make_shared<OperatorNode>("Q_D", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { 
        populate_Q_from_mels(values[0], values[1], values[2], values[3], values[4]);
        return 0;
    });
    Q_D->addChildren({Q1_D, Q2_D, Q3_D, Q4_D, Q5_D});
    auto M12_D = std::make_shared<OperatorNode>("M12_D", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return M_12_NP(values[0]); });
    M12_D->addChildren({m_D, mu_D, Q_D, wilson_cache_D});
    auto n_x_D = std::make_shared<OperatorNode>("x_D", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return x_D(values[0], values[1]); });
    n_x_D->addChildren({M12_D, tau_D});
    roots.emplace(ObservableMapper::to_id(Observables::X_D), n_x_D);

}
