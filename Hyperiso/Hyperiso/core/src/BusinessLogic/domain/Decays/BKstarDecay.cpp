#include "BKstarDecay.h"
#include "BKsllDecay.h"

void BKstarDecay::load_params() {
    ObsParameterProxy p;
    auto run = [this] (double value_1gev, double eta, double gamma) { return value_1gev * pow(eta, gamma / cache.beta_0); };
    auto gamma_perp = [this] (int n) { return 4. * cache.C_F * (psi(n + 1) + GAMMA - 1. + 1. / (n + 1)); };
    auto gamma_par = [this] (int n) { return 4. * cache.C_F * (psi(n + 2) + GAMMA - .75 - 1. /(2. * (n + 1) * (n + 2))); };
    auto get_C = [this] (WCoef id) { return w_proxy->getFR(WGroup::B, id, w_config.order); };
    auto get_Cp = [this] (WCoef id) { return w_proxy->getFR(WGroup::BPrime, id, QCDOrder::LO); };

    cache.mu_b = w_config.hadronic_scale;
    cache.lambda_hat_u = std::conj(p(ParamId{ParameterType::SM, "VCKM", {0, 1}})) * p(ParamId{ParameterType::SM, "VCKM", {0, 2}}) 
                            / (std::conj(p(ParamId{ParameterType::SM, "VCKM", {1, 1}})) * p(ParamId{ParameterType::SM, "VCKM", {1, 2}}));
    cache.m_b_1S = p(ParamId{ParameterType::SM, "QCD", {5, 3}});
    cache.m_b_mu_b = ObsQCDProxy()(MassConfig(5, cache.mu_b, MassType::MSBAR, MassType::POLE));
    cache.s_c = std::pow(ObsQCDProxy()(MassConfig(4, cache.mu_b, MassType::MSBAR, MassType::POLE)) / cache.m_b_mu_b, 2);
    cache.f_Ks_par = p(ParamId{ParameterType::FLAVOR, "FCONST", {323, 1}});
    cache.f_B = p(ParamId{ParameterType::FLAVOR, "FCONST", {521, 1}});
    cache.m_B = p(ParamId{ParameterType::FLAVOR, "FMASS", 521});
    cache.m_Ks = p(ParamId{ParameterType::FLAVOR, "FMASS", 323});
    cache.zeta_3_A = p(ParamId{ParameterType::DECAY, "B_Ks", 5});
    cache.zeta_3_V = p(ParamId{ParameterType::DECAY, "B_Ks", 6});
    cache.omega_10_A = p(ParamId{ParameterType::DECAY, "B_Ks", 7});
    cache.delta_t_p = p(ParamId{ParameterType::DECAY, "B_Ks", 8});
    cache.delta_t_m = p(ParamId{ParameterType::DECAY, "B_Ks", 9});
    cache.T1_B_Ks = p(ParamId{ParameterType::DECAY, "B_Ks", 11});
    cache.Lambda_h = p(ParamId{ParameterType::DECAY, "B_Ks", 12});
    cache.mu_0 = p(ParamId{ParameterType::DECAY, "B_Ks", 13});
    if (fpeq(cache.mu_0, -1.)) cache.mu_0 = cache.mu_b;
    double mu_h = std::sqrt(cache.mu_b * cache.Lambda_h);
    cache.C_F = ObsQCDProxy().get_constants()->C_F;
    cache.Nc = ObsQCDProxy().get_constants()->Nc;
    cache.beta_0 = ObsQCDProxy().get_constants()->beta[5][0];
    cache.n_f = 5.0; // TODO : Link with get_nf vs. hard-coded ?
    cache.alpha_s_mu_b = ObsQCDProxy()(AlphasConfig(cache.mu_b, MassType::POLE, MassType::POLE));
    cache.alpha_s_mu_h = ObsQCDProxy()(AlphasConfig(mu_h, MassType::POLE, MassType::POLE));
    double eta = cache.alpha_s_mu_b / ObsQCDProxy()(AlphasConfig(1.0, MassType::POLE, MassType::POLE));
    cache.f_Ks_perp = run(p(ParamId{ParameterType::FLAVOR, "FCONST", {323, 2}}), eta, cache.C_F);
    cache.a_1_perp = run(p(ParamId{ParameterType::DECAY, "B_Ks", 1}), eta, gamma_perp(1));
    cache.a_2_perp = run(p(ParamId{ParameterType::DECAY, "B_Ks", 2}), eta, gamma_perp(2));
    cache.a_1_par = run(p(ParamId{ParameterType::DECAY, "B_Ks", 3}), eta, gamma_par(1));
    cache.a_2_par = run(p(ParamId{ParameterType::DECAY, "B_Ks", 4}), eta, gamma_par(2));
    cache.lambda_B = p(ParamId{ParameterType::DECAY, "B_Ks", 10}) / (1. - cache.alpha_s_mu_h * log(pow(mu_h, 2)) * 1.8 / (3. * PI));
    
    cache.C_b.emplace(WCoef::C1, get_C(WCoef::C1) /* + get_Cp(WCoef::CP1) */ );
    cache.C_b.emplace(WCoef::C2, get_C(WCoef::C2) /* + get_Cp(WCoef::CP2) */ );
    cache.C_b.emplace(WCoef::C3, get_C(WCoef::C3) /* + get_Cp(WCoef::CP3) */ );
    cache.C_b.emplace(WCoef::C4, get_C(WCoef::C4) /* + get_Cp(WCoef::CP4) */ );
    cache.C_b.emplace(WCoef::C5, get_C(WCoef::C5) /* + get_Cp(WCoef::CP5) */ );
    cache.C_b.emplace(WCoef::C6, get_C(WCoef::C6) /* + get_Cp(WCoef::CP6) */ );
    cache.C_b.emplace(WCoef::C7, get_C(WCoef::C7) + get_Cp(WCoef::CP7));
    cache.C_b.emplace(WCoef::C8, get_C(WCoef::C8) + get_Cp(WCoef::CP8));

    ObsParameterMutator().set(ParamId{ParameterType::WILSON, "B_SCALE", 1}, mu_h);
    cache.C2_h = get_C(WCoef::C2) /* + get_Cp(WCoef::CP2) */;
    cache.C8_h = get_C(WCoef::C8) + get_Cp(WCoef::CP8);
}

complex_t BKstarDecay::h(double u) {
    complex_t rt = sqrt((u - 4 * cache.s_c + I * EPSILON) / u);
    return (4 * cache.s_c * (CLi2(2. / (1. - rt)) + CLi2(2. / (1. + rt))) / u - 2.) / u;
}

double BKstarDecay::phi_perp(double u) {
    double ubar = 1 - u;
    double xi = u - ubar;
    return 6 * u * ubar * (1 + 3 * cache.a_1_perp * xi + 1.5 * cache.a_2_perp * (5 * xi * xi - 1));
}

double BKstarDecay::gv_dga_4(double u) {
    double a1 = -60 * cache.zeta_3_A * (cache.omega_10_A + 4) + 1680 * cache.zeta_3_V;
    double a2 = 30 * cache.zeta_3_A * (15 * cache.omega_10_A + 32) - 12600 * cache.zeta_3_V + 36 * cache.a_1_par - 72 * cache.a_2_par - 12;
    double a3 = -100 * cache.zeta_3_A * (9 * cache.omega_10_A + 8) + 25200 * cache.zeta_3_V - 48 * cache.a_1_par + 240 * cache.a_2_par;
    double a4 = 525 * cache.zeta_3_A * cache.omega_10_A - 14700 * cache.zeta_3_V - 180 * cache.a_2_par;
    return -u * (a1 + u * (a2 + u * (a3 + u * a4))) / 4 + cache.delta_t_p * (9 * u - 1.5) + cache.delta_t_m * 6 * u + 3 * (cache.delta_t_p + cache.delta_t_m) * log(1 - u);
}

complex_t BKstarDecay::G(double xbar) {
    complex_t alpha = complex_t(cache.s_c, -EPSILON);

    if (xbar < EPSILON)
        return -2. * log(alpha) / 3.;

    complex_t A = (xbar + 2. * alpha) * sqrt(xbar - 4. * alpha) / pow(xbar, 1.5);
    return 10. / 9 - 2. * log(alpha) / 3. + 8. * alpha / (3 * xbar) + 2. * A * log(-(xbar * (1. - A) + 2. * alpha) / (xbar * (1. + A) + 2. * alpha)) / 3.;
}

complex_t BKstarDecay::G_perp() {
    auto iG_perp = [this] (double x) {
        double xbar = 1 - x;
        return phi_perp(x) * G(xbar) / (3 * xbar);
    };

    return c_integrate(iG_perp, 0, 1, 1e-3);
}

complex_t BKstarDecay::G2() {
    double ls = log(cache.s_c);
    double ls2 = ls * ls;
    double ls3 = ls * ls2;
    complex_t a_0 {-833. / 162., -20. * PI / 27.};
    complex_t a_1 = 48. - 5. * PI2 - 36. * ZETA3 + I * (30. * PI - 2. * PI3) + (36. - 9. * PI2 + 6. * PI * I) * ls + (3. + 6. * PI * I) * ls2 + ls3;
    complex_t a_2 = 18. + 2. * PI2 - 2. * PI3 * I + (12. - 6. * PI2) * ls + 6. * PI * I * ls2 + ls3;
    complex_t a_3 = -9. - 14. * PI2 + 112. * PI * I + (182. - 48. * PI * I) * ls - 126. * ls2;
    complex_t g_2 = a_0 + 2. * (cache.s_c * a_1 + cache.s_c * cache.s_c * a_2) / 9. + cache.s_c * cache.s_c * cache.s_c * a_3 / 27. + 8. * PI2 * pow(cache.s_c, 1.5) / 9.;
    
    return 8 * std::log(cache.mu_b / cache.m_b_1S) / 3 + g_2;
}

complex_t BKstarDecay::G8() {
    complex_t g_8 = 11. / 3 - 2 * PI2 / 9 + 2. * I * PI / 3.;
    return -104 * std::log(cache.mu_b / cache.m_b_1S) / 27 + g_8;
}

complex_t BKstarDecay::H_perp() {
    auto iH_perp = [this] (double x) {
        return gv_dga_4(x) * G(1 - x);
    }; 
    return c_integrate(iH_perp, 0, 1, 1e-3);
}

complex_t BKstarDecay::H2() {
    auto iH2_perp = [this] (double x) {
        return h(1 - x) * phi_perp(x);
    }; 

    return c_integrate(iH2_perp, 0, 1, 1e-4);
}

complex_t BKstarDecay::K1() {
    double F_perp = 1 + cache.a_1_perp + cache.a_2_perp;
    return -(cache.C_b[WCoef::C6] + cache.C_b[WCoef::C5] / cache.Nc) * F_perp
           + cache.C_F * cache.alpha_s_mu_b / (4 * cache.Nc * PI) * (
                pow(cache.m_b_mu_b / cache.m_B, 2) * cache.C_b[WCoef::C8] * X_perp() 
              - cache.C_b[WCoef::C2] * ((4 * log(cache.m_b_1S / cache.mu_b) + 2) * F_perp / 3 - G_perp()) 
              + F_perp * log(cache.mu_b / cache.mu_0) * (
                    8. * cache.C_b[WCoef::C3] / 3. 
                  + 4 * cache.n_f * (cache.C_b[WCoef::C4] + cache.C_b[WCoef::C6]) / 3. 
                  - 8. * (cache.Nc * cache.C_b[WCoef::C6] + cache.C_b[WCoef::C5]))
            );
}

complex_t BKstarDecay::K2(int q) {
    complex_t k2 = cache.C_b[WCoef::C4] + cache.C_b[WCoef::C3] / cache.Nc 
                    + cache.C_F * cache.alpha_s_mu_b / (4 * cache.Nc * PI) * (
                        cache.C_b[WCoef::C2] * ((2 - 4 * log(cache.m_b_1S / cache.mu_b)) / 3. - H_perp()) 
                      + log(cache.mu_b / cache.mu_0) * (
                        -44. * cache.C_b[WCoef::C3] / 3. 
                        -4. * cache.n_f * (cache.C_b[WCoef::C4] + cache.C_b[WCoef::C6]) / 3.)
                    );
    if (q == 2) 
        k2 += cache.lambda_hat_u * (cache.C_b[WCoef::C2] + cache.C_b[WCoef::C1] / cache.Nc);
    return k2;
}

double BKstarDecay::X_perp() {
    double cutoff = cache.Lambda_h / cache.m_B;
    return -2 * (1 + 3 * cache.a_1_perp + 6 * cache.a_2_perp) * log(cutoff) - (1 + 11 * cache.a_1_perp + 31 * cache.a_2_perp) + 12 * cutoff * (cache.a_1_perp + 5 * cache.a_2_perp);
}

double BKstarDecay::delta_0() {    
    double prefactor = PI * cache.C_F * cache.alpha_s_mu_h * cache.f_B * cache.f_Ks_perp / (6. * cache.Nc * cache.T1_B_Ks * cache.m_B * cache.lambda_B);
    double H_8 = 3 * (1 - cache.a_1_perp + cache.a_2_perp);
    complex_t a7c = cache.C_b[WCoef::C7] + cache.alpha_s_mu_b * cache.C_F * (cache.C_b[WCoef::C2] * G2() + cache.C_b[WCoef::C8] * G8()) / (4. * PI) 
                    + prefactor * (2. * cache.C8_h * H_8 - cache.C2_h * H2());

    complex_t pref = 4 * PI2 * cache.f_B / (cache.m_b_mu_b * cache.T1_B_Ks * a7c);
    complex_t t1 = cache.f_Ks_perp * K1() / cache.m_b_mu_b;
    complex_t f2 = cache.f_Ks_par * cache.m_Ks / (6 * cache.lambda_B * cache.m_B);
    complex_t bd = -pref * (t1 + f2 * K2(1));
    complex_t bu = 2. * pref * (t1 + f2 * K2(2));
    return real(bd - bu);
}

std::vector<ObservableValue> BKstarDecay::compute_observable(Observables obs) {
    double value;
    switch (obs) {
    case Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA:   
        value = delta_0();
        break;
    default:
        LOG_ERROR("IndexError", "Observable", ObservableMapper::str(obs), "doesn't belong to the decay", DecayMapper::str(this->id));
    }

    return {ObservableValue(ObservableMapper::to_id(obs), value)};
}

std::vector<ObservableValue> BKstarDecay::compute_observable(ObservableId obs) {
    return compute_observable(ObservableMapper::enum_of(obs).value());
}