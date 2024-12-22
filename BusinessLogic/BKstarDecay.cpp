#include "BKstarDecay.h"

double BKstarDecay::alpha_s(double mu) {
    auto sm_p = Parameters::GetInstance(ParameterType::SM);
    return sm_p->alpha_s(mu); 
}

double BKstarDecay::beta_0(double mu) {
    auto QCD = Parameters::GetInstance(ParameterType::SM)->QCDaddress();
    auto [b0, _, __] = QCD->getBetas(QCD->getNf(mu));
    return b0;
}

double BKstarDecay::sc(double m_c_m_c, double m_b_m_b) {
    auto sm_p = Parameters::GetInstance(ParameterType::SM);
    double m_b_mu_b = sm_p->running_mass(m_b_m_b, m_b_m_b, this->winfo.hadronic_scale, "pole");
    double m_c_mu_b = sm_p->running_mass(m_c_m_c, m_c_m_c, this->winfo.hadronic_scale, "pole");
    return std::pow(m_c_mu_b / m_b_mu_b, 2);
}

double BKstarDecay::run(double initial_value, double eta, double gamma, double beta) {
    return initial_value * pow(eta, gamma / beta);
}

double BKstarDecay::a_n_perp(int n, double a_1_gev, double beta_0, double eta) {
    double C_F = Parameters::GetInstance(ParameterType::SM)->QCDaddress()->C_F;
    double gamma = 4 * C_F * (psi(n + 1) + GAMMA - 1. + 1. / (n + 1));
    return run(a_1_gev, eta, gamma, beta_0);
}

double BKstarDecay::a_n_par(int n, double a_1_gev, double beta_0, double eta) {
    double C_F = Parameters::GetInstance(ParameterType::SM)->QCDaddress()->C_F;
    double gamma = 4 * C_F * (psi(n + 2) + GAMMA - .75 - 1. /(2 * (n + 1) * (n + 2)));
    return run(a_1_gev, eta, gamma, beta_0);
}

double BKstarDecay::lambda_B(double lam_B_1_gev, double mu_h, double alpha_s_mu_h) {
    return lam_B_1_gev / (1 - alpha_s_mu_h * std::log(std::pow(mu_h, 2)) * 1.8 / (3 * PI));;
}

double BKstarDecay::f_Ks_perp(double f_1_gev, double beta_0, double eta) {
    double C_F = Parameters::GetInstance(ParameterType::SM)->QCDaddress()->C_F;
    return run(f_1_gev, eta, C_F, beta_0);
}

complex_t BKstarDecay::h(double s, double u) {
    complex_t rt = std::sqrt((u - 4 * s + I * EPSILON) / u);
    return (4 * s * (CLi2(2. / (1. - rt)) + CLi2(2. / (1. + rt))) / u - 2.) / u;
}

complex_t BKstarDecay::g_2(double s) {
    double ls = std::log(s);
    double ls2 = ls * ls;
    double ls3 = ls * ls2;
    complex_t a_0 = -833 / 162. - 20 * PI * I / 27.;
    complex_t a_1 = 48 - 5 * PI2 - 36 * ZETA3 + I * (30 * PI - 2 * PI3) + (36 - 9 * PI2 + 6 * PI * I) * ls + (3. + 6 * PI * I) * ls2 + ls3;
    complex_t a_2 = 18 + 2 * PI2 - 2 * PI3 * I + (12 - 6 * PI2) * ls + 6 * PI * I * ls2 + ls3;
    complex_t a_3 = -9 - 14 * PI2 + 112 * PI * I + (182. - 48 * PI * I) * ls - 126 * ls2;
    return a_0 + 2. * (s * a_1 + s * s * a_2) / 9. + s * s * s * a_3 / 27. + 8 * PI2 * std::pow(s, 1.5) / 9;
}

double BKstarDecay::phi_perp(double a1, double a2, double u) {
    double ubar = 1 - u;
    double xi = u - ubar;
    return 6 * u * ubar * (1 + 3 * a1 * xi + 1.5 * a2 * (5 * xi * xi - 1));
}

double BKstarDecay::gv_dga_4(double a1p, double a2p, double z3a, double z3v, double w10a, double dtp, double dtm, double u) {
    double a1 = -60 * z3a * (w10a + 4) + 1680 * z3v;
    double a2 = 30 * z3a * (15 * w10a + 32) - 12600 * z3v + 36 * a1p - 72 * a2p - 12;
    double a3 = -100 * z3a * (9 * w10a + 8) + 25200 * z3v - 48 * a1p + 240 * a2p;
    double a4 = 525 * z3a * w10a - 14700 * z3v - 180 * a2p;

    return -u * (a1 + u * (a2 + u * (a3 + u * a4))) / 4 + dtp * (9 * u - 1.5) + dtm * 6 * u + 3 * (dtp + dtm) * std::log(1 - u);
}

complex_t BKstarDecay::G(double s, double xbar) {
    complex_t alpha = complex_t(s, -EPSILON);

    if (xbar < EPSILON)
        return -2. * std::log(alpha) / 3.;

    complex_t A = (xbar + 2. * alpha) * std::sqrt(xbar - 4. * alpha) / std::pow(xbar, 1.5);
    return 10. / 9 - 2. * std::log(alpha) / 3. + 8. * alpha / (3 * xbar) + 2. * A * std::log(-(xbar * (1. - A) + 2. * alpha) / (xbar * (1. + A) + 2. * alpha)) / 3.;
}

double BKstarDecay::F_perp(double a1, double a2) {
    return 1 + a1 + a2;
}

complex_t BKstarDecay::G_perp(double s, double a1, double a2) {
    auto iG_perp = [this, s, a1, a2] (double x) {
        double xbar = 1 - x;
        return phi_perp(a1, a2, x) * G(s, xbar) / (3 * xbar);
    };

    return c_integrate(iG_perp, 0, 1, 1e-3);
}

complex_t BKstarDecay::G2(double s) {
    double rb = this->winfo.hadronic_scale / Parameters::GetInstance(ParameterType::SM)->get_QCD_masse("mb_1S"); 
    return 8 * std::log(rb) / 3 + g_2(s);
}

complex_t BKstarDecay::G8() {
    double rb = this->winfo.hadronic_scale / Parameters::GetInstance(ParameterType::SM)->get_QCD_masse("mb_1S"); 
    complex_t g_8 = 11 / 3 - 2 * PI2 / 9 + 2. * I * PI / 3.;
    return -104 * std::log(rb) / 27 + g_8;
}

complex_t BKstarDecay::H_perp(double s, double a1par, double a2par, double z3a, double z3v, double w10a, double dtp, double dtm) {
    auto iH_perp = [this, s, a1par, a2par, z3a, z3v, w10a, dtp, dtm] (double x) {
        return gv_dga_4(a1par, a2par, z3a, z3v, w10a, dtp, dtm, x) * G(s, 1 - x);
    }; 

    return c_integrate(iH_perp, 0, 1, 1e-3);
}

complex_t BKstarDecay::H2(double s, double a1, double a2) {
    auto iH2_perp = [this, s, a1, a2] (double x) {
        double xbar = 1 - x;
        return -h(s, xbar) * phi_perp(a1, a2, x);
    }; 

    return c_integrate(iH2_perp, 0, 1, 1e-3);
}

double BKstarDecay::H8(double a1, double a2) {
    return 3 * (1 - a1 + a2);
}

complex_t BKstarDecay::a7c_h(double mu_h,
                             double alpha_s_mu_h,
                             double f_B,
                             double f_Ks_perp,
                             double T1,
                             double m_B,
                             double lambda_B,
                             complex_t h2,
                             complex_t h8) {
    int N = Parameters::GetInstance(ParameterType::SM)->QCDaddress()->Nc;
    double C_F = Parameters::GetInstance(ParameterType::SM)->QCDaddress()->C_F;
    auto manager = compute_wilsons();
    manager->setGroupScale(GroupMapper::str(WilsonGroups::BCoefficients), mu_h);
    manager->setGroupScale(GroupMapper::str(WilsonGroups::BPrimeCoefficients), mu_h);
    manager->switchbasis(GroupMapper::str(WilsonGroups::BCoefficients));
    complex_t C2_h = manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BCoefficients), "C2", "NNLO") + manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BPrimeCoefficients), "CP2", "LO");
    complex_t C8_h = manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BCoefficients), "C8",  "NNLO") + manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BPrimeCoefficients), "CP8", "LO");

    return PI * C_F * alpha_s_mu_h * f_B * f_Ks_perp * (2. * C8_h * h8 - C2_h * h2) / (6 * N * T1 * m_B * lambda_B);
}

complex_t BKstarDecay::a7c_b(double alpha_s_mu_b, complex_t g2, complex_t g8) {
    auto manager = compute_wilsons();
    complex_t C2 = manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BCoefficients), "C2", "NNLO") + manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BPrimeCoefficients), "CP2", "LO");
    complex_t C7 = manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BCoefficients), "C7", "NNLO") + manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BPrimeCoefficients), "CP7", "LO");
    complex_t C8 = manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BCoefficients), "C8",  "NNLO") + manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BPrimeCoefficients), "CP8", "LO");
    double C_F = Parameters::GetInstance(ParameterType::SM)->QCDaddress()->C_F;
    return C7 + alpha_s_mu_b * C_F * (C2 * g2 + C8 * g8) / (4 * PI);
}

complex_t BKstarDecay::r1(double mu_0, double F_p) {
    if (fpeq(mu_0, -1.)) mu_0 = winfo.hadronic_scale;
    auto manager = compute_wilsons();
    complex_t C3 = manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BCoefficients), "C3", "NNLO") + manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BPrimeCoefficients), "CP3", "LO");
    complex_t C4 = manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BCoefficients), "C4", "NNLO") + manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BPrimeCoefficients), "CP4", "LO");
    complex_t C5 = manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BCoefficients), "C5",  "NNLO") + manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BPrimeCoefficients), "CP5", "LO");
    complex_t C6 = manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BCoefficients), "C6",  "NNLO") + manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BPrimeCoefficients), "CP6", "LO");
    double N = Parameters::GetInstance(ParameterType::SM)->QCDaddress()->Nc;
    double nf = Parameters::GetInstance(ParameterType::SM)->QCDaddress()->getNf(winfo.hadronic_scale);
    return (8. * C3 / 3. + 4 * nf * (C4 + C6) / 3. - 8. * (N * C6 + C5)) * F_p * std::log(winfo.hadronic_scale / mu_0);
}

complex_t BKstarDecay::r2(double mu_0) {
    if (fpeq(mu_0, -1.)) mu_0 = winfo.hadronic_scale;
    auto manager = compute_wilsons();
    complex_t C3 = manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BCoefficients), "C3", "NNLO") + manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BPrimeCoefficients), "CP3", "LO");
    complex_t C4 = manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BCoefficients), "C4", "NNLO") + manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BPrimeCoefficients), "CP4", "LO");
    complex_t C6 = manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BCoefficients), "C6",  "NNLO") + manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BPrimeCoefficients), "CP6", "LO");
    double nf = Parameters::GetInstance(ParameterType::SM)->QCDaddress()->getNf(winfo.hadronic_scale);
    return (-44. * C3 / 3. - 4 * nf * (C4 + C6) / 3.) * std::log(winfo.hadronic_scale / mu_0);
}

complex_t BKstarDecay::K1(double mb_mb,
                          double m_B,
                          double alpha_s_mu_b,
                          double F_p,
                          complex_t G_p,
                          complex_t X_p,
                          complex_t r1)
{
    auto manager = compute_wilsons();
    complex_t C2 = manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BCoefficients), "C2", "NNLO") + manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BPrimeCoefficients), "CP2", "LO");
    complex_t C5 = manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BCoefficients), "C5",  "NNLO") + manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BPrimeCoefficients), "CP5", "LO");
    complex_t C6 = manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BCoefficients), "C6",  "NNLO") + manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BPrimeCoefficients), "CP6", "LO");
    complex_t C8 = manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BCoefficients), "C8", "NNLO") + manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BPrimeCoefficients), "CP8", "LO");
    double N = Parameters::GetInstance(ParameterType::SM)->QCDaddress()->Nc;
    double C_F = Parameters::GetInstance(ParameterType::SM)->QCDaddress()->C_F;

    return -(C6 + C5 / N) * F_p + C_F * alpha_s_mu_b / (4 * N * PI) * (std::pow(mb_mb / m_B, 2) * C8 * X_p - C2 * ((4 * std::log(mb_mb / winfo.hadronic_scale) + 2) * F_p / 3 - G_p) + r1);
}

complex_t BKstarDecay::K2d(double mb_mb,
                           double alpha_s_mu_b,
                           complex_t H_p,
                           complex_t r2)
{
    auto manager = compute_wilsons();
    complex_t C2 = manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BCoefficients), "C2", "NNLO") + manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BPrimeCoefficients), "CP2", "LO");
    complex_t C3 = manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BCoefficients), "C3",  "NNLO") + manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BPrimeCoefficients), "CP3", "LO");
    complex_t C4 = manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BCoefficients), "C4",  "NNLO") + manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BPrimeCoefficients), "CP4", "LO");
    double N = Parameters::GetInstance(ParameterType::SM)->QCDaddress()->Nc;
    double C_F = Parameters::GetInstance(ParameterType::SM)->QCDaddress()->C_F;
    return C4 + C3 / N + C_F * alpha_s_mu_b * (C2 * ((2 - 4 * std::log(mb_mb / winfo.hadronic_scale)) / 3. - H_p) + r2) / (4 * N * PI);
}

complex_t BKstarDecay::ckm_factor(double Vus_r,
                                  double Vus_i,
                                  double Vub_r,
                                  double Vub_i,
                                  double Vcs_r,
                                  double Vcs_i,
                                  double Vcb_r,
                                  double Vcb_i)
{
    complex_t V_us = complex_t(Vus_r, Vus_i);
    complex_t V_ub = complex_t(Vub_r, Vub_i);
    complex_t V_cs = complex_t(Vcs_r, Vcs_i);
    complex_t V_cb = complex_t(Vcb_r, Vcb_i);
    return std::conj(V_us) * V_ub / (std::conj(V_cs) * V_cb);
}

complex_t BKstarDecay::K2u(complex_t ckm, complex_t K2d) {
    auto manager = compute_wilsons();
    complex_t C1 = manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BCoefficients), "C1", "NNLO") + manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BPrimeCoefficients), "CP1", "LO");
    complex_t C2 = manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BCoefficients), "C2",  "NNLO") + manager->getFullRunCoefficient(GroupMapper::str(WilsonGroups::BPrimeCoefficients), "CP2", "LO");
    double N = Parameters::GetInstance(ParameterType::SM)->QCDaddress()->Nc;
    return ckm * (C2 + C1 / N) + K2d;
}

double BKstarDecay::delta_0(double f_B,
                            double mb_mb,
                            double T1,
                            double f_Ks_perp,
                            double f_Ks_par,
                            double m_Ks,
                            double m_B,
                            double lambda_B,
                            complex_t a7c,
                            complex_t K1,
                            complex_t K2d,
                            complex_t K2u)
{
    complex_t pref = 4 * PI2 * f_B / (mb_mb * T1 * a7c);
    complex_t t1 = f_Ks_perp * K1 / mb_mb;
    complex_t f2 = f_Ks_par * m_Ks / (6 * lambda_B * m_B);
    complex_t bd = -pref * (t1 + f2 * K2d);
    complex_t bu = 2. * pref * (t1 + f2 * K2u);
    return std::real(bd - bu);
}

double BKstarDecay::X_perp(double a1, double a2, double m_B, double Lambda_h) {
    double cutoff = Lambda_h / m_B;
    return -2 * (1 + 3 * a1 + 6 * a2) * std::log(cutoff) - (1 + 11 * a1 + 31 * a2) + 12 * cutoff * (a1 + 5 * a2);
}

void BKstarDecay::build_op_tree() {
    // Formfactors and decay-specific parameters
    auto a_1_perp   = std::make_shared<ParameterNode>(ParamId(ParameterType::FF, "B_Ks", 1));
    auto a_2_perp   = std::make_shared<ParameterNode>(ParamId(ParameterType::FF, "B_Ks", 2));
    auto a_1_par    = std::make_shared<ParameterNode>(ParamId(ParameterType::FF, "B_Ks", 3));
    auto a_2_par    = std::make_shared<ParameterNode>(ParamId(ParameterType::FF, "B_Ks", 4));
    auto zeta_3_A   = std::make_shared<ParameterNode>(ParamId(ParameterType::FF, "B_Ks", 5));
    auto zeta_3_V   = std::make_shared<ParameterNode>(ParamId(ParameterType::FF, "B_Ks", 6));
    auto w_10_A     = std::make_shared<ParameterNode>(ParamId(ParameterType::FF, "B_Ks", 7));
    auto delta_t_p  = std::make_shared<ParameterNode>(ParamId(ParameterType::FF, "B_Ks", 8));
    auto delta_t_m  = std::make_shared<ParameterNode>(ParamId(ParameterType::FF, "B_Ks", 9));
    auto lambda_B   = std::make_shared<ParameterNode>(ParamId(ParameterType::FF, "B_Ks", 10));
    auto T1_B_Ks    = std::make_shared<ParameterNode>(ParamId(ParameterType::FF, "B_Ks", 11));
    auto Lambda_h   = std::make_shared<ParameterNode>(ParamId(ParameterType::FF, "B_Ks", 12));
    auto mu_0       = std::make_shared<ParameterNode>(ParamId(ParameterType::FF, "B_Ks", 13));
    
    // Flavor parameters
    auto f_Ks_par   = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FCONST", 32301));
    auto f_Ks_perp  = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FCONST", 32302));
    auto f_B        = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FCONST", 52101));
    auto m_B        = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FMASS", 521));
    auto m_Ks       = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FMASS", 323));

    // SM parameters
    auto V_us_r = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "RECKM", 01));
    auto V_us_i = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "IMCKM", 01));
    auto V_cs_r = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "RECKM", 11));
    auto V_cs_i = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "IMCKM", 11));
    auto V_ub_r = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "RECKM", 02));
    auto V_ub_i = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "IMCKM", 02));
    auto V_cb_r = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "RECKM", 12));
    auto V_cb_i = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "IMCKM", 12));
    auto m_c    = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 4));
    auto m_b    = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 5));

    // Operator nodes
    auto alpha_s_mu_b   = std::make_shared<OperatorNode>("alpha_s_mu_b",    [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->alpha_s(this->winfo.hadronic_scale); });
    auto alpha_s_1_gev  = std::make_shared<OperatorNode>("alpha_s_1_gev",   [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->alpha_s(1.); });
    auto eta            = std::make_shared<OperatorNode>("eta",             [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] / values[1]; });
    eta->addChildren({alpha_s_mu_b, alpha_s_1_gev});
    auto beta_0         = std::make_shared<OperatorNode>("beta_0",          [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->beta_0(this->winfo.hadronic_scale); });
    auto mu_h           = std::make_shared<OperatorNode>("mu_h",            [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return std::sqrt(values[0] * this->winfo.hadronic_scale); });
    mu_h->addChild(Lambda_h);
    auto alpha_s_mu_h   = std::make_shared<OperatorNode>("alpha_s_mu_h",    [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->alpha_s(values[0]); });
    alpha_s_mu_h->addChild(mu_h);
    auto lambda_B_mu_h  = std::make_shared<OperatorNode>("lambda_B_mu_h",   [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->lambda_B(values[0], values[1], values[2]); });
    lambda_B_mu_h->addChildren({lambda_B, mu_h, alpha_s_mu_h});
    auto f_Ks_perp_mu_b = std::make_shared<OperatorNode>("f_Ks_perp_mu_b",  [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->f_Ks_perp(values[0], values[1], values[2]); });
    f_Ks_perp_mu_b->addChildren({f_Ks_perp, beta_0, eta});
    auto a_1_perp_mu_b  = std::make_shared<OperatorNode>("a_1_perp_mu_b",   [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->a_n_perp(1, values[0], values[1], values[2]); });
    a_1_perp_mu_b->addChildren({a_1_perp, beta_0, eta});
    auto a_2_perp_mu_b  = std::make_shared<OperatorNode>("a_2_perp_mu_b",   [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->a_n_perp(2, values[0], values[1], values[2]); });
    a_2_perp_mu_b->addChildren({a_2_perp, beta_0, eta});
    auto a_1_par_mu_b   = std::make_shared<OperatorNode>("a_1_par_mu_b",    [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->a_n_par(1, values[0], values[1], values[2]); });
    a_1_par_mu_b->addChildren({a_1_par, beta_0, eta});
    auto a_2_par_mu_b   = std::make_shared<OperatorNode>("a_2_par_mu_b",    [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->a_n_par(2, values[0], values[1], values[2]); });
    a_2_par_mu_b->addChildren({a_2_par, beta_0, eta});
    auto s_c            = std::make_shared<OperatorNode>("s_c",             [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->sc(values[0], values[1]); });
    s_c->addChildren({m_c, m_b});
    auto H2             = std::make_shared<OperatorNode>("H2",              [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->H2(values[0], values[1], values[2]); });
    H2->addChildren({s_c, a_1_perp_mu_b, a_2_perp_mu_b});
    auto H8             = std::make_shared<OperatorNode>("H8",              [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->H8(values[0], values[1]); });
    H8->addChildren({a_1_perp_mu_b, a_2_perp_mu_b});
    auto G2             = std::make_shared<OperatorNode>("G2",              [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->G2(values[0]); });
    G2->addChildren({s_c});
    auto G8             = std::make_shared<OperatorNode>("G8",              [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->G8(); });
    auto F_perp         = std::make_shared<OperatorNode>("F_perp",          [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->F_perp(values[0], values[1]); });
    F_perp->addChildren({a_1_perp_mu_b, a_2_perp_mu_b});
    auto G_perp         = std::make_shared<OperatorNode>("G_perp",          [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->G_perp(values[0], values[1], values[2]); });
    G_perp->addChildren({s_c, a_1_perp_mu_b, a_2_perp_mu_b});
    auto H_perp         = std::make_shared<OperatorNode>("H_perp",          [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->H_perp(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7]); });
    H_perp->addChildren({s_c, a_1_par_mu_b, a_2_par_mu_b, zeta_3_A, zeta_3_V, w_10_A, delta_t_p, delta_t_m});
    auto X_perp         = std::make_shared<OperatorNode>("X_perp",          [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->X_perp(values[0], values[1], values[2], values[3]); });
    X_perp->addChildren({a_1_perp_mu_b, a_2_perp_mu_b, m_B, Lambda_h});
    auto a7c_h          = std::make_shared<OperatorNode>("a7c_h",           [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->a7c_h(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8]); });
    a7c_h->addChildren({mu_h, alpha_s_mu_h, f_B, f_Ks_perp_mu_b, T1_B_Ks, lambda_B_mu_h, H2, H8});
    auto a7c_b          = std::make_shared<OperatorNode>("a7c_b",           [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->a7c_b(values[0], values[1], values[2]); });
    a7c_b->addChildren({alpha_s_mu_b, G2, G8});
    auto a7c            = std::make_shared<OperatorNode>("a7c",             [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] + values[1]; });
    a7c->addChildren({a7c_b, a7c_h});
    auto r1             = std::make_shared<OperatorNode>("r1",              [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->r1(values[0], values[1]); });
    r1->addChildren({mu_0, F_perp});
    auto r2             = std::make_shared<OperatorNode>("r2",              [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->r2(values[0]); });
    r2->addChildren({mu_0});
    auto K1             = std::make_shared<OperatorNode>("K1",              [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K1(values[0], values[1], values[2], values[3], values[4], values[5], values[6]); });
    K1->addChildren({m_b, m_B, alpha_s_mu_b, F_perp, G_perp, X_perp, r1});
    auto K2d            = std::make_shared<OperatorNode>("K2d",             [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K2d(values[0], values[1], values[2], values[3]); });
    K2d->addChildren({m_b, alpha_s_mu_b, H_perp, r2});
    auto ckm            = std::make_shared<OperatorNode>("ckm",             [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->ckm_factor(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7]); });
    ckm->addChildren({V_us_r, V_us_i, V_ub_r, V_ub_i, V_cs_r, V_cs_i, V_cb_r, V_cb_i});
    auto K2u            = std::make_shared<OperatorNode>("K2u",             [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K2u(values[0], values[1]); });
    K2u->addChildren({ckm, K2d});
    auto delta_0        = std::make_shared<OperatorNode>("delta_0",         [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->delta_0(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9], values[10], values[11]); });
    delta_0->addChildren({f_B, m_b, T1_B_Ks, f_Ks_perp_mu_b, f_Ks_par, m_Ks, m_B, lambda_B_mu_h, a7c, K1, K2d, K2u});

    roots.emplace(Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA, delta_0);
}