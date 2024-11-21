#include "B__Kstar_gamma.h"


double Delta0_B__Kstar_gamma::eval() const {
    auto sm_p = Parameters::GetInstance(0); // SM params
    auto flav_p = Parameters::GetInstance(3); // Flavor params
    auto manager = computeWilsons();

    double N = 3;
    double nf = 5;
    double C_F = (N * N - 1) / (2 * N);
    double b0 = 11 - 2 * nf / 3.;

    double lambda_h = 0.5;
    double mu_h = std::sqrt(lambda_h * scale);
    double alphas_mu_b = sm_p->alpha_s(scale);
    double alphas_mu_h = sm_p->alpha_s(mu_h);
    double alphas_1GeV = sm_p->alpha_s(1.);
    double eta = alphas_mu_b / alphas_1GeV;
    
    double m_B = (*flav_p)("FMASS", 521);
    double m_Ks = (*flav_p)("FMASS", 323);
    double m_b = (*sm_p)("MASS", 5);
    double f_B = (*flav_p)("FCONST", 52101);
    double f_Ks = (*flav_p)("FCONST", 32301);
    double f_Ks_perp = (*flav_p)("FCONST", 32302) * pow(eta, FFInput::gamma_n_perp(0, N) / b0);
    double CKM_factor = std::real(std::conj(Parameters::get_c_CKM_entry(12)) * Parameters::get_c_CKM_entry(13) 
                                  / (std::conj(Parameters::get_c_CKM_entry(22)) * Parameters::get_c_CKM_entry(23))); // Vus*Vub/Vcs*Vcb
    double a_1_perp = FFInput::a_1_perp * std::pow(eta, (FFInput::gamma_n_perp(1, N) - C_F) / b0);
    double a_2_perp = FFInput::a_2_perp * std::pow(eta, (FFInput::gamma_n_perp(2, N) - C_F) / b0);
    double a_1_par = FFInput::a_1_par * std::pow(eta, (FFInput::gamma_n_par(1, N) - C_F) / b0);
    double a_2_par = FFInput::a_2_par * std::pow(eta, (FFInput::gamma_n_par(2, N) - C_F) / b0);
    double lambda_B = FFInput::lambda_B / (1 - alphas_mu_h * std::log(std::pow(mu_h, 2)) * 1.8 / (3 * PI));
    
    complex_t C1 = manager->getFullRunCoefficient("BCoefficient", "C1", "NLO") + manager->getFullRunCoefficient("BPrimeCoefficient", "CP1", "NLO");
    complex_t C2 = manager->getFullRunCoefficient("BCoefficient", "C2", "NLO") + manager->getFullRunCoefficient("BPrimeCoefficient", "CP2", "LO");
    complex_t C3 = manager->getFullRunCoefficient("BCoefficient", "C3", "NLO");
    complex_t CP3 = manager->getFullRunCoefficient("BPrimeCoefficient", "CP3", "LO");
    complex_t C4 = manager->getFullRunCoefficient("BCoefficient", "C4", "NLO");
    complex_t CP4 = manager->getFullRunCoefficient("BPrimeCoefficient", "CP4", "LO");
    complex_t C5 = manager->getFullRunCoefficient("BCoefficient", "C5", "NLO");
    complex_t CP5 = manager->getFullRunCoefficient("BPrimeCoefficient", "CP5", "LO");
    complex_t C6 = manager->getFullRunCoefficient("BCoefficient", "C6", "NLO");
    complex_t CP6 = manager->getFullRunCoefficient("BPrimeCoefficient", "CP6", "LO");
    complex_t C7 = manager->getFullRunCoefficient("BCoefficient", "C7", "NLO") + manager->getFullRunCoefficient("BPrimeCoefficient", "CP7", "LO");
    complex_t C8 = manager->getFullRunCoefficient("BCoefficient", "C8", "NLO") + manager->getFullRunCoefficient("BPrimeCoefficient", "CP8", "LO");

    manager->setGroupScale("BCoefficientGroup", mu_h);
    manager->setGroupScale("BPrimeCoefficientGroup", mu_h);
    manager->setGroupScale("BScalarCoefficientGroup", mu_h);
    manager->switchbasis("BCoefficientGroup");

    complex_t C2_h = manager->getFullRunCoefficient("BCoefficient", "C2", "NNLO") + manager->getFullRunCoefficient("BPrimeCoefficient", "CP2", "NNLO");
    complex_t C8_h = manager->getFullRunCoefficient("BCoefficient", "C2",  "NNLO") + manager->getFullRunCoefficient("BPrimeCoefficient", "CP8", "NNLO");
    
    double m_b_mu_b = sm_p->running_mass(m_b, m_b, scale, "pole");
    double m_b_1S = sm_p->get_QCD_masse("mb_1S");
    double m_c_mu_b = sm_p->running_mass((*sm_p)("MASS", 4), (*sm_p)("MASS", 4), scale, "pole");
    double sc = std::pow(m_c_mu_b / m_c_mu_b, 2);

    double mu_0 = scale;
    complex_t r1 = (8. * C3 / 3. + 4 * nf * (C4 + C6) / 3. - 8. * (N * C6 + C5)) * F_perp(a_1_perp, a_2_perp) * std::log(scale / mu_0);
    complex_t r2 = (-44. * C3 / 3. - 4 * nf * (C4 + C6) / 3.) * std::log(scale / mu_0);

    double rb = scale / m_b_1S;
    double cutoff = 1 - lambda_h / m_B;
    double F = F_perp(a_1_perp, a_2_perp);
    complex_t G = G_perp(sc, a_1_perp, a_2_perp);
    complex_t H = H_perp(sc, a_1_par, a_2_par);
    double X = X_perp(a_1_perp, a_2_perp, cutoff);
    complex_t G2 = G2_perp(sc, rb);
    complex_t G8 = G8_perp(rb);
    complex_t H2 = 2 * PI2 * f_B * f_Ks_perp / (3 * N * FFInput::T1_B_Kstar * m_B * lambda_B) * H2_perp(sc, a_1_perp, a_2_perp);
    double H8 = 4 * PI2 * f_B * f_Ks_perp / (3 * N * FFInput::T1_B_Kstar * m_B * lambda_B) * H8_perp(a_1_perp, a_2_perp);

    complex_t a7c = C7 + (alphas_mu_b * (C2 * G2 + C8 * G8) + alphas_mu_h * (C2_h * H2 + C8_h * H8)) * C_F / (4 * PI);
    complex_t K1 = -(C6 + CP6 + (C5 + CP5) / N) * F 
                    + C_F * alphas_mu_b * (std::pow(m_b_mu_b / m_B, 2) * C8 * X - C2 * ((2 - 4 * std::log(rb)) * F / 3 - G) + r1) / (4 * N * PI);
    complex_t K2d = C4 + CP4 + (C3 + CP3) / N + C_F * alphas_mu_b * (C2 * ((2 - 4 * std::log(rb)) / 3. - H) + r2) / (4 * N * PI);
    complex_t K2u = CKM_factor * (C2 + C1 / N) + K2d;

    complex_t pref = 4 * PI2 * f_B / (m_b * FFInput::T1_B_Kstar * a7c);
    complex_t t1 = f_Ks_perp * K1 / m_b;
    complex_t f2 = f_Ks * m_Ks / (6 * lambda_B * m_B);
    complex_t bd = -pref * (t1 + f2 * K2d);
    complex_t bu = 2. * pref * (t1 + f2 * K2u);

    return std::real(bd - bu);
}