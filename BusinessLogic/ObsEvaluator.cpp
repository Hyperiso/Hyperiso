#include "ObsEvaluator.h"
// #include "Wilson.h"
// #include "Wilson_susy.h"
// #include "Wilson_THDM.h"
#include "Logger.h"
#include "Math.h"
#include "Parameters.h"
#include "epsilon_calculator.h"
#include "WilsonManager.h"
#include "Wilsonv2.h"
#include "Wilson_THDMv2.h"
#include "Wilson_susyv2.h"

void *ObsEvaluator::computeWilsons(int model, int order, double scale, bool traditional_basis=false) {
    std::cout << "in computewilson1 " << std::endl;
    double m_W = (*Parameters::GetInstance(0))("MASS", 24);
    std::cout << "in computewilson2 " << std::endl;
    // WilsonManager* wm;
    switch (model) {
        case 0:
            auto* manager_sm = CoefficientManager::GetInstance("StandardModel");
            manager_sm->registerCoefficientGroup("BCoefficients", std::make_unique<BCoefficientGroup>(81.0));
            manager_sm->registerCoefficientGroup("ScalarCoefficients", std::make_unique<BScalarCoefficientGroup>(81.0));
            manager_sm->registerCoefficientGroup("PrimeCoefficients", std::make_unique<BPrimeCoefficientGroup>(81.0));
            manager_sm->setQMatch();
            std::cout << "in ComputeWilson " << std::endl;
            switch (order) {
                case 0:
                    manager_sm->setMatchingCoefficient("BCoefficients", "LO");
                    manager_sm->setMatchingCoefficient("ScalarCoefficients", "LO");
                    manager_sm->setMatchingCoefficient("PrimeCoefficients", "LO");
                    // wm = WilsonManager::GetInstance("LO", m_W, std::make_shared<SM_LO_Strategy>());
                    break;
                case 1:
                    manager_sm->setMatchingCoefficient("BCoefficients", "NLO");
                    manager_sm->setMatchingCoefficient("ScalarCoefficients", "NLO");
                    manager_sm->setMatchingCoefficient("PrimeCoefficients", "NLO");
                    // wm = WilsonManager::GetInstance("NLO", m_W, std::make_shared<SM_NLO_Strategy>());
                    break;
                case 2:
                    manager_sm->setMatchingCoefficient("BCoefficients", "NNLO");
                    manager_sm->setMatchingCoefficient("ScalarCoefficients", "NNLO");
                    manager_sm->setMatchingCoefficient("PrimeCoefficients", "NNLO");
                    // wm = WilsonManager::GetInstance("NNLO", m_W, std::make_shared<SM_NNLO_Strategy>());
                    break;
                default:
                    LOG_WARN("Order too high required for SM Wilson coefficients, defaulting to order 2.");
                    // wm = WilsonManager::GetInstance("NNLO", m_W, std::make_shared<SM_NNLO_Strategy>());
                    manager_sm->setMatchingCoefficient("BCoefficients", "NNLO");
                    manager_sm->setMatchingCoefficient("ScalarCoefficients", "NNLO");
                    manager_sm->setMatchingCoefficient("PrimeCoefficients", "NNLO");
            }
            break;
        case 1:
            auto* manager_susy = CoefficientManager::GetInstance("Susy");
            manager_susy->registerCoefficientGroup("BCoefficients", std::make_unique<BCoefficientGroup_SUSY>(81.0));
            manager_susy->registerCoefficientGroup("ScalarCoefficients", std::make_unique<BScalarCoefficientGroup_SUSY>(81.0));
            manager_susy->registerCoefficientGroup("PrimeCoefficients", std::make_unique<BPrimeCoefficientGroup_SUSY>(81.0));
            switch (order) {
                case 0:
                    manager_susy->setMatchingCoefficient("BCoefficients", "LO");
                    manager_susy->setMatchingCoefficient("ScalarCoefficients", "LO");
                    manager_susy->setMatchingCoefficient("PrimeCoefficients", "LO");
                    break;
                case 1:
                    manager_susy->setMatchingCoefficient("BCoefficients", "NLO");
                    manager_susy->setMatchingCoefficient("ScalarCoefficients", "NLO");
                    manager_susy->setMatchingCoefficient("PrimeCoefficients", "NLO");
                    break;
                case 2:
                    manager_susy->setMatchingCoefficient("BCoefficients", "NNLO");
                    manager_susy->setMatchingCoefficient("ScalarCoefficients", "NNLO");
                    manager_susy->setMatchingCoefficient("PrimeCoefficients", "NNLO");
                    break;
                default:
                    LOG_WARN("Order too high required for SUSY Wilson coefficients, defaulting to order 2.");
                    manager_susy->setMatchingCoefficient("BCoefficients", "NNLO");
                    manager_susy->setMatchingCoefficient("ScalarCoefficients", "NNLO");
                    manager_susy->setMatchingCoefficient("PrimeCoefficients", "NNLO");
            }
            break;
        case 2:
            auto* manager_thdm = CoefficientManager::GetInstance("THDM");
            switch (order) {
                case 0:
                    manager_thdm->registerCoefficientGroup("BCoefficients", std::make_unique<BCoefficientGroup_THDM>(81.0));
                    manager_thdm->registerCoefficientGroup("ScalarCoefficients", std::make_unique<BScalarCoefficientGroup_THDM>(81.0));
                    manager_thdm->registerCoefficientGroup("PrimeCoefficients", std::make_unique<BPrimeCoefficientGroup_THDM>(81.0));
                    break;
                case 1:
                    manager_thdm->setMatchingCoefficient("BCoefficients", "NLO");
                    manager_thdm->setMatchingCoefficient("ScalarCoefficients", "NLO");
                    manager_thdm->setMatchingCoefficient("PrimeCoefficients", "NLO");
                    break;
                case 2:
                    manager_thdm->setMatchingCoefficient("BCoefficients", "NNLO");
                    manager_thdm->setMatchingCoefficient("ScalarCoefficients", "NNLO");
                    manager_thdm->setMatchingCoefficient("PrimeCoefficients", "NNLO");
                    break;
                default:
                    LOG_WARN("Order too high required for THDM Wilson coefficients, defaulting to order 2.");
                    manager_thdm->setMatchingCoefficient("BCoefficients", "NNLO");
                    manager_thdm->setMatchingCoefficient("ScalarCoefficients", "NNLO");
                    manager_thdm->setMatchingCoefficient("PrimeCoefficients", "NNLO");
            }
            break;
        default:
            LOG_ERROR("ModelError", "Unknown model requested for Wilson coefficient calculation.");
            // return nullptr;
    }
    // wm->setScale(scale, true, traditional_basis);
    // return 0;
}

complex_t get_c_CKM_entry(int idx) {
    auto p = Parameters::GetInstance(0);
    return complex_t((*p)("RECKM", idx), (*p)("IMCKM", idx));
}

complex_t ObsEvaluator::Evaluate(Observable *o) {
    std::cout << "in Evaluate " << std::endl;
    auto p = Parameters::GetInstance(0);
    std::cout << "in Evaluate2" << std::endl;
    ObsEvaluator::computeWilsons(o->getModel(), o->getOrder(), o->getScale(), o->getWilsonBasis() == 2);
    std::cout << "in Evaluate3 " << std::endl;
    // if (!wm) {
    //     return std::complex<double>(-1);
    // }

    switch (o->getId()) {
        case Observables::BR_BS_MUMU:
            return ObsEvaluator::Bs_mumu(wm, false);
        case Observables::BR_BS_MUMU_UNTAG:
            return ObsEvaluator::Bs_mumu(wm, true);
        case Observables::BR_BD_MUMU:
            return ObsEvaluator::Bd_mumu(wm);
        case Observables::BR_BU_TAUNU:
            return ObsEvaluator::Bu_taunu(o->getModel(), false);
        case Observables::BR_BU_TAUNU_NP_ONLY:
            return ObsEvaluator::Bu_taunu(o->getModel(), true);
        case Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA:
            return 0;
            return ObsEvaluator::Delta_0_B_Kstargamma(wm, o->getScale());
        case Observables::LAST:
            return 0;
        case Observables::FIRST:
            return 0;
        default:
            LOG_ERROR("ValueError", "Unknown observable.");
            return std::complex<double>(-1);
    }
}

complex_t ObsEvaluator::Bs_mumu(WilsonManager* wm, bool untag)
{
    std::cout << "Bs_mummu" << std::endl;
    auto sm_p = Parameters::GetInstance(0); // SM params
    auto flav_p = Parameters::GetInstance(3); // Flavor params

    complex_t C10 = wm->get_full(WilsonCoefficient::C10, 2);
    complex_t CP10 = wm->get_full(WilsonCoefficient::CP10, 0);
    complex_t CQ1 = wm->get_full(WilsonCoefficient::CQ1, 1);
    complex_t CQ2 = wm->get_full(WilsonCoefficient::CQ2, 1);
    complex_t CPQ1 = wm->get_full(WilsonCoefficient::CPQ1, 0);
    complex_t CPQ2 = wm->get_full(WilsonCoefficient::CPQ2, 0);

    double G_F = (*sm_p)("SMINPUTS", 2);
    double inv_alpha_em = (*sm_p)("SMINPUTS", 1);
    inv_alpha_em = 137.;
    G_F = 1.166e-5;
    double V_tbV_ts = std::abs(get_c_CKM_entry(22) * std::conj(get_c_CKM_entry(21))); 

    double m_Bs = (*flav_p)("MASS", 531);
    double f_Bs = flav_p->getFlavorParam(FlavorParamType::DECAY_CONSTANT, "531|1");
    double life_Bs = flav_p->getFlavorParam(FlavorParamType::LIFETIME, "531");
    
    double r = (*sm_p)("MASS", 13) / m_Bs;  // m_mu / m_Bs
    double x = m_Bs / (sm_p->get_QCD_masse("mb_pole") + (*sm_p)("MASS", 3)); // m_Bs / (m_b_pole + m_s)

    double untag_factor = 1;
    if (untag) {
        WilsonManager* wm_sm = ObsEvaluator::computeWilsons(0, 2, m_Bs);
        complex_t C10_SM = wm_sm->get_full(WilsonCoefficient::C10, 2);
        complex_t S = std::sqrt(1 - 4 * r * r) * x / 2 / r * (CQ1 - CPQ1) / C10_SM;
        complex_t P = (C10 - CP10 + x * (CQ2 - CPQ2) / (2 * r)) / C10_SM;
        double magn_S = std::pow(std::abs(S), 2);
        double magn_P = std::pow(std::abs(P), 2);
        double phi_S = std::arg(S);
        double phi_P = std::arg(P);
        double A = (magn_P * std::cos(2 * phi_P) - magn_S * std::cos(2 * phi_S)) / (magn_P + magn_S);
        double ys_Bs = 0.068;
        untag_factor = (1 + A * ys_Bs) / (1 - ys_Bs * ys_Bs);
    }

    return untag_factor * std::pow(G_F * f_Bs * V_tbV_ts / inv_alpha_em, 2) / (64 * HBAR) * std::pow(m_Bs, 3) * INV_PI3 * life_Bs * std::sqrt(1 - 4 * r * r) 
            * ((1 - 4 * r * r) * pow(x * std::abs(CQ1 - CPQ1), 2) + pow(std::abs(x * (CQ2 - CPQ2) + 2 * r * (C10 - CP10)), 2));
}

complex_t ObsEvaluator::Bd_mumu(WilsonManager* wm) {
    std::cout << "Bdmumu " << std::endl;
    auto sm_p = Parameters::GetInstance(0); // SM params
    auto flav_p = Parameters::GetInstance(3); // Flavor params

    complex_t C10 = wm->get_full(WilsonCoefficient::C10, 2);
    complex_t CQ1 = wm->get_full(WilsonCoefficient::CQ1, 1);
    complex_t CQ2 = wm->get_full(WilsonCoefficient::CQ2, 1);

    std::cout << "C1000000" << std::real(C10) << std::endl;
    std::cout << "CQ111111" << std::real(CQ1) << std::endl;
    std::cout << "CQ222222" << std::real(CQ2) << std::endl;

    double G_F = (*sm_p)("SMINPUTS", 2);
    double inv_alpha_em = (*sm_p)("SMINPUTS", 1);
    double V_tbV_td = std::abs(get_c_CKM_entry(33) * std::conj(get_c_CKM_entry(31))); 
    double m_Bd = (*flav_p)("MASS", 511);
    double f_Bd = flav_p->getFlavorParam(FlavorParamType::DECAY_CONSTANT, "511|1");
    double life_Bd = flav_p->getFlavorParam(FlavorParamType::LIFETIME, "511");

    std::cout << "G_F : " << G_F << std::endl;
    std::cout << "inv_alpha_em : " << inv_alpha_em << std::endl;
    std::cout << "V_tbV_td : " << V_tbV_td << std::endl;
    std::cout << "m_Bd : " << m_Bd << std::endl;
    std::cout << "f_Bd : " << f_Bd << std::endl;
    std::cout << "life_Bd : " << life_Bd<< std::endl;

    double r = (*sm_p)("MASS", 13) / m_Bd;  // m_mu / m_Bd
    double x = m_Bd / ((*sm_p)("SMINPUTS", 5) + (*sm_p)("MASS", 2)); // m_Bd / (m_b_pole + m_d)

    return std::pow(G_F * f_Bd * V_tbV_td / inv_alpha_em, 2) / (64 * HBAR) * std::pow(m_Bd, 3) * INV_PI3 * life_Bd * std::sqrt(1 - 4 * r * r) 
        * ((1 - 4 * r * r) * pow(x * std::abs(CQ1), 2) + pow(std::abs(x * CQ2 + 2 * r * C10), 2));
}

complex_t ObsEvaluator::Bu_taunu(int model, bool np_only) {
    std::cout << "Bu_taunu " << std::endl;
    auto sm_p = Parameters::GetInstance(0); // SM params
    auto flav_p = Parameters::GetInstance(3); // Flavor params
    
    double m_B = (*flav_p)("MASS", 521);
    double life_B = (*flav_p)("FLIFE", 521);
    double f_B = flav_p->getFlavorParam(FlavorParamType::DECAY_CONSTANT, "521|1");
    double m_tau = (*sm_p)("MASS", 15);
    double V_ub = std::abs(get_c_CKM_entry(13)); 
    double G_F = (*sm_p)("SMINPUTS", 2);
    
    double BR_SM = std::pow(G_F * f_B * V_ub * m_tau * (1 - std::pow(m_tau / m_B, 2)), 2) * life_B * m_B;
    
    double np_fact = 1; 
    if (model == 1) {
        double m_Hp = (*Parameters::GetInstance(model))("MASS", 37);
        double tan_b = (*Parameters::GetInstance(model))("HMIX", 2);
        double eps_0 = EpsilonCalculator::GetInstance()->epsilon_0();
        np_fact = std::pow(1 - std::pow(m_B * tan_b / m_Hp, 2) / (1 + eps_0 * tan_b), 2);
    } else if (model == 2) {
        double m_Hp = (*Parameters::GetInstance(model))("MASS", 37);
        double tan_b = (*Parameters::GetInstance(model))("HMIX", 2);
        double l_bb = (*Parameters::GetInstance(model))("YD", 33);
        double l_tt = (*Parameters::GetInstance(model))("YU", 33);
        np_fact = std::pow(1 - std::pow(m_B / m_Hp, 2) * l_bb * l_tt, 2);
    }

    return np_only ? np_fact : np_fact * BR_SM;
}

complex_t ObsEvaluator::Delta_0_B_Kstargamma(WilsonManager* wm, double mu_b) {
    std::cout << "Delta_0 " << std::endl;
    auto sm_p = Parameters::GetInstance(0); // SM params
    auto flav_p = Parameters::GetInstance(3); // Flavor params

    double N = 3;
    double nf = 5;
    double C_F = (N * N - 1) / (2 * N);
    double b0 = 11 - 2 * nf / 3.;

    double lambda_h = 0.5;
    double mu_h = std::sqrt(lambda_h * mu_b);
    double alphas_mu_b = sm_p->alpha_s(mu_b);
    double alphas_mu_h = sm_p->alpha_s(mu_h);
    double alphas_1GeV = sm_p->alpha_s(1.);
    double eta = alphas_mu_b / alphas_1GeV;
    
    double m_B = (*flav_p)("MASS", 521);
    double m_Ks = (*flav_p)("MASS", 323);
    double m_b = (*sm_p)("MASS", 5);
    double f_B = flav_p->getFlavorParam(FlavorParamType::DECAY_CONSTANT, "521|1");
    double f_Ks = flav_p->getFlavorParam(FlavorParamType::DECAY_CONSTANT, "323|1");
    double f_Ks_perp = flav_p->getFlavorParam(FlavorParamType::DECAY_CONSTANT, "323|2") * pow(eta, FFInput::gamma_n_perp(0, N) / b0);
    double CKM_factor = std::real(std::conj(get_c_CKM_entry(12)) * get_c_CKM_entry(13) / (std::conj(get_c_CKM_entry(22)) * get_c_CKM_entry(23))); // Vus*Vub/Vcs*Vcb
    double a_1_perp = FFInput::a_1_perp * std::pow(eta, (FFInput::gamma_n_perp(1, N) - C_F) / b0);
    double a_2_perp = FFInput::a_2_perp * std::pow(eta, (FFInput::gamma_n_perp(2, N) - C_F) / b0);
    double a_1_par = FFInput::a_1_par * std::pow(eta, (FFInput::gamma_n_par(1, N) - C_F) / b0);
    double a_2_par = FFInput::a_2_par * std::pow(eta, (FFInput::gamma_n_par(2, N) - C_F) / b0);
    double lambda_B = FFInput::lambda_B / (1 - alphas_mu_h * std::log(std::pow(mu_h, 2)) * 1.8 / (3 * PI));
    
    complex_t C1 = wm->get_full(WilsonCoefficient::C1, 1) + wm->get_full(WilsonCoefficient::CP1, 1);
    complex_t C2 = wm->get_full(WilsonCoefficient::C2, 1) + wm->get_full(WilsonCoefficient::CP2, 0);
    complex_t C3 = wm->get_full(WilsonCoefficient::C3, 1);
    complex_t CP3 = wm->get_full(WilsonCoefficient::CP3, 0);
    complex_t C4 = wm->get_full(WilsonCoefficient::C4, 1);
    complex_t CP4 = wm->get_full(WilsonCoefficient::CP4, 0);
    complex_t C5 = wm->get_full(WilsonCoefficient::C5, 1);
    complex_t CP5 = wm->get_full(WilsonCoefficient::CP5, 0);
    complex_t C6 = wm->get_full(WilsonCoefficient::C6, 1);
    complex_t CP6 = wm->get_full(WilsonCoefficient::CP6, 0);
    complex_t C7 = wm->get_full(WilsonCoefficient::C7, 1) + wm->get_full(WilsonCoefficient::CP7, 0);
    complex_t C8 = wm->get_full(WilsonCoefficient::C8, 1) + wm->get_full(WilsonCoefficient::CP8, 0);

    wm->setScale(mu_h, true, true);
    complex_t C2_h = wm->get_full(WilsonCoefficient::C2, 2) + wm->get_full(WilsonCoefficient::CP2, 2);
    complex_t C8_h = wm->get_full(WilsonCoefficient::C8, 2) + wm->get_full(WilsonCoefficient::CP8, 2);
    
    double m_b_mu_b = sm_p->running_mass(m_b, m_b, mu_b, "pole");
    double m_b_1S = sm_p->get_QCD_masse("mb_1S");
    double m_c_mu_b = sm_p->running_mass((*sm_p)("MASS", 4), (*sm_p)("MASS", 4), mu_b, "pole");
    double sc = std::pow(m_c_mu_b / m_c_mu_b, 2);

    double mu_0 = mu_b;
    complex_t r1 = (8. * C3 / 3. + 4 * nf * (C4 + C6) / 3. - 8. * (N * C6 + C5)) * F_perp(a_1_perp, a_2_perp) * std::log(mu_b / mu_0);
    complex_t r2 = (-44. * C3 / 3. - 4 * nf * (C4 + C6) / 3.) * std::log(mu_b / mu_0);

    double rb = mu_b / m_b_1S;
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