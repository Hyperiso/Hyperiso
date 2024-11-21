#include "Bs_mumu.h"


double BR_Bs_mumu::eval() const {
    auto sm_p = Parameters::GetInstance(0); // SM params
    auto flav_p = Parameters::GetInstance(3); // Flavor params
    auto manager = computeWilsons();

    complex_t C10 = manager->getFullRunCoefficient("BCoefficient", "C10", "NNLO");
    complex_t CP10 = manager->getFullRunCoefficient("BPrimeCoefficient", "CP10", "LO");
    complex_t CQ1 = manager->getFullRunCoefficient("BScalarCoefficient", "CQ1", "NLO");
    complex_t CQ2 = manager->getFullRunCoefficient("BScalarCoefficient", "CQ2", "NLO");
    complex_t CPQ1 = manager->getFullRunCoefficient("BPrimeCoefficient", "CPQ1", "LO");
    complex_t CPQ2 = manager->getFullRunCoefficient("BPrimeCoefficient", "CPQ2", "LO");

    double G_F = (*sm_p)("SMINPUTS", 2);
    double inv_alpha_em = (*sm_p)("SMINPUTS", 1);
    double V_tbV_ts = std::abs(Parameters::get_c_CKM_entry(22) * std::conj(Parameters::get_c_CKM_entry(21))); 

    double m_Bs = (*flav_p)("FMASS", 531);
    double f_Bs = (*flav_p)("FCONST", 53101);
    double life_Bs = (*flav_p)("FLIFE", 531);
    
    double r = (*sm_p)("MASS", 13) / m_Bs;  // m_mu / m_Bs
    double x = m_Bs / (sm_p->get_QCD_masse("mb_pole") + (*sm_p)("MASS", 3)); // m_Bs / (m_b_pole + m_s)

    return std::pow(G_F * f_Bs * V_tbV_ts / inv_alpha_em, 2) / (64 * HBAR) * std::pow(m_Bs, 3) * INV_PI3 * life_Bs * std::sqrt(1 - 4 * r * r) 
            * ((1 - 4 * r * r) * pow(x * std::abs(CQ1 - CPQ1), 2) + pow(std::abs(x * (CQ2 - CPQ2) + 2 * r * (C10 - CP10)), 2));
}

double BR_Bs_mumu_untag::eval() const {
    auto sm_p = Parameters::GetInstance(0); // SM params
    auto flav_p = Parameters::GetInstance(3); // Flavor params
    auto manager = computeWilsons();

    complex_t C10 = manager->getFullRunCoefficient("BCoefficient", "C10", "NNLO");
    complex_t CP10 = manager->getFullRunCoefficient("BPrimeCoefficient", "CP10", "LO");
    complex_t CQ1 = manager->getFullRunCoefficient("BScalarCoefficient", "CQ1", "NLO");
    complex_t CQ2 = manager->getFullRunCoefficient("BScalarCoefficient", "CQ2", "NLO");
    complex_t CPQ1 = manager->getFullRunCoefficient("BPrimeCoefficient", "CPQ1", "LO");
    complex_t CPQ2 = manager->getFullRunCoefficient("BPrimeCoefficient", "CPQ2", "LO");

    double G_F = (*sm_p)("SMINPUTS", 2);
    double inv_alpha_em = (*sm_p)("SMINPUTS", 1);
    double V_tbV_ts = std::abs(Parameters::get_c_CKM_entry(22) * std::conj(Parameters::get_c_CKM_entry(21))); 

    double m_Bs = (*flav_p)("FMASS", 531);
    double f_Bs = (*flav_p)("FCONST", 53101);
    double life_Bs = (*flav_p)("FLIFE", 531);
    
    double r = (*sm_p)("MASS", 13) / m_Bs;  // m_mu / m_Bs
    double x = m_Bs / (sm_p->get_QCD_masse("mb_pole") + (*sm_p)("MASS", 3)); // m_Bs / (m_b_pole + m_s)

    auto manager_sm = computeWilsons(0, 2, m_Bs);
    complex_t C10_SM = manager_sm->getFullRunCoefficient("BCoefficient", "C10", "NNLO");
    complex_t S = std::sqrt(1 - 4 * r * r) * x / 2 / r * (CQ1 - CPQ1) / C10_SM;
    complex_t P = (C10 - CP10 + x * (CQ2 - CPQ2) / (2 * r)) / C10_SM;
    double magn_S = std::pow(std::abs(S), 2);
    double magn_P = std::pow(std::abs(P), 2);
    double phi_S = std::arg(S);
    double phi_P = std::arg(P);
    double A = (magn_P * std::cos(2 * phi_P) - magn_S * std::cos(2 * phi_S)) / (magn_P + magn_S);
    double ys_Bs = 0.068;
    double untag_factor = (1 + A * ys_Bs) / (1 - ys_Bs * ys_Bs);

    return untag_factor * std::pow(G_F * f_Bs * V_tbV_ts / inv_alpha_em, 2) / (64 * HBAR) * std::pow(m_Bs, 3) * INV_PI3 * life_Bs * std::sqrt(1 - 4 * r * r) 
            * ((1 - 4 * r * r) * pow(x * std::abs(CQ1 - CPQ1), 2) + pow(std::abs(x * (CQ2 - CPQ2) + 2 * r * (C10 - CP10)), 2));
}
