#include "Bd_mumu.h"

double BR_Bd_mumu::eval() const {
    auto sm_p = Parameters::GetInstance(ParameterType::SM); // SM params
    auto flav_p = Parameters::GetInstance(ParameterType::FLAVOR); // Flavor params
    auto manager = computeWilsons();

    complex_t C10 = manager->getFullRunCoefficient("BCoefficient", "C10", "NNLO");
    complex_t CQ1 = manager->getFullRunCoefficient("BScalarCoefficient", "CQ1", "NLO");
    complex_t CQ2 = manager->getFullRunCoefficient("BScalarCoefficient", "CQ2", "NLO");

    double G_F = (*sm_p)("SMINPUTS", 2);
    double inv_alpha_em = (*sm_p)("SMINPUTS", 1);
    double V_tbV_td = std::abs(Parameters::get_c_CKM_entry(22) * std::conj(Parameters::get_c_CKM_entry(20))); 
    double m_Bd = (*flav_p)("FMASS", 511);
    double f_Bd = (*flav_p)("FCONST", 51101);
    double life_Bd = (*flav_p)("FLIFE", 511);

    double r = (*sm_p)("MASS", 13) / m_Bd;  // m_mu / m_Bd
    double x = m_Bd / (sm_p->get_QCD_masse("mb_pole") + (*sm_p)("MASS", 1)); // m_Bd / (m_b_pole + m_d)

    return std::pow(G_F * f_Bd * V_tbV_td / inv_alpha_em, 2) / (64 * HBAR) * std::pow(m_Bd, 3) * INV_PI3 * life_Bd * std::sqrt(1 - 4 * r * r) 
        * ((1 - 4 * r * r) * pow(x * std::abs(CQ1), 2) + pow(std::abs(x * CQ2 + 2 * r * C10), 2));
}