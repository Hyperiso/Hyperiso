#include "Bu_taunu.h"


double BR_Bu_taunu::eval() const {
    auto sm_p = Parameters::GetInstance(0); // SM params
    auto flav_p = Parameters::GetInstance(3); // Flavor params
    
    double m_B = (*flav_p)("FMASS", 521);
    double life_B = (*flav_p)("FLIFE", 521);
    double f_B = (*flav_p)("FCONST", 52101);
    double m_tau = (*sm_p)("MASS", 15);
    double V_ub = std::abs(Parameters::get_c_CKM_entry(02)); 
    double G_F = (*sm_p)("SMINPUTS", 2);
    
    double BR_SM = std::pow(G_F * f_B * V_ub * m_tau * (1 - std::pow(m_tau / m_B, 2)), 2) * life_B * m_B / (8 * M_PI * HBAR);
    
    double np_fact = 1;
    if (model == Model::SUSY) {
        double m_Hp = (*Parameters::GetInstance(1))("MASS", 37);
        double tan_b = (*Parameters::GetInstance(1))("HMIX", 2);
        double eps_0 = EpsilonCalculator::GetInstance()->epsilon_0();
        np_fact = std::pow(1 - std::pow(m_B * tan_b / m_Hp, 2) / (1 + eps_0 * tan_b), 2);
    } else if (model == Model::THDM) {
        double m_Hp = (*Parameters::GetInstance(2))("MASS", 37);
        double tan_b = (*Parameters::GetInstance(2))("HMIX", 2);
        double l_bb = (*Parameters::GetInstance(2))("YD", 33);
        double l_tt = (*Parameters::GetInstance(2))("YU", 33);
        np_fact = std::pow(1 - std::pow(m_B / m_Hp, 2) * l_bb * l_tt, 2);
    }

    return np_fact * BR_SM;
}