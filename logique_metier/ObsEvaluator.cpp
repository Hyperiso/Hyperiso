#include "ObsEvaluator.h"
#include "Wilson.h"
#include "QCDParameters.h"
#include "Logger.h"
#include "Math.h"
#include "Parameters.h"

complex_t ObsEvaluator::Evaluate(Observable *o)
{
    auto p = Parameters::GetInstance(0);
    WilsonManager::GetInstance("LO", (*p)("MASS", 24), std::make_shared<SM_LO_Strategy>())->setScale(o->getScale());

    switch (o->getId()) {
    case Observables::BR_BS_MUMU:
        return ObsEvaluator::Bs_mumu();
    case Observables::BR_BD_MUMU:
        return ObsEvaluator::Bd_mumu();
    default:
        Logger::getInstance()->error("Unknown observable.");
        return std::complex<double>(-1);
    }
}


complex_t ObsEvaluator::Bs_mumu() {
    auto sm_p = Parameters::GetInstance(0); // SM params
    auto flav_p = Parameters::GetInstance(2); // Flavor params

    auto wm = WilsonManager::GetInstance("LO", (*sm_p)("MASS", 24), std::make_shared<SM_LO_Strategy>());
    complex_t C10 = wm->get(WilsonCoefficient::C10, 2);
    complex_t CP10 = wm->get(WilsonCoefficient::CP10, 0);
    complex_t CQ1 = wm->get(WilsonCoefficient::CQ1, 1);
    complex_t CQ2 = wm->get(WilsonCoefficient::CQ2, 1);
    complex_t CPQ1 = wm->get(WilsonCoefficient::CPQ1, 0);
    complex_t CPQ2 = wm->get(WilsonCoefficient::CPQ2, 0);

    double G_F = (*sm_p)("SMINPUT", 2);
    double inv_alpha_em = (*sm_p)("SMINPUT", 1);
    double V_tbV_ts = std::abs((*sm_p)("VCKM", 33) * std::conj((*sm_p)("VCKM", 32))); 
    double m_Bs = (*flav_p)("FMASS", 531);
    double f_Bs = (*flav_p)("FCONST", 531);
    double life_Bs = (*flav_p)("FLIFE", 531);

    double r = (*sm_p)("MASS", 13) / m_Bs;  // m_mu / m_Bs
    double x = m_Bs / ((*sm_p)("SMINPUT", 5) + (*sm_p)("MASS", 3)); // m_Bs / (m_b_pole + m_s)

    return std::pow(G_F * f_Bs * V_tbV_ts / inv_alpha_em, 2) / (64 * HBAR) * std::pow(m_Bs, 3) * INV_PI3 * life_Bs * std::sqrt(1 - 4 * r * r) 
            * ((1 - 4 * r * r) * pow(x * std::abs(CQ1 - CPQ1), 2) + pow(std::abs(x * (CQ2 - CPQ2) + 2 * r * (C10 - CP10)), 2));
}


complex_t ObsEvaluator::Bd_mumu() {
    auto sm_p = Parameters::GetInstance(0); // SM params
    auto flav_p = Parameters::GetInstance(2); // Flavor params

    auto wm = WilsonManager::GetInstance("LO", (*sm_p)("MASS", 24), std::make_shared<SM_LO_Strategy>());
    complex_t C10 = wm->get(WilsonCoefficient::C10, 2);
    complex_t CQ1 = wm->get(WilsonCoefficient::CQ1, 1);
    complex_t CQ2 = wm->get(WilsonCoefficient::CQ2, 1);

    double G_F = (*sm_p)("SMINPUT", 2);
    double inv_alpha_em = (*sm_p)("SMINPUT", 1);
    double V_tbV_td = std::abs((*sm_p)("VCKM", 33) * std::conj((*sm_p)("VCKM", 31))); 
    double m_Bd = (*flav_p)("FMASS", 511);
    double f_Bd = (*flav_p)("FCONST", 511);
    double life_Bd = (*flav_p)("FLIFE", 511);

    double r = (*sm_p)("MASS", 13) / m_Bd;  // m_mu / m_Bd
    double x = m_Bd / ((*sm_p)("SMINPUT", 5) + (*sm_p)("MASS", 2)); // m_Bd / (m_b_pole + m_d)

    return std::pow(G_F * f_Bd * V_tbV_td / inv_alpha_em, 2) / (64 * HBAR) * std::pow(m_Bd, 3) * INV_PI3 * life_Bd * std::sqrt(1 - 4 * r * r) 
        * ((1 - 4 * r * r) * pow(x * std::abs(CQ1), 2) + pow(std::abs(x * CQ2 + 2 * r * C10), 2));
}