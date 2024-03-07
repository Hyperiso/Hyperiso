#include "ObsEvaluator.h"
#include "Wilson.h"
#include "QCDParameters.h"
#include "Logger.h"
#include "Math.h"
#include "Parameters.h"

complex_t ObsEvaluator::Evaluate(Observable *o)
{
    WilsonManager::GetInstance()->setScale(o->getScale());
    Parameters::GetInstance()->setScale(o->getScale());

    switch (o->getId())
    {
    case Observables::BR_BS_MUMU:
        return ObsEvaluator::Bs_mumu();
    case Observables::BR_BD_MUMU:
        return ObsEvaluator::Bd_mumu();
    default:
        Logger::getInstance()->error("Unknown observable.");
        return std::complex<double>(-1);
    }
}


complex_t ObsEvaluator::Bs_mumu()
{
    auto wm = WilsonManager::GetInstance();
    auto p = Parameters::GetInstance();
    complex_t C10 = wm->get(WilsonCoefficient::C10, 2);
    complex_t CP10 = wm->get(WilsonCoefficient::CP10, 0);
    complex_t CQ1 = wm->get(WilsonCoefficient::CQ1, 1);
    complex_t CQ2 = wm->get(WilsonCoefficient::CQ2, 1);
    complex_t CPQ1 = wm->get(WilsonCoefficient::CPQ1, 0);
    complex_t CPQ2 = wm->get(WilsonCoefficient::CPQ2, 0);

    double r = (*p)("MASS", 13) / (*p)("MASS", 531);  // m_mu / m_Bs
    double x = (*p)("MASS", 531) / ((*p)("SMINPUT", 5) + (*p)("MASS", 3)); // m_Bs / (m_b_pole + m_s)

    return std::pow((*p)("SMINPUT", 2 /* G_F */) * (*p)("FCONST", 531 /* f_Bs */) * std::abs((*p)("VCKM", 33) * std::conj((*p)("VCKM", 32)) /* V_tb.V_ts* */) / (*p)("SMINPUT", 1 /* inv_alpha_em */), 2) / (64 * HBAR) 
        * std::pow((*p)("MASS", 531 /* m_Bs */), 3) * INV_PI3 * (*p)("FLIFE", 531 /* tau_Bs */) * std::sqrt(1 - 4 * r * r) 
        * ((1 - 4 * r * r) * pow(x * std::abs(CQ1 - CPQ1), 2) + pow(std::abs(x * (CQ2 - CPQ2) + 2 * r * (C10 - CP10)), 2));
}


complex_t ObsEvaluator::Bd_mumu()
{
    auto wm = WilsonManager::GetInstance();
    auto p = Parameters::GetInstance();
    complex_t C10 = wm->get(WilsonCoefficient::C10, 2);
    complex_t CQ1 = wm->get(WilsonCoefficient::CQ1, 1);
    complex_t CQ2 = wm->get(WilsonCoefficient::CQ2, 1);

    double r = (*p)("MASS", 13) / (*p)("MASS", 511);  // m_mu / m_Bd
    double x = (*p)("MASS", 511) / ((*p)("SMINPUT", 5) + (*p)("MASS", 2)); // m_Bd / (m_b_pole + m_d)

    return std::pow((*p)("SMINPUT", 2 /* G_F */) * (*p)("FCONST", 511 /* f_Bd */) * std::abs((*p)("VCKM", 33) * std::conj((*p)("VCKM", 31)) /* V_tb.V_td* */) / (*p)("SMINPUT", 1 /* inv_alpha_em */), 2) / (64 * HBAR) 
        * std::pow((*p)("MASS", 511 /* m_Bd */), 3) * INV_PI3 * (*p)("FLIFE", 511 /* tau_Bd */) * std::sqrt(1 - 4 * r * r) 
        * ((1 - 4 * r * r) * pow(x * std::abs(CQ1), 2) + pow(std::abs(x * CQ2 + 2 * r * C10), 2));
}