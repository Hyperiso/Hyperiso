#include "ObsEvaluator.h"
#include "./Physical_Model/Wilson.h"
#include "./Physical_Model/QCDParameters.h"

#define pi    3.1415926535897932385
#define hbar  6.58211889e-25 /* in GeV.s */

complex_t ObsEvaluator::Evaluate(Observable *o)
{
    auto wm = WilsonManager::GetInstance();
    wm->setScale(o->getScale());
    QCDParameters QCDParams;
    ObsEvaluator::alpha_s = QCDParams.runningAlphasCalculation(o->getScale());

    switch (o->getId())
    {
    case Observables::BR_Bs_mumu:
        return ObsEvaluator::Bs_mumu();
    case Observables::BR_Bd_mumu:
        return ObsEvaluator::Bd_mumu();
    default:
        // Log an error, unknown observable
    }
}


complex_t ObsEvaluator::Bs_mumu()
{
    auto wm = WilsonManager::GetInstance();
    double m_mu, G_F, inv_alpha_em, m_Bs, f_Bs, life_Bs, m_b_pole, m_s;
    complex_t Vtb, Vts;
    complex_t C10 = wm->get(WilsonCoefficient::C10, 2);
    complex_t CP10 = wm->get(WilsonCoefficient::CP10, 0);
    complex_t CQ1 = wm->get(WilsonCoefficient::CQ1, 1);
    complex_t CQ2 = wm->get(WilsonCoefficient::CQ2, 1);
    complex_t CPQ1 = wm->get(WilsonCoefficient::CPQ1, 0);
    complex_t CPQ2 = wm->get(WilsonCoefficient::CPQ2, 0);

    double r = m_mu / m_Bs;
    double x = m_Bs / (m_b_pole + m_s);

    return std::pow(G_F * f_Bs * std::abs(Vtb * std::conj(Vts)) / inv_alpha_em, 2) / (64 * hbar) * std::pow(m_Bs / pi, 3) * life_Bs 
        * std::sqrt(1 - 4 * r * r) * ((1 - 4 * r * r) * pow(x * std::abs(CQ1 - CPQ1), 2) + pow(std::abs(x * (CQ2 - CPQ2) + 2 * r * (C10 - CP10)), 2));
}


complex_t ObsEvaluator::Bd_mumu()
{
    auto wm = WilsonManager::GetInstance();
    double m_mu, G_F, inv_alpha_em, m_Bd, f_Bd, life_Bd, m_b_pole, m_d;
    complex_t Vtb, Vtd;
    complex_t C10 = wm->get(WilsonCoefficient::C10, 2);
    complex_t CQ1 = wm->get(WilsonCoefficient::CQ1, 1);
    complex_t CQ2 = wm->get(WilsonCoefficient::CQ2, 1);

    double r = m_mu / m_Bd;
    double x = m_Bd / (m_b_pole + m_d);

    return std::pow(G_F * f_Bd * std::abs(Vtb * std::conj(Vtd)) / inv_alpha_em, 2) / (64 * hbar) * std::pow(m_Bd / pi, 3) * life_Bd 
        * std::sqrt(1 - 4 * r * r) * ((1 - 4 * r * r) * pow(x * std::abs(CQ1), 2) + pow(std::abs(x * CQ2 + 2 * r * C10), 2));
}