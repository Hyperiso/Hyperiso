#include "ObsEvaluator.h"
#include "Wilson.h"
// #include "Wilson_susy.h"
// #include "Wilson_THDM.h"
#include "Logger.h"
#include "Math.h"
#include "Parameters.h"


WilsonManager *ObsEvaluator::computeWilsons(int model, int order, double scale) {
    double m_W = (*Parameters::GetInstance(0))("MASS", 24);
    WilsonManager* wm;
    switch (model) {
        case 0:
            switch (order) {
                case 0:
                    wm = WilsonManager::GetInstance("LO", m_W, std::make_shared<SM_LO_Strategy>());
                    break;
                case 1:
                    wm = WilsonManager::GetInstance("NLO", m_W, std::make_shared<SM_NLO_Strategy>());
                    break;
                case 2:
                    wm = WilsonManager::GetInstance("NNLO", m_W, std::make_shared<SM_NNLO_Strategy>());
                    break;
                default:
                    Logger::getInstance()->warn("Order too high required for SM Wilson coefficients, defaulting to order 2.");
                    wm = WilsonManager::GetInstance("NNLO", m_W, std::make_shared<SM_NNLO_Strategy>());
            }
            break;
        // case 1:
        //     switch (order) {
        //         case 0:
        //             wm = WilsonManager::GetInstance("LO", m_W, std::make_shared<SUSY_LO_Strategy>());
        //             break;
        //         case 1:
        //             wm = WilsonManager::GetInstance("NLO", m_W, std::make_shared<SUSY_NLO_Strategy>());
        //             break;
        //         case 2:
        //             wm = WilsonManager::GetInstance("NNLO", m_W, std::make_shared<SUSY_NNLO_Strategy>());
        //             break;
        //         default:
        //             Logger::getInstance()->warn("Order too high required for SUSY Wilson coefficients, defaulting to order 2.");
        //             wm = WilsonManager::GetInstance("NNLO", m_W, std::make_shared<SUSY_NNLO_Strategy>());
        //     }
        //     break;
        // case 2:
        //     switch (order) {
        //         case 0:
        //             wm = WilsonManager::GetInstance("LO", m_W, std::make_shared<THDM_LO_Strategy>());
        //             break;
        //         case 1:
        //             wm = WilsonManager::GetInstance("NLO", m_W, std::make_shared<THDM_NLO_Strategy>());
        //             break;
        //         case 2:
        //             wm = WilsonManager::GetInstance("NNLO", m_W, std::make_shared<THDM_NNLO_Strategy>());
        //             break;
        //         default:
        //             Logger::getInstance()->warn("Order too high required for THDM Wilson coefficients, defaulting to order 2.");
        //             wm = WilsonManager::GetInstance("NNLO", m_W, std::make_shared<THDM_NNLO_Strategy>());
        //     }
        //     break;
        default:
            Logger::getInstance()->error("Unknown model requested for Wilson coefficient calculation.");
            return nullptr;
    }
    wm->setScale(scale, true);
    return wm;
}

complex_t get_c_CKM_entry(int idx) {
    auto p = Parameters::GetInstance(0);
    return (*p)("RECKM", idx) + complex_t(0, 1) * (*p)("IMCKM", idx);
}

complex_t ObsEvaluator::Evaluate(Observable *o) {
    auto p = Parameters::GetInstance(0);
    WilsonManager* wm = ObsEvaluator::computeWilsons(o->getModel(), o->getOrder(), o->getScale());

    if (!wm) {
        return std::complex<double>(-1);
    }

    switch (o->getId()) {
        case Observables::BR_BS_MUMU:
            return ObsEvaluator::Bs_mumu(wm);
        case Observables::BR_BD_MUMU:
            return ObsEvaluator::Bd_mumu(wm);
        default:
            Logger::getInstance()->error("Unknown observable.");
            return std::complex<double>(-1);
    }
}

complex_t ObsEvaluator::Bs_mumu(WilsonManager* wm)
{
    auto sm_p = Parameters::GetInstance(0); // SM params
    auto flav_p = Parameters::GetInstance(2); // Flavor params

    complex_t C10 = wm->get_full(WilsonCoefficient::C10, 2);
    complex_t CP10 = wm->get_full(WilsonCoefficient::CP10, 0);
    complex_t CQ1 = wm->get_full(WilsonCoefficient::CQ1, 1);
    complex_t CQ2 = wm->get_full(WilsonCoefficient::CQ2, 1);
    complex_t CPQ1 = wm->get_full(WilsonCoefficient::CPQ1, 0);
    complex_t CPQ2 = wm->get_full(WilsonCoefficient::CPQ2, 0);

    auto logger = Logger::getInstance();
    logger->info("C10 at all order is : " + std::to_string(std::real(C10)));
    logger->info("CP10 at all order is : " + std::to_string(std::real(CP10)));
    logger->info("CQ1 at all order is : " + std::to_string(std::real(CQ1)));
    logger->info("CQ2 at all order is : " + std::to_string(std::real(CQ2)));
    logger->info("CPQ1 at all order is : " + std::to_string(std::real(CPQ1)));
    logger->info("CPQ2 at all order is : " + std::to_string(std::real(CPQ2)));

    double G_F = (*sm_p)("SMINPUT", 2);
    double inv_alpha_em = (*sm_p)("SMINPUT", 1);
    inv_alpha_em = 137.;
    G_F = 1.166e-5;
    double V_tbV_ts = std::abs(get_c_CKM_entry(22) * std::conj(get_c_CKM_entry(21))); 

    double m_Bs = (*flav_p)("MASS", 531);
    double f_Bs = flav_p->getFlavorParam(FlavorParamType::DECAY_CONSTANT, "531|1");
    double life_Bs = flav_p->getFlavorParam(FlavorParamType::LIFETIME, "531");
    
    double r = (*sm_p)("MASS", 13) / m_Bs;  // m_mu / m_Bs
    double x = m_Bs / (sm_p->QCDRunner.get_mb_pole() + (*sm_p)("MASS", 3)); // m_Bs / (m_b_pole + m_s)

    return std::pow(G_F * f_Bs * V_tbV_ts / inv_alpha_em, 2) / (64 * HBAR) * std::pow(m_Bs, 3) * INV_PI3 * life_Bs * std::sqrt(1 - 4 * r * r) 
            * ((1 - 4 * r * r) * pow(x * std::abs(CQ1 - CPQ1), 2) + pow(std::abs(x * (CQ2 - CPQ2) + 2 * r * (C10 - CP10)), 2));
}

complex_t ObsEvaluator::Bd_mumu(WilsonManager* wm) {
    auto sm_p = Parameters::GetInstance(0); // SM params
    auto flav_p = Parameters::GetInstance(2); // Flavor params

    complex_t C10 = wm->get_full(WilsonCoefficient::C10, 2);
    complex_t CQ1 = wm->get_full(WilsonCoefficient::CQ1, 1);
    complex_t CQ2 = wm->get_full(WilsonCoefficient::CQ2, 1);

    double G_F = (*sm_p)("SMINPUT", 2);
    double inv_alpha_em = (*sm_p)("SMINPUT", 1);
    double V_tbV_td = std::abs(get_c_CKM_entry(33) * std::conj(get_c_CKM_entry(31))); 
    double m_Bd = (*flav_p)("FMASS", 511);
    double f_Bd = flav_p->getFlavorParam(FlavorParamType::DECAY_CONSTANT, "511|1");
    double life_Bd = flav_p->getFlavorParam(FlavorParamType::LIFETIME, "511");

    double r = (*sm_p)("MASS", 13) / m_Bd;  // m_mu / m_Bd
    double x = m_Bd / ((*sm_p)("SMINPUT", 5) + (*sm_p)("MASS", 2)); // m_Bd / (m_b_pole + m_d)

    return std::pow(G_F * f_Bd * V_tbV_td / inv_alpha_em, 2) / (64 * HBAR) * std::pow(m_Bd, 3) * INV_PI3 * life_Bd * std::sqrt(1 - 4 * r * r) 
        * ((1 - 4 * r * r) * pow(x * std::abs(CQ1), 2) + pow(std::abs(x * CQ2 + 2 * r * C10), 2));
}