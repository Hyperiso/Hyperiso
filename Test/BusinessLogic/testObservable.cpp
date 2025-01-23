#include "MemoryManager.h"
#include "Parameters.h"
#include "Logger.h"
#include "ModelEvaluator.h"
#include "BKstarDecay.h"
#include "BllDecay.h"
#include "BXsDecay.h"
#include "Observable.h"

#include <iostream>
#include <cassert>

int main() {

    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    auto mm = MemoryManager::GetInstance();  // Initialize program manager with LHA file containing SMINPUTS block
    mm->init("Test/InputFiles/testinput_thdm.lha", Model::THDM);  // Initialize parameters from given LHA file

    auto sm_p = Parameters::GetInstance(ParameterType::SM);
    auto flavp = Parameters::GetInstance(ParameterType::FLAVOR);
    double m_Bs = (*flavp)("FMASS", 531);
    double m_Bd = (*flavp)("FMASS", 511);
    double m_Bu = (*flavp)("FMASS", 521);

    auto B_Ks = std::make_shared<BKstarDecay>(QCDOrder::NNLO, 80, QCDHelper::mass_b_1S() / 2);
    auto B_ll = std::make_shared<BllDecay>(QCDOrder::NNLO, 80, m_Bd);
    auto B_Xs = std::make_shared<BXsDecay>(QCDOrder::NNLO, 80, QCDHelper::mass_b_1S() / 2);

    Observable delta_0 (Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA, B_Ks);
    // Observable br_Bs__mu_mu(Observables::BR_BS_MUMU, B_ll);
    Observable br_Bs__mu_mu__untag(Observables::BR_BS_MUMU_UNTAG, B_ll);
    Observable br_Bd__mu_mu(Observables::BR_BD_MUMU, B_ll);
    Observable br_B__Xs_gamma(Observables::BR_B_XS_GAMMA, B_Xs);

    LOG_INFO("---------- Observables -----------");

    // LOG_INFO("CP-Averaged BR(Bs > mu mu) = ",  br_Bs__mu_mu.eval(), " +- ", std::sqrt(br_Bs__mu_mu.variance()));
    LOG_INFO("Untagged BR(Bs > mu mu) = ",  br_Bs__mu_mu__untag.eval(), " +- ", std::sqrt(br_Bs__mu_mu__untag.variance()));
    LOG_INFO("BR(Bd > mu mu) = ",           br_Bd__mu_mu.eval(),        " +- ", std::sqrt(br_Bd__mu_mu.variance()));
    // LOG_INFO("BR(Bu > tau nu) = ",          br_Bu__tau_nu.eval(),       " +- ", std::sqrt(br_Bu__tau_nu.variance()));
    LOG_INFO("Delta_0(B > K* gamma) = ",       delta_0.eval(),             " +- ", std::sqrt(delta_0.variance()));
    LOG_INFO("BR(B > X_s gamma) = ",       br_B__Xs_gamma.eval(),             " +- ", std::sqrt(br_B__Xs_gamma.variance()));

}