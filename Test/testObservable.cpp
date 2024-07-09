#include "MemoryManager.h"
#include "Parameters.h"
#include "Logger.h"
#include "Observable.h"
#include "Observables.h"

#include <iostream>
#include <cassert>

int main() {

    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    auto mm = MemoryManager::GetInstance("Test/testInput.flha", {0, 3});  // Initialize program manager with LHA file containing SMINPUTS block
    mm->init();  // Initialize parameters from given LHA file

    auto flavp = Parameters::GetInstance(3);
    double m_Bs = (*flavp)("MASS", 531);
    double m_Bd = (*flavp)("MASS", 511);

    auto sm_p = Parameters::GetInstance(0);

    Observable bs_mumu(Observables::BR_BS_MUMU, m_Bs, 2, 0); 
    Observable bs_mumu_untag(Observables::BR_BS_MUMU_UNTAG, m_Bs, 2, 0); 
    Observable bd_mumu(Observables::BR_BD_MUMU, m_Bd, 2, 0);
    Observable bu_taunu(Observables::BR_BU_TAUNU, m_Bd, 2, 0);
    Observable bu_taunu_np(Observables::BR_BU_TAUNU_NP_ONLY, m_Bd, 2, 0);
    Observable delta_0(Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA, sm_p->QCDRunner.mb_1S() / 2, 1, 0, 2);

    LOG_INFO("---------- Observables -----------");
    LOG_INFO("BR(Bs > mu mu) = ", bs_mumu.getValue());
    LOG_INFO("Untagged BR(Bs > mu mu) = ", bs_mumu_untag.getValue());
    LOG_INFO("BR(Bd > mu mu) = ", bd_mumu.getValue());
    LOG_INFO("BR(Bu > tau nu) = ", bu_taunu.getValue());
    LOG_INFO("NP Only BR(Bd > tau mu) = ", bu_taunu_np.getValue());
    LOG_INFO("Delta_0(B > K* gamma) = ", delta_0.getValue());
}