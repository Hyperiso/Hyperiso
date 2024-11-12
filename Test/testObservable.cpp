#include "MemoryManager.h"
#include "Parameters.h"
#include "Logger.h"
#include "Bs_mumu.h"
#include "Bd_mumu.h"
#include "Bu_taunu.h"
#include "B__Kstar_gamma.h"
#include "ModelEvaluator.h"

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

    BR_Bs_mumu          br_Bs__mu_mu(0, 2, m_Bs);
    BR_Bs_mumu_untag    br_Bs__mu_mu__untag(0, 2, m_Bs);
    BR_Bd_mumu          br_Bd__mu_mu(0, 2, m_Bd);
    BR_Bu_taunu         br_Bu__tau_nu(0, 2, m_Bd);
    // Delta0_B__Kstar_gamma delta_0 (0, 1, sm_p->get_QCD_masse("mb_1S") / 2);

    // ModelEvaluator me ({br_Bs__mu_mu, br_Bd__mu_mu, br_Bu__tau_nu});

    LOG_INFO("---------- Observables -----------");
    LOG_INFO("BR(Bs > mu mu) = ", br_Bs__mu_mu.eval());
    LOG_INFO("Untagged BR(Bs > mu mu) = ", br_Bs__mu_mu__untag.eval());
    LOG_INFO("BR(Bd > mu mu) = ", br_Bd__mu_mu.eval());
    LOG_INFO("BR(Bu > tau nu) = ", br_Bu__tau_nu.eval());
    // LOG_INFO("Delta_0(B > K* gamma) = ", delta_0.eval());
    // LOG_INFO("\nChi squared = ", me.chi2());
}