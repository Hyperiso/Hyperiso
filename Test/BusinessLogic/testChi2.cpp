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

    auto mm = MemoryManager::GetInstance();  // Initialize program manager with LHA file containing SMINPUTS block
    mm->init("Test/testInput.flha", Model::SM);  // Initialize parameters from given LHA file

    auto flavp = Parameters::GetInstance(ParameterType::FLAVOR);
    double m_Bs = (*flavp)("FMASS", 531);
    double m_Bd = (*flavp)("FMASS", 511);
    double m_Bu = (*flavp)("FMASS", 521);

    std::shared_ptr<BR_Bs_mumu_untag>   br_Bs__mu_mu = std::make_shared<BR_Bs_mumu_untag>(Model::SM, QCDOrder::NNLO, m_Bs);
    std::shared_ptr<BR_Bd_mumu>         br_Bd__mu_mu = std::make_shared<BR_Bd_mumu>(Model::SM, QCDOrder::NNLO, m_Bd);
    std::shared_ptr<BR_Bu_taunu>        br_Bu__tau_nu = std::make_shared<BR_Bu_taunu>(Model::SM, QCDOrder::NNLO, m_Bu);
    // Delta0_B__Kstar_gamma delta_0 (0, 1, sm_p->get_QCD_masse("mb_1S") / 2);

    ModelEvaluator me {{
        br_Bs__mu_mu,
        br_Bd__mu_mu,
        br_Bu__tau_nu,
    }};

    LOG_INFO("---------- Chi2 calculation -----------");
    LOG_INFO("Chi2 = ", me.chi2());
}