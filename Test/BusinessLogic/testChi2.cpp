#include "MemoryManager.h"
#include "Parameters.h"
#include "Logger.h"
#include "ModelEvaluator.h"
#include "BKstarDecay.h"
#include "Observable.h"

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

    auto B_Ks = std::make_shared<BKstarDecay>(QCDOrder::NNLO, 80, Parameters::GetInstance(ParameterType::SM)->get_QCD_masse("mb_1S") / 2);

    // std::shared_ptr<BR_Bs_mumu_untag>   br_Bs__mu_mu = std::make_shared<BR_Bs_mumu_untag>(Model::SM, QCDOrder::NNLO, m_Bs);
    // std::shared_ptr<BR_Bd_mumu>         br_Bd__mu_mu = std::make_shared<BR_Bd_mumu>(Model::SM, QCDOrder::NNLO, m_Bd);
    // std::shared_ptr<BR_Bu_taunu>        br_Bu__tau_nu = std::make_shared<BR_Bu_taunu>(Model::SM, QCDOrder::NNLO, m_Bu);
    auto delta_0 = std::make_shared<Observable>(Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA, B_Ks);

    ModelEvaluator me {{
        delta_0
    }};

    LOG_INFO("---------- Chi2 calculation -----------");
    LOG_INFO("Chi2 = ", me.chi2());
}