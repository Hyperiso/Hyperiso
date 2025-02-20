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

    auto mm = MemoryManager::GetInstance("Test/testInput.flha", {0, 3});
    mm->init();

    auto flavp = Parameters::GetInstance(3);
    double m_Bs = (*flavp)("MASS", 531);
    double m_Bd = (*flavp)("MASS", 511);

    auto sm_p = Parameters::GetInstance(0);
    BR_Bs_mumu br_Bs__mu_mu(0, 2, m_Bs);
    BR_Bd_mumu br_Bd__mu_mu(0, 2, m_Bd);

    br_Bs__mu_mu.print_gradient(std::cout);

    std::vector<std::shared_ptr<Observable>> obss = {std::make_shared<BR_Bs_mumu>(br_Bs__mu_mu), std::make_shared<BR_Bd_mumu>(br_Bd__mu_mu)};
    ModelEvaluator me (obss);

    LOG_INFO("---------- Observables -----------");
    LOG_INFO("BR(Bs > mu mu) = ", br_Bs__mu_mu.eval(), " +- ", std::sqrt(br_Bs__mu_mu.variance()));
    LOG_INFO("BR(Bd > mu mu) = ", br_Bd__mu_mu.eval(), " +- ", std::sqrt(br_Bs__mu_mu.variance()));
    LOG_INFO("\nChi squared = ", me.chi2());
}