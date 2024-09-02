#include "Wilsonv2.h"
#include "Wilson_susyv2.h"
#include "Wilson_THDMv2.h"
#include "WilsonManager.h"

#include <iostream>
int main() {
    auto mm = MemoryManager::GetInstance("Test/InputFiles/testinput_thdm.lha", {0,2});  // Initialize program manager with LHA file containing SMINPUTS block
    mm->init();

    // auto* manager = CoefficientManager::GetInstance("StandardModel");
    std::map<std::string, std::shared_ptr<CoefficientGroup>> groups = {std::make_pair("BCoefficient", std::make_shared<BCoefficientGroup>(81.)), std::make_pair("BScalarCoefficient", std::make_shared<BScalarCoefficientGroup>(81.)), std::make_pair("BPrimeCoefficient", std::make_shared<BPrimeCoefficientGroup>(81.))};

    auto* manager = CoefficientManager::Builder("StandardModel", groups , 81., 42., "LO");

    // manager->registerCoefficientGroup("BCoefficients", std::make_unique<BCoefficientGroup>(81.0));
    // manager->registerCoefficientGroup("ScalarBCoefficients", std::make_unique<BScalarCoefficientGroup>(81.0));
    // manager->setQMatch("BCoefficients", 81.);
    // manager->setMatchingCoefficient("BCoefficients", "NLO");

    // manager->setGroupScale("BCoefficients", 42.);
    // manager->setRunCoefficient("BCoefficients", "LO");

    // manager->setQMatch("ScalarBCoefficients", 81.);
    // manager->setMatchingCoefficient("ScalarBCoefficients", "NLO");

    // manager->setGroupScale("ScalarBCoefficients", 42.);
    // manager->setRunCoefficient("ScalarBCoefficients", "LO");
    
    // std::cout << "matching coeff LO : " <<  manager->getMatchingCoefficient("BCoefficients", "C7", "LO") << std::endl;
    // std::cout << "matching coeff NLO : " <<  manager->getMatchingCoefficient("BCoefficients", "C7", "NLO") << std::endl;
    // std::cout << "matching coeff full NLO : " <<  manager->getFullMatchingCoefficient("BCoefficients", "C7", "NLO") << std::endl;
    // std::cout << "run coeff : " <<  manager->getRunCoefficient("BCoefficients", "C7", "LO") << std::endl;

    // std::cout << "scalar matching coeff LO : " <<  manager->getMatchingCoefficient("ScalarBCoefficients", "CQ1", "LO") << std::endl;
    // std::cout << "scalar matching coeff NLO : " <<  manager->getMatchingCoefficient("ScalarBCoefficients", "CQ1", "NLO") << std::endl;
    // std::cout << "scalar matching coeff full NLO : " <<  manager->getFullMatchingCoefficient("ScalarBCoefficients", "CQ1", "NLO") << std::endl;
    // std::cout << "scalar run coeff : " <<  manager->getRunCoefficient("ScalarBCoefficients", "CQ1", "LO") << std::endl;
    // manager->setGroupScale("BCoefficients", 100.0);
    // manager->printGroupCoefficients("BCoefficients");

    // C4 C4_test{81.};
    // C4_test.LO_calculation();
    // C4_test.NLO_calculation();
    // C4_test.NNLO_calculation();

    // BCoefficientGroup bcoeff{81.};
    // bcoeff.set_Q_run(81.);
    // bcoeff.init_LO();
    // bcoeff.set_base_1_LO();

    // bcoeff.init_NLO();
    // bcoeff.set_base_1_NLO();

    // bcoeff.init_NNLO();
    // bcoeff.set_base_1_NNLO();

    // std::cout <<  bcoeff << std::endl;

    // std::cout << ".................................................................................................." << std::endl;
    
    // BCoefficientGroup_THDM bcoeff_thdm{81.};

    // // for (auto& coeff : bcoeff_thdm) {
    // //     std::cout << coeff.first << std::endl;
    // //     std:: cout << coeff.second << std::endl;

    // // }
    // bcoeff_thdm.set_Q_run(81.);
    // bcoeff_thdm.init_LO();
    // bcoeff_thdm.set_base_1_LO();

    // bcoeff_thdm.init_NLO();
    // bcoeff_thdm.set_base_1_NLO();

    // bcoeff_thdm.init_NNLO();
    // bcoeff_thdm.set_base_1_NNLO();

    // std::cout <<  bcoeff_thdm << std::endl;
    return 0;
}