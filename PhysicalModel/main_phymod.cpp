#include "Wilsonv2.h"
#include "Wilson_susyv2.h"
#include "Wilson_THDMv2.h"
#include "WilsonManager.h"

#include <iostream>
int main() {
    auto mm = MemoryManager::GetInstance("Test/InputFiles/testinput_thdm.lha", {0,2});  // Initialize program manager with LHA file containing SMINPUTS block
    mm->init();

    auto* manager = CoefficientManager::GetInstance("StandardModel");
    manager->registerCoefficientGroup("BCoefficients", std::make_unique<BCoefficientGroup>(81.0));
    manager->registerCoefficientGroup("ScalarCoefficients", std::make_unique<BScalarCoefficientGroup>(81.0));
    manager->setGroupScale("BCoefficients", 42.);
    manager->setQMatch();
    manager->setMatchingCoefficient("BCoefficients", "LO");
    manager->getRunCoefficient("BCoefficients", "C1", "LO");
    manager->setGroupScale("BCoefficients", 100.0);
    manager->printGroupCoefficients("BCoefficients");
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