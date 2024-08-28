#include "Wilsonv2.h"
#include "Wilson_susyv2.h"
#include "Wilson_THDMv2.h"
#include <iostream>
int main() {
    MemoryManager::GetInstance()->init();

    C4 C4_test{81.};
    C4_test.LO_calculation();
    C4_test.NLO_calculation();
    C4_test.NNLO_calculation();

    BCoefficientGroup bcoeff{81.};
    bcoeff.set_Q_run(81.);
    bcoeff.init_LO();
    bcoeff.set_base_1_LO();

    bcoeff.init_NLO();
    bcoeff.set_base_1_NLO();

    bcoeff.init_NNLO();
    bcoeff.set_base_1_NNLO();

    std::cout <<  bcoeff << std::endl;

    std::cout << ".................................................................................................." << std::endl;
    return 0;
}