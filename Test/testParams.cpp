#include "MemoryManager.h"
#include "Parameters.h"

#include <iostream>
#include <cassert>

int main() {
    auto mm = MemoryManager::GetInstance("../testInput.slha", {0, 1});  // Initialize program manager with LHA file containing SMINPUTS block
    mm->init();  // Initialize parameters from given LHA file

    auto sm_params = Parameters::GetInstance(0); 
    double alpha_s_MZ = std::pow((*sm_params)("GAUGE", 3), 2) / (4 * M_PI);
    assert(alpha_s_MZ == 0.12);  // gauge[3] is g_s

    auto susy_params = Parameters::GetInstance(1);
    double tan_beta = (*susy_params)("HMIX", 2);
    std::cout << tan_beta << std::endl;
    assert(tan_beta == 9.75139550e+00);
}