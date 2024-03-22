#include "MemoryManager.h"
#include "Parameters.h"

#include <iostream>
#include <cassert>

int main() {
    auto mm = MemoryManager::GetInstance("../testInput.flha", {0});  // Initialize program manager with LHA file containing SMINPUTS block
    mm->init();  // Initialize parameters from given LHA file

    auto sm_params = Parameters::GetInstance(0); 
    double alpha_s_MZ = std::pow((*sm_params)("Coupling", 3), 2) / (4 * M_PI);
    assert(alpha_s_MZ == 0.12);  // gauge[3] is g_s
}