#include "MemoryManager.h"
#include "Parameters.h"
#include "Logger.h"

#include <iostream>
#include <cassert>

int main() {

    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    auto mm = MemoryManager::GetInstance("Test/testInput.slha", {0, 1});  // Initialize program manager with LHA file containing SMINPUTS block
    mm->init();  // Initialize parameters from given LHA file

    auto sm_params = Parameters::GetInstance(0); 
    double alpha_s_MZ = std::pow((*sm_params)("GAUGE", 3), 2) / (4 * M_PI);
    assert(std::abs(alpha_s_MZ - 0.1172) < 1e-5);  // gauge[3] is g_s

    auto susy_params = Parameters::GetInstance(1);
    double tan_beta = (*susy_params)("HMIX", 2);
    assert(std::abs(tan_beta - 10) < 1e-5);
}