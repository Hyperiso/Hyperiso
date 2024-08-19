#include "MemoryManager.h"
#include "Parameters.h"
#include "Logger.h"
#include <iostream>
#include <cassert>

int main() {

    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    auto mm = MemoryManager::GetInstance("Test/InputFiles/testInput.flha", {0});  // Initialize program manager with LHA file containing SMINPUTS block
    mm->init();  // Initialize parameters from given LHA file
    LOG_INFO("tout va bien");
    auto sm_params = Parameters::GetInstance(0); 
    double alpha_s_MZ = std::pow((*sm_params)("GAUGE", 3), 2) / (4 * M_PI);
    std::cout << alpha_s_MZ << std::endl;
    assert(std::abs(alpha_s_MZ - 0.11999) < 1e-5);  // gauge[3] is g_s

}