#include "MemoryManager.h"
#include "Parameters.h"
#include "Logger.h"
#include <iostream>
#include <cassert>

int main() {

    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    auto mm = MemoryManager::GetInstance("Test/InputFiles/testInput.slha", {1});  // Initialize program manager with LHA file containing SMINPUTS block
    mm->init();  // Initialize parameters from given LHA file

    auto sm_params = Parameters::GetInstance(0); 
    double mass_n1 = (*sm_params)("MASS", 1000022);
    std::cout << mass_n1 << std::endl;
    assert(mass_n1 == 105.39);  // gauge[3] is g_s

}