#include "MemoryManager.h"
#include "Parameters.h"
#include "Logger.h"
#include "config.hpp"
#include <iostream>
#include <cassert>

int main() {

    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    std::string root_file = project_root.data();

    auto mm = MemoryManager::GetInstance();  // Initialize program manager with LHA file containing SMINPUTS block
    mm->init(root_file + "Test/InputFiles/testInput.slha", Model::SUSY);  // Initialize parameters from given LHA file

    auto sm_params = Parameters::GetInstance(ParameterType::SUSY); 
    double mass_n1 = (*sm_params)("MASS", 1000022);
    std::cout << mass_n1 << std::endl;
    assert(mass_n1 == 105.39);  // gauge[3] is g_s

}