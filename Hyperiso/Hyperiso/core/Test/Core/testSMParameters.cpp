#include "MemoryManager.h"
#include "Parameters.h"
#include "Logger.h"
#include "config.hpp"
#include <iostream>
#include <cassert>

int main() {

    Logger::getInstance()->setLevel(Logger::LogLevel::DEBUG);

    std::string root_file = project_root.data();

    auto mm = MemoryManager::GetInstance();  // Initialize program manager with LHA file containing SMINPUTS block
    mm->init(root_file + "Test/InputFiles/testInput.flha", Model::SM);  // Initialize parameters from given LHA file

    auto sm_params = Parameters::GetInstance(ParameterType::SM); 
    double alpha_s_MZ = std::pow((*sm_params)("GAUGE", 3), 2) / (4 * M_PI);
    std::cout << alpha_s_MZ << std::endl;
    assert(std::abs(alpha_s_MZ - 0.1172) < 1e-5);  // gauge[3] is g_s
}