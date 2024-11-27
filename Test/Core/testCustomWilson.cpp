#include "MemoryManager.h"
#include "Parameters.h"
#include "Logger.h"
#include "General.h"
#include <iostream>
#include <cassert>

int main() {

    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    auto mm = MemoryManager::GetInstance();  // Initialize program manager with LHA file containing SMINPUTS block
    mm->init("Test/InputFiles/testInput.flha", Model::SM, false, true);  // Initialize parameters from given LHA file

    auto wc = Parameters::GetInstance(ParameterType::WILSON); 
    double scale = (*wc)("REWCOEF", -1);
    double C7_NLO = (*wc)("REWCOEF", (int)BWilsonCoefficients::C7 * 100 + (int)QCDOrder::NLO * 10 + 2);

    std::cout << scale << std::endl;
    std::cout << C7_NLO << std::endl;
}