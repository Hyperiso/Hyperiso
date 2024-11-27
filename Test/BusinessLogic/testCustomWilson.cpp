#include "MemoryManager.h"
#include "Parameters.h"
#include "Logger.h"
#include "General.h"
#include "Bs_mumu.h"
#include <iostream>
#include <cassert>

int main() {

    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    auto mm = MemoryManager::GetInstance();  // Initialize program manager with LHA file containing SMINPUTS block
    mm->init("Test/InputFiles/testInput.flha", Model::SM, false, true);  // Initialize parameters from given LHA file

    auto wc = Parameters::GetInstance(ParameterType::WILSON); 
    double scale = (*wc)("REWCOEF", -1);
    double C7_NLO = (*wc)("REWCOEF", (int)BWilsonCoefficients::C7 * 10 + (int)QCDOrder::NLO-1);

    std::cout << scale << std::endl;
    std::cout << C7_NLO << std::endl;

    BR_Bs_mumu br_Bs__mu_mu(Model::SM, QCDOrder::NNLO, 5); 
    LOG_INFO(br_Bs__mu_mu.eval());
}