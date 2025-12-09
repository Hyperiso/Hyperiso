#include "MemoryManager.h"
#include "Parameters.h"
#include "Logger.h"
#include "Include.h"
#include <iostream>
#include <cassert>

int main() {

    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    auto mm = MemoryManager::GetInstance(); 
    mm->init("Test/InputFiles/testInput.flha", Model::SM, false, true); 

    auto wc = Parameters::GetInstance(ParameterType::WILSON); 
    double scale = (*wc)("REWCOEF", -1);
    double C7_NLO = (*wc)("REWCOEF", (int)WCoef::C7 * 10 + (int)QCDOrder::NLO-1);

    std::cout << scale << std::endl;
    std::cout << C7_NLO << std::endl;


}