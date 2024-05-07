#include "MemoryManager.h"
#include "Parameters.h"
#include "Logger.h"
#include "Observable.h"
#include "Observables.h"

#include <iostream>
#include <cassert>

int main() {

    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    auto mm = MemoryManager::GetInstance("Test/testInput.flha", {0, 2});  // Initialize program manager with LHA file containing SMINPUTS block
    mm->init();  // Initialize parameters from given LHA file

    auto flavp = Parameters::GetInstance(2);
    double m_Bs = (*flavp)("MASS", 531);
    double m_Bd = (*flavp)("MASS", 511);

    Observable bs_mumu(Observables::BR_BS_MUMU, m_Bs, 2, 0); 
    Observable bd_mumu(Observables::BR_BS_MUMU, m_Bd, 2, 0);

    std::cout << bs_mumu.getValue() << std::endl; 
    std::cout << bd_mumu.getValue() << std::endl; 
}