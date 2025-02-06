#include "MemoryManager.h"
#include "Logger.h"
#include <functional>
#include <iostream>
#include <cassert>
#include "ObservableInterface.h"

int main() {

    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    auto mm = MemoryManager::GetInstance();  // Initialize program manager with LHA file containing SMINPUTS block
    mm->init("Test/InputFiles/testinput_thdm.lha", Model::THDM);  // Initialize parameters from given LHA file

    auto interface = ObservableInterface();
    interface.add_observable(Observables::BR_B_XS_GAMMA, QCDOrder::NNLO, false);
    LOG_INFO("Observable added to manager");

    LOG_INFO(interface.compute_observable(Observables::BR_B_XS_GAMMA)/*, "+-", interface.compute_uncertainty(Observables::BR_B_XS_GAMMA)*/);
    LOG_INFO(interface.compute_chi2());

}