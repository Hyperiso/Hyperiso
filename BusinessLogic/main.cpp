#include "MemoryManager.h"
#include "Logger.h"
#include <functional>
#include <iostream>
#include <cassert>
#include "ObservableInterface.h"

int main() {

    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    auto mm = MemoryManager::GetInstance();  // Initialize program manager with LHA file containing SMINPUTS block
    mm->init("Test/InputFiles/testInput.flha", Model::SM);  // Initialize parameters from given LHA file

    auto interface = ObservableInterface();
    interface.add_observable(Observables::BR_BS_MUMU, QCDOrder::NNLO)
             .add_observable(Observables::BR_BU_TAU_NU, QCDOrder::LO);

    interface.add_observable_parameters(Observables::BR_BS_MUMU, {
        ParamId{ParameterType::SM, "MASS", 13},
        ParamId{ParameterType::FLAVOR, "FMASS", 531},
        ParamId{ParameterType::FLAVOR, "FCONST", 53101}
    });

    interface.add_observable_parameters(Observables::BR_BU_TAU_NU, {
        ParamId{ParameterType::SM, "MASS", 15},
        ParamId{ParameterType::FLAVOR, "FMASS", 521},
        ParamId{ParameterType::FLAVOR, "FCONST", 52101}
    });

    LOG_INFO(interface.compute_observable(Observables::BR_BS_MUMU), "+-", interface.compute_uncertainty(Observables::BR_BS_MUMU));
    LOG_INFO(interface.compute_observable(Observables::BR_BU_TAU_NU), "+-", interface.compute_uncertainty(Observables::BR_BU_TAU_NU));
    LOG_INFO(interface.compute_chi2());

}