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

    interface.add_observable(Observables::BR_B__D_TAU_NU, QCDOrder::LO)
    .add_observable_parameters(Observables::BR_B__D_TAU_NU, {
        ParamId{ParameterType::SM, "RECKM", 12},
        ParamId{ParameterType::FF, "B_Dlnu", 1},
        ParamId{ParameterType::FF, "B_Dlnu", 2},
        ParamId{ParameterType::FF, "B_Dlnu", 3},
    });

    interface.add_observable(Observables::A_FB_B__D_TAU_NU, QCDOrder::LO)
    .add_observable_parameters(Observables::A_FB_B__D_TAU_NU, {
        ParamId{ParameterType::FF, "B_Dlnu", 2},
        ParamId{ParameterType::FF, "B_Dlnu", 3},
    });

    interface.add_observable(Observables::P_TAU_B__D_TAU_NU, QCDOrder::LO)
    .add_observable_parameters(Observables::P_TAU_B__D_TAU_NU, {
        ParamId{ParameterType::FF, "B_Dlnu", 2},
        ParamId{ParameterType::FF, "B_Dlnu", 3},
    });

    interface.add_observable(Observables::R_D, QCDOrder::LO)
    .add_observable_parameters(Observables::R_D, {
        ParamId{ParameterType::FF, "B_Dlnu", 2},
        ParamId{ParameterType::FF, "B_Dlnu", 3},
    });

    interface.add_observable(Observables::BR_B__DSTAR_TAU_NU, QCDOrder::LO)
    .add_observable_parameters(Observables::BR_B__DSTAR_TAU_NU, {
        ParamId{ParameterType::SM, "RECKM", 12},
        ParamId{ParameterType::FF, "B_Dslnu", 1},
        ParamId{ParameterType::FF, "B_Dslnu", 2},
        ParamId{ParameterType::FF, "B_Dslnu", 3},
        ParamId{ParameterType::FF, "B_Dslnu", 4},
    });

    interface.add_observable(Observables::A_FB_B__DSTAR_TAU_NU, QCDOrder::LO)
    .add_observable_parameters(Observables::A_FB_B__DSTAR_TAU_NU, {
        ParamId{ParameterType::FF, "B_Dslnu", 2},
        ParamId{ParameterType::FF, "B_Dslnu", 3},
        ParamId{ParameterType::FF, "B_Dslnu", 4},
    });

    interface.add_observable(Observables::P_TAU_B__DSTAR_TAU_NU, QCDOrder::LO)
    .add_observable_parameters(Observables::P_TAU_B__DSTAR_TAU_NU, {
        ParamId{ParameterType::FF, "B_Dslnu", 2},
        ParamId{ParameterType::FF, "B_Dslnu", 3},
        ParamId{ParameterType::FF, "B_Dslnu", 4},
    });

    interface.add_observable(Observables::P_D_B__DSTAR_TAU_NU, QCDOrder::LO)
    .add_observable_parameters(Observables::P_D_B__DSTAR_TAU_NU, {
        ParamId{ParameterType::FF, "B_Dslnu", 2},
        ParamId{ParameterType::FF, "B_Dslnu", 3},
        ParamId{ParameterType::FF, "B_Dslnu", 4},
    });

    interface.add_observable(Observables::R_DSTAR, QCDOrder::LO)
    .add_observable_parameters(Observables::R_DSTAR, {
        ParamId{ParameterType::FF, "B_Dslnu", 2},
        ParamId{ParameterType::FF, "B_Dslnu", 3},
        ParamId{ParameterType::FF, "B_Dslnu", 4},
    });

    LOG_INFO("BR(B > D tau nu)\t\t",    interface.compute_observable(Observables::BR_B__D_TAU_NU),    "\t+-", interface.compute_uncertainty(Observables::BR_B__D_TAU_NU));
    LOG_INFO("R(D)\t\t\t",              interface.compute_observable(Observables::R_D),               "\t+-", interface.compute_uncertainty(Observables::R_D));
    LOG_INFO("A_FB(B > D tau nu)\t",    interface.compute_observable(Observables::A_FB_B__D_TAU_NU),  "\t+-", interface.compute_uncertainty(Observables::A_FB_B__D_TAU_NU));
    LOG_INFO("P_tau(B > D tau nu)\t",   interface.compute_observable(Observables::P_TAU_B__D_TAU_NU), "\t+-", interface.compute_uncertainty(Observables::P_TAU_B__D_TAU_NU));

    LOG_INFO("BR(B > D* tau nu)\t",     interface.compute_observable(Observables::BR_B__DSTAR_TAU_NU),      "\t+-", interface.compute_uncertainty(Observables::BR_B__DSTAR_TAU_NU));
    LOG_INFO("R(D*)\t\t\t",             interface.compute_observable(Observables::R_DSTAR),                 "\t+-", interface.compute_uncertainty(Observables::R_DSTAR));
    LOG_INFO("A_FB(B > D* tau nu)\t",   interface.compute_observable(Observables::A_FB_B__DSTAR_TAU_NU),    "\t+-", interface.compute_uncertainty(Observables::A_FB_B__DSTAR_TAU_NU));
    LOG_INFO("P_tau(B > D* tau nu)\t",  interface.compute_observable(Observables::P_TAU_B__DSTAR_TAU_NU),   "\t+-", interface.compute_uncertainty(Observables::P_TAU_B__DSTAR_TAU_NU));
    LOG_INFO("P_D(B > D* tau nu)\t",    interface.compute_observable(Observables::P_D_B__DSTAR_TAU_NU),     "\t+-", interface.compute_uncertainty(Observables::P_D_B__DSTAR_TAU_NU));

}