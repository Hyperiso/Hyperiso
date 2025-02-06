#include "MemoryManager.h"
#include "Logger.h"
#include <functional>
#include <iostream>
#include <cassert>
#include "ObservableInterface.h"

int main() {

    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    auto mm = MemoryManager::GetInstance();  // Initialize program manager with LHA file containing SMINPUTS block
    mm->init("Test/InputFiles/testInput.lha", Model::THDM);  // Initialize parameters from given LHA file

    auto interface = ObservableInterface();
    interface.add_observable(Observables::BR_BU_TAU_NU, QCDOrder::LO, true);

    const size_t n = 100;
    double tanb_min = 1;
    double tanb_max = 100;
    double dtanb = (tanb_max - tanb_min) / n;
    double m_Hp = 50;
    double m_Hp_max = 1000;
    double dmHp = (m_Hp_max - m_Hp) / n;

    std::ofstream fs {"./BR_B_Xs_gamma.csv"};

    fs << "tan_beta,m_Hp,BR,uncert,chi2\n";

    for (size_t i = 0; i < n; i++) {
        double tanb = tanb_min;
        for (size_t j = 0; j < n; j++) {
            Parameters::GetInstance(ParameterType::THDM)->setBlockValue("YU", 22, 1 / tanb, true);
            Parameters::GetInstance(ParameterType::THDM)->setBlockValue("YD", 22, -tanb, true);
            Parameters::GetInstance(ParameterType::THDM)->setBlockValue("YL", 00, -tanb, true);
            Parameters::GetInstance(ParameterType::THDM)->setBlockValue("YL", 11, -tanb, true);
            Parameters::GetInstance(ParameterType::THDM)->setBlockValue("YL", 22, -tanb, true);
            Parameters::GetInstance(ParameterType::THDM)->setBlockValue("MASS", 37, m_Hp, true);
            Parameters::GetInstance(ParameterType::FLAVOR)->setBlockValue("FMASS", 511, i + j, true);
            
            fs << tanb << ","
               << m_Hp << ","
               << interface.compute_observable(Observables::BR_BU_TAU_NU) << ","
               << interface.compute_uncertainty(Observables::BR_BU_TAU_NU) << ","
               << interface.compute_chi2() << "\n"; 

            tanb += dtanb;
        }
        m_Hp += dmHp;
    }

    fs.close();

    // LOG_INFO(interface.compute_observable(Observables::BR_B_XS_GAMMA), "+-", interface.compute_uncertainty(Observables::BR_B_XS_GAMMA));
    // LOG_INFO(interface.compute_chi2());

}