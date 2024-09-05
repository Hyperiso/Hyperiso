#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath> // Pour std::abs
#include <memory>
#include "Wilson.h"
#include "MemoryManager.h"
#include "Logger.h"
#include "Wilson_THDM.h"
#include "config.hpp"

void writeCoefficientsToFile(const std::string& strat_name, const std::string& fileName, const std::shared_ptr<InitializationStrategy>& strategy, double Q_match) {
    std::ofstream file(fileName);

    file << "Q,alpha_s";
    for (int i = 1; i <= 10; ++i) {
        file << ",C" << i << "_real,C" << i << "_imag";
    }
    file << "\n";

    MemoryManager::GetInstance("Test/testinput_thdm.lha", {0, 2})->init();
    Parameters* sm = Parameters::GetInstance();
    WilsonManager* wm = WilsonManager::GetInstance(strat_name, 81.0, strategy);

    double alpha_s = (*sm).alpha_s(Q_match);

    file << Q_match << "," << alpha_s;

            for (int i = 0; i <= 9; ++i) {
        complex_t C = {0,0};

        if (strat_name == "LO")
            C = wm->get_matchs(static_cast<WilsonCoefficient>(i), 0);
        else if (strat_name == "NLO")
            C = wm->get_matchs(static_cast<WilsonCoefficient>(i), 1);
        else if (strat_name == "NNLO")
            C = wm->get_matchs(static_cast<WilsonCoefficient>(i), 2);
        file << "," << C.real() << "," << C.imag();
    }

    file << "\n"; 

    file.close();

}

int main() {

    Logger* logger = Logger::getInstance();
    logger->setLevel(Logger::LogLevel::INFO);

    auto loStrategy = std::make_shared<THDM_LO_Strategy>();
    auto nloStrategy = std::make_shared<THDM_NLO_Strategy>();
    auto nnloStrategy = std::make_shared<THDM_NNLO_Strategy>();

    std::string root_file = project_root.data();

    writeCoefficientsToFile("NLO", root_file + "/Test/csv/thdm/WilsonCoefficients_NLO.csv", nloStrategy, 81);
    writeCoefficientsToFile("LO", root_file + "/Test/csv/thdm/WilsonCoefficients_LO.csv", loStrategy, 81);
    writeCoefficientsToFile("NNLO", root_file + "/Test/csv/thdm/WilsonCoefficients_NNLO.csv", nnloStrategy, 81);
    
    WilsonManager::Cleanup();
    return 0;
}