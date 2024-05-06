#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath> // Pour std::abs
#include <memory>
#include "Wilson.h"
#include "MemoryManager.h"
#include "Logger.h"
#include "Wilson_susy.h"

void writeCoefficientsToFile(const std::string& strat_name, const std::string& fileName, const std::shared_ptr<InitializationStrategy>& strategy, double Q_match) {
    std::ofstream file(fileName);

    file << "Q,alpha_s";
    for (int i = 1; i <= 10; ++i) {
        file << ",C" << i << "_real,C" << i << "_imag";
    }
    file << "\n";

    MemoryManager::GetInstance("Test/testInput.slha", {0, 1})->init();
    Parameters* sm = Parameters::GetInstance();
    WilsonManager* wm = WilsonManager::GetInstance(strat_name, 81.0, strategy);


    double alpha_s = (*sm).QCDRunner.runningAlphasCalculation(Q_match);

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

    auto loStrategy = std::make_shared<SUSY_LO_Strategy>();
    auto nloStrategy = std::make_shared<SUSY_NLO_Strategy>();
    auto nnloStrategy = std::make_shared<SUSY_NNLO_Strategy>();

    
    writeCoefficientsToFile("NLO", "../csv/susy/WilsonCoefficients_NLO.csv", nloStrategy, 81);
    writeCoefficientsToFile("LO", "../csv/susy/WilsonCoefficients_LO.csv", loStrategy, 81);
    // writeCoefficientsToFile("NNLO", "../csv/susy/WilsonCoefficients_NNLO.csv", nnloStrategy, 81);
    
    WilsonManager::Cleanup();
    return 0;
}