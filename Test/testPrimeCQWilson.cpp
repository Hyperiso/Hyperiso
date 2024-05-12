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

    std::vector<std::string> name {"C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "CQ1", "CQ2", "CP7", "CP8", "CP9", "CP10", "CPQ1", "CPQ2"};

    file << "Q,alpha_s";
    for (int i = 10; i <= 17; ++i) {
        file << "," << name[i] << "_real," << name[i] << "_imag";
    }
    file << "\n";

    MemoryManager::GetInstance("Test/testInput.slha", {0, 1})->init();
    Parameters* sm = Parameters::GetInstance();
    WilsonManager* wm = WilsonManager::GetInstance(strat_name, 81.0, strategy);

    wm->setScale(Q_match);

    double alpha_s = (*sm).QCDRunner.runningAlphasCalculation(Q_match);

    file << Q_match << "," << alpha_s;

            for (int i = 10; i <= 17; ++i) {
        complex_t C = {0,0};

        if (strat_name == "LO")
            C = wm->get(static_cast<WilsonCoefficient>(i), 0);
        else if (strat_name == "NLO")
            C = wm->get(static_cast<WilsonCoefficient>(i), 1);
        else if (strat_name == "NNLO")
            C = wm->get(static_cast<WilsonCoefficient>(i), 2);
        file << "," << C.real() << "," << C.imag();
    }

    file << "\n"; 

    file.close();

}

int main() {

    Logger* logger = Logger::getInstance();
    logger->setLevel(Logger::LogLevel::INFO);

    auto loStrategy = std::make_shared<SM_LO_Strategy>();
    auto susyloStrategy = std::make_shared<SUSY_LO_Strategy>();
    auto susynloStrategy = std::make_shared<SUSY_NLO_Strategy>();

    
    // writeCoefficientsToFile("NLO", "../csv/susy/WilsonCoefficients_NLO.csv", nloStrategy, 81);
    writeCoefficientsToFile("LO", "../csv/sm/WilsonCoefficients_PRIMECQ_LO.csv", loStrategy, 81);
    // writeCoefficientsToFile("NNLO", "../csv/susy/WilsonCoefficients_NNLO.csv", nnloStrategy, 81);
    
    WilsonManager::Cleanup();
    return 0;
}