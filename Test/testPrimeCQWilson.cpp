#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath> // Pour std::abs
#include <memory>
#include "Wilson.h"
#include "MemoryManager.h"
#include "Logger.h"
#include "Wilson_susy.h"
#include "config.hpp"

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

    double answer = 42.;
    wm->setScale(answer);

    double alpha_s = (*sm).QCDRunner.runningAlphasCalculation(answer);

    file << answer << "," << alpha_s;

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
    auto nloStrategy = std::make_shared<SM_NLO_Strategy>();
    auto nnloStrategy = std::make_shared<SM_NNLO_Strategy>();

    std::string root_file = project_root.data();
    
    writeCoefficientsToFile("NLO", root_file + "/Test/csv/sm/WilsonCoefficients_PRIMECQ_NLO.csv", nloStrategy, 81);
    writeCoefficientsToFile("LO", root_file + "/Test/csv/sm/WilsonCoefficients_PRIMECQ_LO.csv", loStrategy, 81);
    writeCoefficientsToFile("NNLO", root_file + "/Test/csv/sm/WilsonCoefficients_PRIMECQ_NNLO.csv", nnloStrategy, 81);
    
    WilsonManager::Cleanup();
    return 0;
}