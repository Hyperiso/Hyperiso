#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath> // Pour std::abs
#include <memory>
#include "Wilson.h"
#include "MemoryManager.h"

void writeCoefficientsToFile(const std::string& strat_name, const std::string& fileName, const std::shared_ptr<InitializationStrategy>& strategy) {
    std::ofstream file(fileName);

    file << "Q,alpha_s";
    for (int i = 1; i <= 9; ++i) {
        file << ",C" << i << "_real,C" << i << "_imag";
    }
    file << "\n";

    double Q_initial = 1.0;
    double Q_final = 100.0;
    double Q_step = 1.0;

    MemoryManager::GetInstance()->init();

    WilsonManager* wm = WilsonManager::GetInstance(strat_name, 81.0, strategy);
    Parameters* sm = Parameters::GetInstance();
    for (double Q = Q_initial; Q <= Q_final; Q += Q_step) {
        wm->setScale(Q);

        double alpha_s = (*sm).QCDRunner.runningAlphasCalculation(Q);

        file << Q << "," << alpha_s;

                for (int i = 1; i <= 9; ++i) {
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
    }

    file.close();

    std::cout << "Done writing " << fileName << std::endl;
}

int main() {

    auto loStrategy = std::make_shared<SM_LO_Strategy>();
    auto nloStrategy = std::make_shared<SM_NLO_Strategy>();
    auto nnloStrategy = std::make_shared<SM_NNLO_Strategy>();

    
    writeCoefficientsToFile("NLO", "WilsonCoefficients_NLO.csv", nloStrategy);
    writeCoefficientsToFile("LO", "WilsonCoefficients_LO.csv", loStrategy);
    writeCoefficientsToFile("NNLO", "WilsonCoefficients_NNLO.csv", nnloStrategy);

    WilsonManager::Cleanup();
    return 0;
}
