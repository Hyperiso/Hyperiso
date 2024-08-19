#include "WilsonUtils.h"
#include <fstream>
#include "MemoryManager.h"
#include "Logger.h"
#include "CompareCsv.h"

void writeCoefficientsToFile(const std::string& strat_name, const std::string& fileName, const std::shared_ptr<InitializationStrategy>& strategy, double Q_match, const std::string& model) {
    std::ofstream file(fileName);

    file << "Q,alpha_s";
    for (int i = 1; i <= 10; ++i) {
        file << ",C" << i << "_real,C" << i << "_imag";
    }
    file << "\n";

    if (model == "SM") {
        MemoryManager::GetInstance("Test/InputFiles/testinput_thdm.lha", {0})->init();
    }
    else if (model == "THDM") {
        MemoryManager::GetInstance("Test/InputFiles/testinput_thdm.lha", {0,2})->init();
    }
    else if (model == "SUSY") {
        MemoryManager::GetInstance("Test/InputFiles/testInput.slha", {0,1})->init();
    }
    else {
        LOG_ERROR("ModelError", "MODEL not known");
    }
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

    wm->Cleanup();
}

void writeCoefficientsPrimeCQToFile(const std::string& strat_name, const std::string& fileName, const std::shared_ptr<InitializationStrategy>& strategy, double Q_match, const std::string& model) {
    std::ofstream file(fileName);

    std::vector<std::string> name {"C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "CQ1", "CQ2", "CP1", "CP2", "CP3", "CP4", "CP5", "CP6", "CP7", "CP8", "CP9", "CP10", "CPQ1", "CPQ2"};

    file << "Q,alpha_s";
    for (int i = 10; i <= 23; ++i) {
        file << "," << name[i] << "_real," << name[i] << "_imag";
    }
    file << "\n";

    if (model == "SM") {
        MemoryManager::GetInstance("Test/InputFiles/testinput_thdm.lha", {0})->init();
    }
    else if (model == "THDM") {
        MemoryManager::GetInstance("Test/InputFiles/testinput_thdm.lha", {0,2})->init();
    }
    else if (model == "SUSY") {
        MemoryManager::GetInstance("Test/InputFiles/testInput.slha", {0,1})->init();
    }
    else {
        LOG_ERROR("ModelError", "MODEL not known");
    }
    Parameters* sm = Parameters::GetInstance();
    WilsonManager* wm = WilsonManager::GetInstance(strat_name, 81.0, strategy);

    double answer = 42.;
    wm->setScale(answer);

    double alpha_s = (*sm).QCDRunner.runningAlphasCalculation(answer);

    file << answer << "," << alpha_s;

            for (int i = 10; i <= 23; ++i) {
        complex_t C = {0,0};

        if (strat_name == "LO"){
            LOG_INFO(i, wm->get(static_cast<WilsonCoefficient>(i), 0));
            C = wm->get(static_cast<WilsonCoefficient>(i), 0);}
        else if (strat_name == "NLO")
            C = wm->get(static_cast<WilsonCoefficient>(i), 1);
        else if (strat_name == "NNLO")
            C = wm->get(static_cast<WilsonCoefficient>(i), 2);
        file << "," << C.real() << "," << C.imag();
    }

    file << "\n"; 

    file.close();

    wm->Cleanup();
}

void runTest(const std::string& strategyName, const std::shared_ptr<InitializationStrategy>& strategy, const std::string& testFile, const std::string& referenceFile,const std::string& model, double tolerance, bool primeCQ) {
    if (primeCQ) {
        writeCoefficientsPrimeCQToFile(strategyName, testFile, strategy, 81, model);
    }
    else {
        writeCoefficientsToFile(strategyName, testFile, strategy, 81, model);
    }
    
    if (!compareCSV(testFile, referenceFile, tolerance)) {
        std::cerr << "Test failed for " << strategyName << std::endl;
        exit(EXIT_FAILURE);
    } else {
        std::cout << "Test passed for " << strategyName << std::endl;
    }
}