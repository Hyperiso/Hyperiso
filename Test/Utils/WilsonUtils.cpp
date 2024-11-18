#include "WilsonUtils.h"
#include <fstream>
#include "MemoryManager.h"
#include "Logger.h"
#include "CompareCsv.h"

void writeCoefficientsToFile(const std::string& strat_name, const std::string& fileName, double Q_match, const std::string& model) {
    std::ofstream file(fileName);

    file << "Q,alpha_s";
    for (int i = 1; i <= 10; ++i) {
        file << ",C" << i << "_real,C" << i << "_imag";
    }
    file << "\n";

    CoefficientManager* wm = CoefficientManager::GetInstance(model);

    if (model == "SM") {
        MemoryManager::GetInstance("Test/InputFiles/testinput_thdm.lha", {0})->init();
        wm->registerCoefficientGroup("BCoefficient", std::make_shared<BCoefficientGroup>());
    }
    else if (model == "THDM") {
        MemoryManager::GetInstance("Test/InputFiles/testinput_thdm.lha", {0,2})->init();
        wm->registerCoefficientGroup("BCoefficient", std::make_shared<BCoefficientGroup_THDM>());
    }
    else if (model == "SUSY") {
        MemoryManager::GetInstance("Test/InputFiles/testInput.slha", {0,1})->init();
        wm->registerCoefficientGroup("BCoefficient", std::make_shared<BCoefficientGroup_susy>());
    }
    else {
        LOG_ERROR("ModelError", "MODEL not known");
    }
    Parameters* sm = Parameters::GetInstance();
    // WilsonManager* wm = WilsonManager::GetInstance(strat_name, 81.0, strategy);
    

    
    wm->setQMatch("BCoefficient", Q_match);
    wm->setMatchingCoefficient("BCoefficient", strat_name);
    double alpha_s = (*sm).alpha_s(Q_match);

    file << Q_match << "," << alpha_s;
    std::vector<std::string> name {"C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"};
    for (auto& coeff : name) {
        complex_t C = {0,0};
        C = wm->getMatchingCoefficient("BCoefficient", coeff, strat_name);
        std::cout << coeff << " " << C << std::endl;
        file << "," << C.real() << "," << C.imag();
    }

    file << "\n"; 

    file.close();

    wm->Cleanup();
}

void writeCoefficientsPrimeCQToFile(const std::string& strat_name, const std::string& fileName, double Q_match, const std::string& model) {
    std::ofstream file(fileName);

    std::vector<std::string> name {"C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "CQ1", "CQ2", "CP1", "CP2", "CP3", "CP4", "CP5", "CP6", "CP7", "CP8", "CP9", "CP10", "CPQ1", "CPQ2"};

    file << "Q,alpha_s";
    for (int i = 10; i <= 23; ++i) {
        file << "," << name[i] << "_real," << name[i] << "_imag";
    }
    file << "\n";

    CoefficientManager* wm = CoefficientManager::GetInstance(model);

    if (model == "SM") {
        MemoryManager::GetInstance("Test/InputFiles/testinput_thdm.lha", {0})->init();
        wm->registerCoefficientGroup("BPrimeCoefficient", std::make_shared<BPrimeCoefficientGroup>());
        wm->registerCoefficientGroup("BScalarCoefficient", std::make_shared<BScalarCoefficientGroup>());    
    }
    else if (model == "THDM") {
        MemoryManager::GetInstance("Test/InputFiles/testinput_thdm.lha", {0,2})->init();
        wm->registerCoefficientGroup("BPrimeCoefficient", std::make_shared<BPrimeCoefficientGroup_THDM>());
        wm->registerCoefficientGroup("BScalarCoefficient", std::make_shared<BScalarCoefficientGroup_THDM>());
    }
    else if (model == "SUSY") {
        MemoryManager::GetInstance("Test/InputFiles/testInput.slha", {0,1})->init();
        wm->registerCoefficientGroup("BPrimeCoefficient", std::make_shared<BPrimeCoefficientGroup_susy>());
        wm->registerCoefficientGroup("BScalarCoefficient", std::make_shared<BScalarCoefficientGroup_susy>());
    }
    else {
        LOG_ERROR("ModelError", "MODEL not known");
    }
    Parameters* sm = Parameters::GetInstance();
    // WilsonManager* wm = WilsonManager::GetInstance(strat_name, 81.0, strategy);
    
    wm->setQMatch("BPrimeCoefficient", Q_match);
    wm->setMatchingCoefficient("BPrimeCoefficient", strat_name);
    wm->setQMatch("BScalarCoefficient", Q_match);
    wm->setMatchingCoefficient("BScalarCoefficient", strat_name);
    double answer = 42.;
    wm->setGroupScale("BPrimeCoefficient", answer);
    wm->setGroupScale("BScalarCoefficient", answer);
    wm->setRunCoefficient("BPrimeCoefficient", strat_name);
    wm->setRunCoefficient("BScalarCoefficient", strat_name);
    // wm->setScale(answer);

    double alpha_s = (*sm).alpha_s(answer);

    file << answer << "," << alpha_s;

    std::vector<std::string> name_scalar {"CQ1", "CQ2"};
    std::vector<std::string> name_prime {"CP1", "CP2", "CP3", "CP4", "CP5", "CP6", "CP7", "CP8", "CP9", "CP10", "CPQ1", "CPQ2"};

    for (auto& coeff : name_scalar) {
        complex_t C = {0.,0.};
        C = wm->getRunCoefficient("BScalarCoefficient", coeff, strat_name);
        file << "," << C.real() << "," << C.imag();
    }
    for (auto& coeff : name_prime) {
        complex_t C = {0.,0.};
        C = wm->getRunCoefficient("BPrimeCoefficient", coeff, strat_name);
        file << "," << C.real() << "," << C.imag();
    }

    file << "\n"; 

    file.close();

    wm->Cleanup();
}

void writeRunCoefficientsToFile(const std::string& strat_name, const std::string& fileName, double Q_match, double Q, const std::string& model, int base) {
    std::ofstream file(fileName);

    file << "Q,alpha_s";
    for (int i = 1; i <= 10; ++i) {
        file << ",C" << i << "_real,C" << i << "_imag";
    }
    file << "\n";

    CoefficientManager* wm;

    if (model == "SM") {
        MemoryManager::GetInstance("Test/InputFiles/testinput_thdm.lha", {0})->init();
        std::map<std::string, std::shared_ptr<CoefficientGroup>> temp_map;
        temp_map["BCoefficient"] = std::make_shared<BCoefficientGroup>();
        wm = CoefficientManager::Builder(model, temp_map, Q_match, Q, strat_name);
    }
    else if (model == "THDM") {
        MemoryManager::GetInstance("Test/InputFiles/testinput_thdm.lha", {0,2})->init();
        std::map<std::string, std::shared_ptr<CoefficientGroup>> temp_map;
        temp_map["BCoefficient"] = std::make_shared<BCoefficientGroup_THDM>();
        wm = CoefficientManager::Builder(model, temp_map, Q_match, Q, strat_name);
    }
    else if (model == "SUSY") {
        MemoryManager::GetInstance("Test/InputFiles/testInput.slha", {0,1})->init();
        std::map<std::string, std::shared_ptr<CoefficientGroup>> temp_map;
        temp_map["BCoefficient"] = std::make_shared<BCoefficientGroup_susy>();
        wm = CoefficientManager::Builder(model, temp_map, Q_match, Q, strat_name);
    }
    else {
        LOG_ERROR("ModelError", "MODEL not known");
    }
    Parameters* sm = Parameters::GetInstance();
    // WilsonManager* wm = WilsonManager::GetInstance(strat_name, 81.0, strategy);

    if (base==2){
        
        wm->switchbasis("BCoefficient");
    }
    
    double alpha_s = (*sm).alpha_s(Q_match);

    file << Q << "," << alpha_s;
    std::vector<std::string> name {"C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"};
    for (auto& coeff : name) {
        complex_t C = {0,0};
        C = wm->getRunCoefficient("BCoefficient", coeff, strat_name);
        std::cout << coeff << " " << C << std::endl;
        file << "," << C.real() << "," << C.imag();
    }

    file << "\n"; 

    file.close();

    wm->Cleanup();
}

void runTest(const std::string& strategyName, const std::string& testFile, const std::string& referenceFile,const std::string& model, double tolerance, bool primeCQ, int run) {
    
    
    if (run==1 || run==2) {
        writeRunCoefficientsToFile(strategyName, testFile, 81, 42, model, run);
    }
    else if (run == 0) {
        if (primeCQ) {
            writeCoefficientsPrimeCQToFile(strategyName, testFile, 81, model);
        }
        else {
            writeCoefficientsToFile(strategyName, testFile, 81, model);
        }
    }
    else {
        LOG_ERROR("ValueError", "Value error for base choice");
    }

    if (!compareCSV(testFile, referenceFile, tolerance)) {
        std::cerr << "Test failed for " << strategyName << std::endl;
        exit(EXIT_FAILURE);
    } else {
        std::cout << "Test passed for " << strategyName << std::endl;
    }
}