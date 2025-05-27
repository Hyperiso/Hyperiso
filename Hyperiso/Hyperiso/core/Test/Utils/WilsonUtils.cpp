#include "WilsonUtils.h"
#include <fstream>
#include "MemoryManager.h"
#include "WilsonInterface.h"
#include "HyperisoMaster.h"
#include "Logger.h"
#include "config.hpp"
#include "CompareCsv.h"
#include "Configs.h"
#include <unordered_set>


void writeCoefficientsToFile(const std::string& strat_name, const std::string& fileName, double Q_match, const std::string& model) {
    std::ofstream file(fileName);

    std::string root_data_file = project_assets_root.data();

    file << "Q,alpha_s";
    for (int i = 1; i <= 10; ++i) {
        file << ",C" << i << "_real,C" << i << "_imag";
    }
    file << "\n";
    std::shared_ptr<HyperisoMaster> hyperiso = std::make_shared<HyperisoMaster>();
    std::shared_ptr<WilsonInterface> wi;
    Config config;
    if (model == "SM") {
        config.model = Model::SM;
        hyperiso->init(root_data_file + "Test/InputFiles/testinput_thdm.lha", config);
    }
    else if (model == "THDM") {
        config.model = Model::THDM;
        hyperiso->init(root_data_file + "Test/InputFiles/testinput_thdm.lha", config);
    }
    else if (model == "SUSY") {
        config.model = Model::SUSY;
        hyperiso->init(root_data_file + "Test/InputFiles/testInput.slha", config);
    }
    else {
        LOG_ERROR("ModelError", "MODEL not known");
    }
    wi = std::make_shared<WilsonInterface>();
    std::unordered_set<WGroup> groups = {WGroup::B};
    WilsonBuildConfig wbc = {groups, Q_match, Q_match, OrderMapper::enum_elt(strat_name)};
    wi->build(wbc);
    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
    double alpha_s = QCDHelper::alpha_s(Q_match);

    file << Q_match << "," << alpha_s;
    std::vector<std::string> name {"C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"};
    for (auto& coeff : name) {
        complex_t C = {0,0};
        if (model == "SM") {
            C = wi->getMatchingCoefficient(WGroup::B,WCoefMapper::enum_elt(coeff), OrderMapper::enum_elt(strat_name), ContributionType::SM );
        } else {
            C = wi->getMatchingCoefficient(WGroup::B,WCoefMapper::enum_elt(coeff), OrderMapper::enum_elt(strat_name), ContributionType::BSM );
        }
        std::cout << coeff << " " << C << std::endl;
        file << "," << C.real() << "," << C.imag();
    }

    file << "\n"; 

    file.close();
}

void writeCoefficientsPrimeCQToFile(const std::string& strat_name, const std::string& fileName, double Q_match, const std::string& model) {
    std::ofstream file(fileName);

    std::string root_data_file = project_assets_root.data();

    std::vector<std::string> name {"C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "CQ1", "CQ2", "CP1", "CP2", "CP3", "CP4", "CP5", "CP6", "CP7", "CP8", "CP9", "CP10", "CPQ1", "CPQ2"};

    file << "Q,alpha_s";
    for (int i = 10; i <= 23; ++i) {
        file << "," << name[i] << "_real," << name[i] << "_imag";
    }
    file << "\n";

    std::shared_ptr<HyperisoMaster> hyperiso = std::make_shared<HyperisoMaster>();
    std::shared_ptr<WilsonInterface> wi;

    Config config;
    if (model == "SM") {
        config.model = Model::SM;
        hyperiso->init(root_data_file + "Test/InputFiles/testinput_thdm.lha", config);
    }
    else if (model == "THDM") {
        config.model = Model::THDM;
        hyperiso->init(root_data_file + "Test/InputFiles/testinput_thdm.lha", config);
    }
    else if (model == "SUSY") {
        config.model = Model::SUSY;
        hyperiso->init(root_data_file + "Test/InputFiles/testInput.slha", config);
        wi = std::make_shared<WilsonInterface>();
    }
    else {
        LOG_ERROR("ModelError", "MODEL not known");
    }
    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
    wi = std::make_shared<WilsonInterface>();
    std::unordered_set<WGroup> groups = {WGroup::BPrime, WGroup::BScalar};
    double answer = 42.;
    WilsonBuildConfig wbc = {groups, Q_match, answer, OrderMapper::enum_elt(strat_name)};
    wi->build(wbc);

    double alpha_s = QCDHelper::alpha_s(answer);

    file << answer << "," << alpha_s;

    std::vector<std::string> name_scalar {"CQ1", "CQ2"};
    std::vector<std::string> name_prime {"CP1", "CP2", "CP3", "CP4", "CP5", "CP6", "CP7", "CP8", "CP9", "CP10", "CPQ1", "CPQ2"};

    for (auto& coeff : name_scalar) {
        complex_t C = {0.,0.};
        if (model == "SM") {
            C = wi->getRunCoefficient(WGroup::BScalar, WCoefMapper::enum_elt(coeff), OrderMapper::enum_elt(strat_name), ContributionType::SM);
        } else {
            C = wi->getRunCoefficient(WGroup::BScalar, WCoefMapper::enum_elt(coeff), OrderMapper::enum_elt(strat_name), ContributionType::BSM);
        }
        file << "," << C.real() << "," << C.imag();
    }
    for (auto& coeff : name_prime) {
        complex_t C = {0.,0.};
        if (model == "SM") {
            C = wi->getRunCoefficient(WGroup::BPrime, WCoefMapper::enum_elt(coeff), OrderMapper::enum_elt(strat_name), ContributionType::SM);
        } else {
            C = wi->getRunCoefficient(WGroup::BPrime, WCoefMapper::enum_elt(coeff), OrderMapper::enum_elt(strat_name), ContributionType::BSM);
        }
        
        file << "," << C.real() << "," << C.imag();
    }

    file << "\n"; 

    file.close();
}

void writeRunCoefficientsToFile(const std::string& strat_name, const std::string& fileName, double Q_match, double Q, const std::string& model, int base) {
    std::ofstream file(fileName);
    std::string root_data_file = project_assets_root.data();
    file << "Q,alpha_s";
    for (int i = 1; i <= 10; ++i) {
        file << ",C" << i << "_real,C" << i << "_imag";
    }
    file << "\n";

    std::shared_ptr<HyperisoMaster> hyperiso = std::make_shared<HyperisoMaster>();
    std::shared_ptr<WilsonInterface> wi;

    Config config;

    std::cout << "wtf " << std::endl;
    if (model == "SM") {
        config.model = Model::SM;
        hyperiso->init(root_data_file + "Test/InputFiles/testinput_thdm.lha", config);
    }
    else if (model == "THDM") {
        config.model = Model::THDM;
        hyperiso->init(root_data_file + "Test/InputFiles/testinput_thdm.lha", config);
    }
    else if (model == "SUSY") {
        config.model = Model::SUSY;
        hyperiso->init(root_data_file + "Test/InputFiles/testInput.slha", config);
    }
    else {
        LOG_ERROR("ModelError", "MODEL not known");
    }
    wi = std::make_shared<WilsonInterface>();
    wi->build({{WGroup::B}, Q_match, Q, OrderMapper::enum_elt(strat_name)});
    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
    WilsonBasis basis = WilsonBasis::B_STANDARD;
    if (base==2){
        basis = WilsonBasis::B_TRADITIONAL;
    }
    
    double alpha_s = QCDHelper::alpha_s(Q_match);

    file << Q << "," << alpha_s;
    std::vector<std::string> name {"C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"};
    for (auto& coeff : name) {
        complex_t C = {0,0};
        if (model == "SM") {
            C = wi->getRunCoefficient(WGroup::B, WCoefMapper::enum_elt(coeff), OrderMapper::enum_elt(strat_name), ContributionType::SM, basis);
        } else {
            C = wi->getRunCoefficient(WGroup::B, WCoefMapper::enum_elt(coeff), OrderMapper::enum_elt(strat_name), ContributionType::BSM, basis);
        }
        std::cout << coeff << " " << C << std::endl;
        file << "," << C.real() << "," << C.imag();
    }

    file << "\n"; 

    file.close();
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