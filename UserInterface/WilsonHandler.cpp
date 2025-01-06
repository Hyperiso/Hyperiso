#include "WilsonHandler.h"
#include <iostream>
#include <string>
#include <complex>
#include <map>
#include <memory>
#include <vector>
#include "WilsonManager.h"
#include "MartyWilson.h"

void print_wilson_usage() {
    std::cout << "Usage: ./main wilson [options]\n"
              << "\nOptions:\n"
              << "  --model/-m <model_name>               : Specify the model (SM, THDM, MSSM, ...)\n"
              << "  --wilson/-w <coefficient_name>        : Specify the Wilson coefficient (e.g., C1, C2, ..., CQ1)\n"
              << "  --group/-g <group_name>               : Specify the group if multiple coefficients are provided (e.g., BCoefficientGroup)\n"
              << "  --Q_match/-q <value>                  : Set Q_match scale\n"
              << "  --Q/-Q <value>                        : Set Q scale\n"
              << "  --n_flavor/-f <value> (optional)      : Set the number of flavors (optional)\n"
              << "  --Marty/-M <true|false>               : Use Marty groups (default: false)\n"
              << "  --input_file/-if <slha_name>               : input file for parameters spectrum\n"
              << "  --order/-o <LO|NLO|NNLO>              : Specify the calculation order (default: LO)\n"
              << "  --help/-h                             : Display this help message\n";
}

int handleWilsonOptions(int argc, char* argv[]) {
    std::string model_name;
    std::string coefficient_name;
    std::string group_name;
    double Q_match = 0.0;
    double Q = 0.0;
    int n_flavor = 5;
    bool use_marty = false;
    std::string order = "LO";
    std::string input_file = "Test/InputFiles/testinput_thdm.lha";
    std::vector<std::string> Bcoefficient = {"C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"};
    std::vector<std::string> BPrimecoefficient = {"CP1", "CP2", "CP3", "CP4", "CP5", "CP6", "CP7", "CP8", "CP9", "CP10", "CPQ1", "CPQ2"};
    std::vector<std::string> BScalarcoefficient = {"CQ1", "CQ2"};
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--model" || arg == "-m") && i + 1 < argc) {
            model_name = argv[++i];
        } else if ((arg == "--wilson" || arg == "-w") && i + 1 < argc) {
            coefficient_name = argv[++i];
        } else if ((arg == "--group" || arg == "-g") && i + 1 < argc) {
            group_name = argv[++i];
        } else if ((arg == "--Q_match" || arg == "-q") && i + 1 < argc) {
            Q_match = std::stod(argv[++i]);
        } else if ((arg == "--Q" || arg == "-Q") && i + 1 < argc) {
            Q = std::stod(argv[++i]);
        } else if ((arg == "--n_flavor" || arg == "-f") && i + 1 < argc) {
            n_flavor = std::stoi(argv[++i]);
        } else if ((arg == "--Marty" || arg == "-M") && i + 1 < argc) {
            use_marty = (std::string(argv[++i]) == "true");
        } else if ((arg == "--input_file" || arg == "-if") && i + 1 < argc) {
            input_file = argv[++i];
        } else if ((arg == "--order" || arg == "-o") && i + 1 < argc) {
            order = argv[++i];
            if (use_marty && order != "LO") {
                std::cerr << "Error: When --Marty is true, --order can only be LO." << std::endl;
                return 1;
            }
        } else if (arg == "--help" || arg == "-h") {
            print_wilson_usage();
            return 0;
        } else {
            std::cerr << "Unknown option: " << arg << std::endl;
            print_wilson_usage();
            return 1;
        }
    }

    try {
        MemoryManager::GetInstance()->init(input_file, Model::SM);
        auto manager = CoefficientManager::GetInstance(model_name);
        
        if (group_name.empty()) {
            if (std::find(Bcoefficient.begin(), Bcoefficient.end(), coefficient_name) != Bcoefficient.end()) {
                group_name = "BCoefficientGroup";
            } else if (std::find(BPrimecoefficient.begin(), BPrimecoefficient.end(), coefficient_name) != BPrimecoefficient.end()) {
                group_name = "BPrimeCoefficientGroup";
            } else if (std::find(BScalarcoefficient.begin(), BScalarcoefficient.end(), coefficient_name) != BScalarcoefficient.end()) {
                group_name = "BScalarCoefficientGroup";
            } else {
                throw std::invalid_argument("Invalid coefficient name: not found in any group.");
            }
        }
        if (use_marty) {
            group_name += "Marty";
        }
        std::cout << group_name << std::endl;
        if (group_name == "BCoefficientGroup" || group_name == "BCoefficientGroupMarty") {
            manager->registerCoefficientGroup(group_name, use_marty ? std::make_shared<BCoefficientGroupMarty>() : std::make_shared<BCoefficientGroup>());
        } else if (group_name == "BPrimeCoefficientGroup" || group_name == "BPrimeCoefficientGroupMarty") {
            manager->registerCoefficientGroup(group_name, use_marty ? std::make_shared<BPrimeCoefficientGroupMarty>() : std::make_shared<BPrimeCoefficientGroup>());
        } else if (group_name == "BScalarCoefficientGroup" || group_name == "BScalarCoefficientGroupMarty") {
            manager->registerCoefficientGroup(group_name, use_marty ? std::make_shared<BScalarCoefficientGroupMarty>() : std::make_shared<BScalarCoefficientGroup>());
        } else {
            throw std::invalid_argument("Invalid group name specified.");
        }

        manager->setQMatch(group_name, Q_match);
        manager->setMatchingCoefficient(group_name, order);
        manager->setGroupScale(group_name, Q);
        manager->setRunCoefficient(group_name, order);
        auto group = manager->getCoefficientGroup(group_name);
        if (coefficient_name.empty()) {
            for (const auto& elem : *group) {
                std::complex<double> coeff_M = elem.second->get_CoefficientMatchingValue(order);
                std::complex<double> coeff_Q = elem.second->get_CoefficientRunValue(order);
                std::cout << "Coefficient " << elem.first << " at Q_match = " << Q_match
                          << ": " << coeff_M.real() << " + " << coeff_M.imag() << "i\n";
                std::cout << "Coefficient " << elem.first << " at Q = " << Q
                          << ": " << coeff_Q.real() << " + " << coeff_Q.imag() << "i\n";
            }
        } else {
            auto it = group->find(coefficient_name);
            if (it != group->end()) {
                std::complex<double> coeff_M = it->second->get_CoefficientMatchingValue(order);
                std::complex<double> coeff_Q = it->second->get_CoefficientRunValue(order);
                std::cout << "Coefficient " << coefficient_name << " at Q_match = " << Q_match
                          << ": " << coeff_M.real() << " + " << coeff_M.imag() << "i\n";
                std::cout << "Coefficient " << coefficient_name << " at Q = " << Q
                          << ": " << coeff_Q.real() << " + " << coeff_Q.imag() << "i\n";
            } else {
                std::cerr << "Error: Specified coefficient not found in the group." << std::endl;
                return 1;
            }
        }
    } catch (const std::exception& ex) {
        std::cerr << "An error occurred: " << ex.what() << std::endl;
        return 1;
    }

    return 0;
}
