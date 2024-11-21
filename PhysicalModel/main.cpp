#include <iostream>
#include <complex>
#include "MartyWilson.h"
#include "CSVReader.h"
#include "DataFrame.h"
#include "WilsonManager.h"

int main() {
    MemoryManager::GetInstance()->init("Test/InputFiles/testinput_thdm.lha", {0});
    std::string csv_path = "../../DataBase/MartyWilson/SM_wilson.csv";

    std::string coefficient_name = "C7";
    double Q_match_value = 201.0;

    try {
        MartyWilson marty_coefficient(Q_match_value, coefficient_name, csv_path);
        std::complex<double> result = marty_coefficient.LO_calculation();
        std::cout << "Coefficient " << coefficient_name << " at Q_match = " << Q_match_value << ": "
                  << result.real() << " + " << result.imag() << "i" << std::endl;
    } catch (const std::exception& ex) {
        std::cerr << "An error occurred: " << ex.what() << std::endl;
    }

    CoefficientManager* wm = CoefficientManager::GetInstance("SM");

    wm->registerCoefficientGroup("BCoefficient", std::make_shared<BCoefficientGroupMarty>());

    wm->setQMatch("BCoefficient", 81);
    std::cout << "mass_top before : " << wm->get_params("MASS", 6) << std::endl;
    wm->setParams("BCoefficient", "MASS", 6, 150);
    std::cout << "mass_top after : " << wm->get_params("MASS", 6) << std::endl;

    wm->setMatchingCoefficient("BCoefficient", "LO");
    auto group = wm->getCoefficientGroup("BCoefficient");
    for (auto& elem : *group) {
        std::cout << elem.first << std::endl;
        std::cout << elem.second->get_CoefficientMatchingValue("LO")<< std::endl;
    }
    return 0;
}
