#include <iostream>
#include <complex>
#include "MartyWilson.h"
#include "CSVReader.h"
#include "DataFrame.h"

int main() {
    MemoryManager::GetInstance("Test/InputFiles/testinput_thdm.lha", {0})->init();
    std::string csv_path = "../../DataBase/MartyWilson/SM_wilson.csv";

    std::string coefficient_name = "C7";
    double Q_match_value = 200.0;

    try {
        MartyWilson marty_coefficient(Q_match_value, coefficient_name, csv_path);
        std::complex<double> result = marty_coefficient.LO_calculation();
        std::cout << "Coefficient " << coefficient_name << " at Q_match = " << Q_match_value << ": "
                  << result.real() << " + " << result.imag() << "i" << std::endl;
    } catch (const std::exception& ex) {
        std::cerr << "An error occurred: " << ex.what() << std::endl;
    }

    return 0;
}
