#include <iostream>
#include <complex>
#include "MartyWilson.h"
#include "CSVReader.h"
#include "DataFrame.h"
#include "WilsonManager.h"

int main() {
    MemoryManager::GetInstance()->init("Test/InputFiles/testinput_thdm.lha", Model::SM,  true);
    std::string csv_path = "../../DataBase/MartyWilson/SM_wilson.csv";

    std::shared_ptr<CoefficientGroup> b_coeff = std::make_shared<BCoefficientGroup>(BCoefficientGroup(160.));

    b_coeff->init_LO();
    std::cout << b_coeff << std::endl;
}
