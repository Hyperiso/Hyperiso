#include <iostream>
#include <cstdlib>
#include "WilsonUtils.h"
#include "config.hpp"



int main() {
    Logger* logger = Logger::getInstance();

    if (std::getenv("DISABLE_LOGGER")) {
        logger->setEnabled(false);
    }
    
    std::string root_file = project_root.data();
    
    double tolerance = 0.01;
    
    runTest("LO", root_file + "Test/csv/sm/WilsonCoefficients_PRIMECQ_LO.csv", root_file + "Test/csv/superiso/sm/WilsonCoefficients_PRIMECQ_LO.csv", "SM", tolerance, true);

    return 0;
}
