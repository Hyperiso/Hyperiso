#include <iostream>
#include <cstdlib>
#include "WilsonUtils.h"
#include "config.hpp"
#include "Logger.h"


int main() {
    Logger* logger = Logger::getInstance();

    if (std::getenv("DISABLE_LOGGER")) {
        logger->setEnabled(false);
    }
    

    std::string root_file = project_root.data();
    LOG_INFO("BONCAMARCHEOUPAS");
    // auto loStrategy = std::make_shared<SUSY_LO_Strategy>();
    // auto nloStrategy = std::make_shared<SUSY_NLO_Strategy>();
    // auto nnloStrategy = std::make_shared<SUSY_NNLO_Strategy>();

    double tolerance = 0.01;

    runTest("LO", root_file + "/Test/csv/susy/WilsonCoefficients_LO.csv", root_file + "/Test/csv/superiso/susy/WilsonCoefficients_LO.csv", "SUSY", tolerance);
    runTest("NLO", root_file + "/Test/csv/susy/WilsonCoefficients_NLO.csv", root_file + "/Test/csv/superiso/susy/WilsonCoefficients_NLO.csv", "SUSY",tolerance);
    runTest("NNLO", root_file + "/Test/csv/susy/WilsonCoefficients_NNLO.csv", root_file + "/Test/csv/superiso/susy/WilsonCoefficients_NNLO.csv", "SUSY",tolerance);


    return 0;
}