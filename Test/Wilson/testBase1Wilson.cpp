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

    // auto loStrategy = std::make_shared<SM_LO_Strategy>();
    // auto nloStrategy = std::make_shared<SM_NLO_Strategy>();
    // auto nnloStrategy = std::make_shared<SM_NNLO_Strategy>();

    double tolerance = 0.01;
    
    runTest("LO", root_file + "/Test/csv/sm/Run1WilsonCoefficients_LO.csv", root_file + "/Test/csv/superiso/sm/Run1WilsonCoefficients_LO.csv", "SM", tolerance, false, 1);
    runTest("NLO", root_file + "/Test/csv/sm/Run1WilsonCoefficients_NLO.csv", root_file + "/Test/csv/superiso/sm/Run1WilsonCoefficients_NLO.csv", "SM",tolerance, false, 1);
    runTest("NNLO", root_file + "/Test/csv/sm/Run1WilsonCoefficients_NNLO.csv", root_file + "/Test/csv/superiso/sm/Run1WilsonCoefficients_NNLO.csv", "SM",tolerance, false, 1);

    return 0;
}