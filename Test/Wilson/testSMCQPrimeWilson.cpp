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
    
    runTest("LO", root_file + "/Test/csv/sm/WilsonCoefficients_PRIMECQ_LO.csv", root_file + "/Test/csv/superiso/sm/WilsonCoefficients_PRIMECQ_LO.csv", "SM", tolerance, true);
    // runTest("NLO", nloStrategy, root_file + "/Test/csv/sm/WilsonCoefficients_PRIMECQ_NLO.csv", root_file + "/Test/csv/superiso/sm/WilsonCoefficients_PRIMECQ_NLO.csv", "SM",tolerance, false);
    // runTest("NNLO", nnloStrategy, root_file + "/Test/csv/sm/WilsonCoefficients_PRIMECQ_NNLO.csv", root_file + "/Test/csv/superiso/sm/WilsonCoefficients_PRIMECQ_NNLO.csv", "SM",tolerance, false);

    return 0;
}
