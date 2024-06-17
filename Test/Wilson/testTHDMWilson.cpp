#include <iostream>
#include "WilsonUtils.h"
#include "config.hpp"



int main() {
    std::string root_file = project_root.data();

    auto loStrategy = std::make_shared<THDM_LO_Strategy>();
    auto nloStrategy = std::make_shared<THDM_NLO_Strategy>();
    auto nnloStrategy = std::make_shared<THDM_NNLO_Strategy>();

    double tolerance = 0.01;

    runTest("LO", loStrategy, root_file + "/Test/csv/thdm/WilsonCoefficients_LO.csv", root_file + "/Test/csv/superiso/thdm/WilsonCoefficients_LO.csv", "THDM", tolerance);
    runTest("NLO", nloStrategy, root_file + "/Test/csv/thdm/WilsonCoefficients_NLO.csv", root_file + "/Test/csv/superiso/thdm/WilsonCoefficients_NLO.csv", "THDM",tolerance);
    runTest("NNLO", nnloStrategy, root_file + "/Test/csv/thdm/WilsonCoefficients_NNLO.csv", root_file + "/Test/csv/superiso/thdm/WilsonCoefficients_NNLO.csv", "THDM",tolerance);

    // Ajoutez d'autres tests pour SM, THDM, etc.

    return 0;
}