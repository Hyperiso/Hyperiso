#include <iostream>
#include "WilsonUtils.h"
#include "config.hpp"



int main() {
    std::string root_file = project_root.data();

    auto loStrategy = std::make_shared<SUSY_LO_Strategy>();
    auto nloStrategy = std::make_shared<SUSY_NLO_Strategy>();
    auto nnloStrategy = std::make_shared<SUSY_NNLO_Strategy>();

    double tolerance = 0.01;

    runTest("LO", loStrategy, root_file + "/Test/csv/susy/WilsonCoefficients_LO.csv", root_file + "/Test/csv/superiso/susy/WilsonCoefficients_LO.csv", "SUSY", tolerance);
    runTest("NLO", nloStrategy, root_file + "/Test/csv/susy/WilsonCoefficients_NLO.csv", root_file + "/Test/csv/superiso/susy/WilsonCoefficients_NLO.csv", "SUSY",tolerance);
    runTest("NNLO", nnloStrategy, root_file + "/Test/csv/susy/WilsonCoefficients_NNLO.csv", root_file + "/Test/csv/superiso/susy/WilsonCoefficients_NNLO.csv", "SUSY",tolerance);

    // Ajoutez d'autres tests pour SM, THDM, etc.

    return 0;
}