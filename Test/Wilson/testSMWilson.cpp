#include <iostream>
#include "CompareCsv.h"
#include "WilsonUtils.h"
#include "config.hpp"

void runTest(const std::string& strategyName, const std::shared_ptr<InitializationStrategy>& strategy, const std::string& testFile, const std::string& referenceFile,const std::string& model, double tolerance) {
    writeCoefficientsToFile(strategyName, testFile, strategy, 81, model);
    if (!compareCSV(testFile, referenceFile, tolerance)) {
        std::cerr << "Test failed for " << strategyName << std::endl;
        exit(EXIT_FAILURE);
    } else {
        std::cout << "Test passed for " << strategyName << std::endl;
    }
}

int main() {
    std::string root_file = project_root.data();

    auto loStrategy = std::make_shared<SM_LO_Strategy>();
    auto nloStrategy = std::make_shared<SM_NLO_Strategy>();
    auto nnloStrategy = std::make_shared<SM_NNLO_Strategy>();

    double tolerance = 0.01;

    runTest("LO", loStrategy, root_file + "/Test/csv/sm/WilsonCoefficients_LO.csv", root_file + "/Test/csv/superiso/sm/WilsonCoefficients_LO.csv", "SM", tolerance);
    runTest("NLO", nloStrategy, root_file + "/Test/csv/sm/WilsonCoefficients_NLO.csv", root_file + "/Test/csv/superiso/sm/WilsonCoefficients_NLO.csv", "SM",tolerance);
    runTest("NNLO", nnloStrategy, root_file + "/Test/csv/sm/WilsonCoefficients_NNLO.csv", root_file + "/Test/csv/superiso/sm/WilsonCoefficients_NNLO.csv", "SM",tolerance);

    // Ajoutez d'autres tests pour SM, THDM, etc.

    return 0;
}
