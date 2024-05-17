#include "2HDMC.h"

int main() {

    std::string commandName = "calculateSpectrum";
    std::string inputFilePath = "Test/LH_sample_input.lha";
    std::string outputFilePath = "Test/output.lha";

    TwoHDMCalculatorFactory::executeCommand(commandName, inputFilePath, outputFilePath);

    return 0;

    return 0;
}