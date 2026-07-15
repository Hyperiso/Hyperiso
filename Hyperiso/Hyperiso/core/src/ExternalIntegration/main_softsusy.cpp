#include <iostream>
#include "Logger.h"
#include "SoftSusy.h"


int main() {
    Logger* logger = Logger::getInstance();
    logger->setLevel(Logger::LogLevel::INFO);

    std::string inputPath = "External_Integration/softsusy_ewinos_example.in";
    std::string outputPath = "output.slha";

    LOG_INFO("Starting spectrum calculation...");
    SoftsusyCalculator().calculateSpectrum( inputPath, outputPath);
    LOG_INFO("Spectrum calculation completed.");

    return 0;
}