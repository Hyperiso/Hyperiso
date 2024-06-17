#include <iostream>
#include "../Core/Logger.h"
#include "SoftSusy.h"


int main() {
    Logger* logger = Logger::getInstance();
    logger->setLevel(Logger::LogLevel::INFO);

    // Les chemins vers les fichiers d'entrée et de sortie
    std::string inputPath = "External_Integration/softsusy_ewinos_example.in";
    std::string outputPath = "output.slha";

    LOG_INFO("Starting spectrum calculation...");
    SoftsusyCalculatorFactory::executeCommand("calculateSpectrum", inputPath, outputPath);
    LOG_INFO("Spectrum calculation completed.");

    return 0;
}