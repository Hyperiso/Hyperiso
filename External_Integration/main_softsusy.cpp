#include <iostream>
#include "../Core/Logger.h"
#include "SoftSusy.h"


int main() {
    Logger* logger = Logger::getInstance();
    logger->setLevel(Logger::LogLevel::INFO);

    // Les chemins vers les fichiers d'entrée et de sortie
    std::string inputPath = "External_Integration/softsusy_ewinos_example.in";
    std::string outputPath = "output.slha";

    logger->info("Starting spectrum calculation...");
    SoftsusyCalculatorFactory::executeCommand("calculateSpectrum", inputPath, outputPath);
    logger->info("Spectrum calculation completed.");

    return 0;
}