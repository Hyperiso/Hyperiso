#include <iostream>
#include "../Core/Logger.h"
#include "SoftSusy.h"


int main() {
    // Initialize the logger (example setup)
    Logger* logger = Logger::getInstance();
    logger->setLevel(Logger::LogLevel::INFO);
    // Optionally, set a log file
    // logger->setLogFile("calculation.log");

    // Factory creates a calculator instance
    ICalculator* calculator = SoftsusyCalculatorFactory::createCalculator();
    
    // Assuming you have prepared "input.slha" and expect output in "output.slha"
    std::string inputPath = "../External_Integration/softsusy_ewinos_example.in";
    std::string outputPath = "output.slha";

    // Create the command to calculate the spectrum
    CalculateSpectrumCommand command(*static_cast<SoftsusyCalculator*>(calculator), inputPath, outputPath);

    // Execute the command
    logger->info("Starting spectrum calculation...");
    command.execute();
    logger->info("Spectrum calculation completed.");

    // Cleanup
    delete calculator;

    return 0;
}