#include <iostream>
#include <fstream>
#include <cstdlib> // For system()
#include "../Core/Logger.h"
#include "SoftSusy.h"

/**
 * Creates a new calculator instance.
 * This method implements the Factory Method pattern, providing a way to instantiate
 * ICalculator objects. It currently defaults to creating instances of SoftsusyCalculator.
 * 
 * @return ICalculator* Pointer to the newly created ICalculator instance.
 */
ICalculator* SoftsusyCalculatorFactory::createCalculator() {
    return new SoftsusyCalculator();
}


/**
 * @brief Executes the spectrum calculation using SOFTSUSY.
 * 
 * @param inputFilePath Path to the input file containing model parameters.
 * @param outputFilePath Path where the calculated spectrum will be written.
 */
void SoftsusyCalculator::calculateSpectrum(const std::string& inputFilePath, const std::string& outputFilePath) {
    Logger* logger = Logger::getInstance();
    // Example system call to SOFTSUSY - replace with actual implementation
    std::string command = "./External_Integration/softsusy-4.1.13/softpoint.x leshouches < " + inputFilePath + " > " + outputFilePath;

    logger->debug("SOFTSUSY COMMAND : " +command);

    int result = system(command.c_str());
    if (result != 0) {
        logger->error("SOFTSUSY execution failed with code " + std::to_string(result));
    } else {
        logger->info("SOFTSUSY execution successful.");
    }
}

// /**
//  * Execute the command to calculate the spectrum.
//  * This method calls upon the SoftsusyCalculator's calculateSpectrum method,
//  * passing in the previously specified input and output file paths, to perform
//  * the actual calculation. This encapsulation allows for the command to be executed
//  * without knowledge of the calculation specifics.
//  */
// void CalculateSpectrumCommand::execute() {
//     calculator.calculateSpectrum(inputFilePath, outputFilePath);
// }