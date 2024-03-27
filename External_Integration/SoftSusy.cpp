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
// ICalculator* SoftsusyCalculatorFactory::createCalculator() {
//     return new SoftsusyCalculator();
// }


/**
 * @brief Executes the spectrum calculation using SOFTSUSY.
 * 
 * @param inputFilePath Path to the input file containing model parameters.
 * @param outputFilePath Path where the calculated spectrum will be written.
 */
void SoftsusyCalculator::calculateSpectrum(const std::string& inputFilePath, const std::string& outputFilePath) {
    Logger* logger = Logger::getInstance();
    // Example system call to SOFTSUSY - replace with actual implementation
    std::string command = "./External_Integration/SOFTSUSY/src/SOFTSUSY/softpoint.x leshouches < " + inputFilePath + " > " + outputFilePath;

    logger->debug("SOFTSUSY COMMAND : " +command);

    int result = system(command.c_str());
    if (result != 0) {
        logger->error("SOFTSUSY execution failed with code " + std::to_string(result));
    } else {
        logger->info("SOFTSUSY execution successful.");
    }
}
