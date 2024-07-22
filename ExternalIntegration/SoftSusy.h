#ifdef BUILD_WITH_SOFTSUSY

#include <string>
#include <functional> // Pour std::function
#include <unordered_map> 
#pragma once
// #include "MemoryManager.h"
#include "Interface.h"
#include "config.hpp"

/**
 * @class SoftsusyCalculator
 * @brief Concrete implementation of ICalculator for SOFTSUSY.
 *
 * Implements the spectrum calculation using the SOFTSUSY library or executable.
 */
class SoftsusyCalculator : public ICalculator {
public:
    /**
     * Calculates the particle spectrum using SOFTSUSY with the given input and output file paths.
     * 
     * @param inputFilePath The path to the input file containing the model parameters.
     * @param outputFilePath The path where the calculated spectrum should be written.
     */
    void calculateSpectrum(const std::string& inputFilePath, const std::string& outputFilePath) override;
};


// /**
//  * @class CalculateSpectrumCommand
//  * @brief Command to calculate the spectrum using a given calculator.
//  *
//  * Encapsulates a request to calculate the particle spectrum, allowing the
//  * command to be parameterized with different calculators and input/output paths.
//  */
// class CalculateSpectrumCommand : public ICommand {
// private:
//     SoftsusyCalculator& calculator;
//     std::string inputFilePath;
//     std::string outputFilePath;

// public:
//     /**
//      * Construct a new Calculate Spectrum Command object.
//      * 
//      * @param calculator Reference to the SoftsusyCalculator.
//      * @param inputFilePath Path to the input file.
//      * @param outputFilePath Path to the output file.
//      */
//     CalculateSpectrumCommand(SoftsusyCalculator& calculator, std::string inputFilePath, std::string outputFilePath)
//         : calculator(calculator), inputFilePath(std::move(inputFilePath)), outputFilePath(std::move(outputFilePath)) {}

//     /**
//      * Execute the command to calculate the spectrum.
//      */
//     void execute() override {
//         calculator.calculateSpectrum(inputFilePath, outputFilePath);
//     }
// };








#endif
