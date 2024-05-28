#include <string>
#include <functional> // Pour std::function
#include <unordered_map> 
#pragma once
#include "MemoryManager.h"
#include "Interface.h"


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


/**
 * @class SoftsusyCalculatorFactory
 * @brief Factory class to create instances of Softsusy calculators.
 *
 * Implements the Factory Method pattern to provide a way to instantiate
 * ICalculator objects. Currently defaults to creating SoftsusyCalculator instances.
 */
class SoftsusyCalculatorFactory {
public:
    /**
     * @brief Creates an ICalculator instance for SOFTSUSY calculations.
     * 
     * @return ICalculator* Pointer to the newly created ICalculator instance.
     */
    static ICalculator* createCalculator() {
        return new SoftsusyCalculator();
    }

    static void executeCommand(const std::string& commandName, const std::string& inputFilePath, const std::string& outputFilePath) {
        // Mapping of command names to their respective execution functions
        static std::unordered_map<std::string, std::function<void(const std::string&, const std::string&)>> commands = {
            {"calculateSpectrum", [](const std::string& input, const std::string& output) {
                SoftsusyCalculator* calculator = static_cast<SoftsusyCalculator*>(createCalculator());
                CalculateSpectrumCommand command(*calculator, input, output);
                command.execute();
                delete calculator; // Clean up the calculator instance after command execution
            }}
        };

        auto commandIterator = commands.find(commandName);
        if (commandIterator != commands.end()) {
            commandIterator->second(inputFilePath, outputFilePath); // Execute the found command
        } else {
            // Handle the case where the command is not found
            std::cerr << "Command '" << commandName << "' not recognized." << std::endl;
        }
    }
};






