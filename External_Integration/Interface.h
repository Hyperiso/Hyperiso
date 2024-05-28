#pragma once
#include <string>

/**
 * @class ICalculator
 * @brief Interface for calculator classes.
 *
 * Defines the interface for calculators, ensuring they can process input
 * and output paths for spectrum calculations.
 */
class ICalculator {
public:
    virtual void calculateSpectrum(const std::string& inputFilePath, const std::string& outputFilePath) = 0;
    virtual ~ICalculator() = default;
    
};


class ICommand {
public:
    virtual void execute() = 0;
    virtual ~ICommand() = default;
};


/**
 * @class CalculateSpectrumCommand
 * @brief Command to calculate the spectrum using a given calculator.
 *
 * Encapsulates a request to calculate the particle spectrum, allowing the
 * command to be parameterized with different calculators and input/output paths.
 */
class CalculateSpectrumCommand : public ICommand {
private:
    ICalculator& calculator;
    std::string inputFilePath;
    std::string outputFilePath;

public:
    /**
     * Construct a new Calculate Spectrum Command object.
     * 
     * @param calculator Reference to the SoftsusyCalculator.
     * @param inputFilePath Path to the input file.
     * @param outputFilePath Path to the output file.
     */
    CalculateSpectrumCommand(ICalculator& calculator, std::string inputFilePath, std::string outputFilePath)
        : calculator(calculator), inputFilePath(std::move(inputFilePath)), outputFilePath(std::move(outputFilePath)) {}

    /**
     * Execute the command to calculate the spectrum.
     */
    void execute() override {
        calculator.calculateSpectrum(inputFilePath, outputFilePath);
    }
};
