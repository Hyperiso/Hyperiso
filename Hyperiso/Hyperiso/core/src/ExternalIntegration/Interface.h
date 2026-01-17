#ifndef INTERFACE_H
#define INTERFACE_H

#include <string>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <functional>

/**
 * @file Interface.h
 * @brief Declares generic calculator interfaces and the general calculator factory.
 *
 * This header provides a small command/factory framework used to wrap
 * external spectrum calculators (e.g. SOFTSUSY, 2HDMC) behind a common
 * interface:
 *
 * - ::ICalculator: abstract interface for any spectrum calculator.
 * - ::ICommand: minimal command interface.
 * - ::CalculateSpectrumCommand: concrete command that runs a calculator.
 * - ::GeneralCalculatorFactory: factory & dispatcher for calculators/commands.
 *
 * @ingroup SpectrumCalculationModule
 */

 /**
 * @enum CalculatorType
 * @ingroup SpectrumCalculationModule
 * @brief Enumerates the available concrete spectrum calculators.
 *
 * The mapping from ::Model to ::CalculatorType is done in ::SpectrumCalculator.
 */
enum class CalculatorType {
    Softsusy,   ///< Use SOFTSUSY–based calculator.
    TwoHDM      ///< Use 2HDMC–based calculator.
};


/**
 * @class ICalculator
 * @ingroup SpectrumCalculationModule
 * @brief Interface for calculator classes.
 *
 * Defines the minimal interface that all spectrum calculators have to
 * implement. A calculator is responsible for:
 *  - reading model parameters from an input file,
 *  - performing the spectrum calculation,
 *  - writing the resulting spectrum to an output file.
 *
 * Concrete implementations include:
 *  - ::SoftsusyCalculator
 *  - ::TwoHDMCalculator
 */
class ICalculator {
public:
    /**
     * @brief Performs the spectrum calculation.
     *
     * Implementations should:
     *  - read the input file at @p inputFilePath,
     *  - run the corresponding external tool / library,
     *  - write the resulting spectrum to @p outputFilePath.
     *
     * @param inputFilePath   Path to the input LHA file (or equivalent format).
     * @param outputFilePath  Path where the calculated spectrum must be stored.
     */
    virtual void calculateSpectrum(const std::string& inputFilePath, const std::string& outputFilePath) = 0;

    /// Virtual destructor.
    virtual ~ICalculator() = default;
    
};

/**
 * @class ICommand
 * @ingroup SpectrumCalculationModule
 * @brief Minimal command interface.
 *
 * Represents an executable operation that can be scheduled, queued, or
 * dispatched in a generic way. In this context, commands typically wrap
 * calls to a specific ::ICalculator.
 */
class ICommand {
public:
    /// Executes the command.
    virtual void execute() = 0;

    /// Virtual destructor.
    virtual ~ICommand() = default;
};


/**
 * @class CalculateSpectrumCommand
 * @ingroup SpectrumCalculationModule
 * @brief Command to calculate the spectrum using a given calculator.
 *
 * Encapsulates a request to calculate a particle spectrum, parameterized by:
 *  - a reference to an ::ICalculator implementation,
 *  - input / output file paths.
 *
 * This decouples the invocation site from the actual calculator type and
 * makes it easier to extend with other commands later.
 */
class CalculateSpectrumCommand : public ICommand {
private:
    ICalculator& calculator;    ///< Calculator instance to be used.
    std::string inputFilePath;  ///< Path to the input file.
    std::string outputFilePath; ///< Path to the output file.

public:
    /**
     * @brief Constructs a new CalculateSpectrumCommand.
     *
     * @param calculator       Reference to the calculator used for the spectrum.
     * @param inputFilePath    Path to the input model file.
     * @param outputFilePath   Path to the file where the spectrum will be written.
     */
    CalculateSpectrumCommand(ICalculator& calculator_, std::string inputFilePath_, std::string outputFilePath_)
        : calculator(calculator_), inputFilePath(std::move(inputFilePath_)), outputFilePath(std::move(outputFilePath_)) {}

    /**
     * @brief Executes the command to calculate the spectrum.
     *
     * Simply forwards the stored input/output paths to the underlying
     * ::ICalculator instance.
     */
    void execute() override {
        calculator.calculateSpectrum(inputFilePath, outputFilePath);
    }
};



/**
 * @class GeneralCalculatorFactory
 * @ingroup SpectrumCalculationModule
 * @brief Factory and dispatcher for spectrum calculators and commands.
 *
 * GeneralCalculatorFactory is responsible for:
 *  - instantiating concrete ::ICalculator implementations from a
 *    ::CalculatorType (see ::createCalculator()),
 *  - exposing a generic command dispatcher (::executeCommand()) that
 *    supports named commands (currently "calculateSpectrum").
 *
 * It is the main integration point between higher-level code
 * (e.g. ::SpectrumCalculator, ::MemoryManager) and the low-level
 * calculator backends.
 */
class GeneralCalculatorFactory {
public:
    /**
     * @brief Creates a concrete calculator for a given type.
     *
     * This method hides the compile-time conditional availability of
     * calculators (e.g. SOFTSUSY, 2HDMC) behind a single entry point.
     *
     * @param type  The requested calculator type.
     * @return Shared pointer to a concrete ::ICalculator instance.
     *
     * @throws std::invalid_argument if @p type is invalid or unsupported.
     */
    static std::shared_ptr<ICalculator> createCalculator(CalculatorType type);

    /**
     * @brief Executes a named command on a chosen calculator type.
     *
     * Currently supported @p commandName values:
     *  - `"calculateSpectrum"`: instantiate the requested calculator and
     *    call ::ICalculator::calculateSpectrum().
     *
     * @param type            Calculator type to use.
     * @param commandName     Name of the command (e.g. "calculateSpectrum").
     * @param inputFilePath   Path to the input file.
     * @param outputFilePath  Path to the output file.
     *
     * If the command name is unknown, an error message is printed to stderr.
     */
    static void executeCommand(CalculatorType type, const std::string& commandName, const std::string& inputFilePath, const std::string& outputFilePath);
};


#endif