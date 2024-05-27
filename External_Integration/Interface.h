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