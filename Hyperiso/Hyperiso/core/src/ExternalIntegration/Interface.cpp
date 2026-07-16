#include "Interface.h"

#include <stdexcept>

#include "2HDMC.h"
#include "SoftSusy.h"

std::shared_ptr<ICalculator> GeneralCalculatorFactory::createCalculator(CalculatorType type) {
    switch(type) {
        case CalculatorType::Softsusy:
            return std::make_shared<SoftsusyCalculator>();
        case CalculatorType::TwoHDM:
            return std::make_shared<TwoHDMCalculator>();
        default:
            throw std::invalid_argument("Invalid Calculator Type");
    }
}

void GeneralCalculatorFactory::executeCommand(CalculatorType type, const std::string& commandName, const std::string& inputFilePath, const std::string& outputFilePath) {
    using CommandFunction = std::function<void(const std::string&, const std::string&)>;
    static std::unordered_map<std::string, CommandFunction> commonCommands = {
        {"calculateSpectrum", [type](const std::string& input, const std::string& output) {
            auto calculator = createCalculator(type);
            calculator->calculateSpectrum(input, output);
        }}
    };

    auto commandIterator = commonCommands.find(commandName);
    if (commandIterator != commonCommands.end()) {
        commandIterator->second(inputFilePath, outputFilePath);
    } else {
        throw std::invalid_argument("Calculator command '" + commandName + "' is not recognized.");
    }
}
