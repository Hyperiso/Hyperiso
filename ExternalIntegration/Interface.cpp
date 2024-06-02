#include "Interface.h"

#ifdef BUILD_WITH_SOFTSUSY
#include "SoftSusy.h"
#endif

#ifdef BUILD_WITH_2HDMC
#include "2HDMC.h"
#endif


std::unique_ptr<ICalculator> GeneralCalculatorFactory::createCalculator(CalculatorType type) {
    switch(type) {
        case CalculatorType::Softsusy:
            #ifdef BUILD_WITH_SOFTSUSY
            return std::make_unique<SoftsusyCalculator>();
            #else
            // throw std::runtime_error("SOFTSUSY support not enabled.");
            #endif
        case CalculatorType::TwoHDM:
            #ifdef BUILD_WITH_2HDMC
            return std::make_unique<TwoHDMCalculator>();
            #else
            // throw std::runtime_error("TwoHDM support not enabled.");
            #endif
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
        std::cerr << "Command '" << commandName << "' not recognized." << std::endl;
    }
}
