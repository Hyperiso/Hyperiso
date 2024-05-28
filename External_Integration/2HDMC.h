#pragma once
#include "Interface.h"
#include <unordered_map>
#include <functional>
#include <iostream>


class TwoHDMCalculator : public ICalculator {
public:
    void calculateSpectrum(const std::string& inputFilePath, const std::string& outputFilePath) override;
};



// class CalculateSpectrumCommand : public ICommand {
// private:
//     TwoHDMCalculator& calculator;
//     std::string inputFilePath;
//     std::string outputFilePath;

// public:
//     CalculateSpectrumCommand(TwoHDMCalculator& calculator, std::string inputFilePath, std::string outputFilePath)
//         : calculator(calculator), inputFilePath(std::move(inputFilePath)), outputFilePath(std::move(outputFilePath)) {}

//     void execute() override {
//         calculator.calculateSpectrum(inputFilePath, outputFilePath);
//     }
// };



class TwoHDMCalculatorFactory {
public:
    static ICalculator* createCalculator() {
        return new TwoHDMCalculator();
    }

    static void executeCommand(const std::string& commandName, const std::string& inputFilePath, const std::string& outputFilePath) {
        static std::unordered_map<std::string, std::function<void(const std::string&, const std::string&)>> commands = {
            {"calculateSpectrum", [](const std::string& input, const std::string& output) {
                TwoHDMCalculator* calculator = static_cast<TwoHDMCalculator*>(createCalculator());
                CalculateSpectrumCommand command(*calculator, input, output);
                command.execute();
                delete calculator;
            }}
        };

        auto commandIterator = commands.find(commandName);
        if (commandIterator != commands.end()) {
            commandIterator->second(inputFilePath, outputFilePath);
        } else {
            std::cerr << "Command '" << commandName << "' not recognized." << std::endl;
        }
    }
};