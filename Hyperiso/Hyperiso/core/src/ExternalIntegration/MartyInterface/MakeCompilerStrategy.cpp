#include "MakeCompilerStrategy.h"
#include <cstdlib>
#include <iostream>
#include <string>

void MakeCompilerStrategy::compile_run(const std::string& sourceFile, const std::string& outputBinary) {
    if (!this->check_if_compile(outputBinary)) {
        this->compile(sourceFile, outputBinary);
    }
    std::string command_run =  outputBinary + " -Q " + std::to_string(this->Q_match);
    executeCommand(command_run);
}

void MakeCompilerStrategy::compile(const std::string& sourceFile, const std::string& outputBinary) {
    std::string command = "cd " + sourceFile + " && make";
    executeCommand(command);
}