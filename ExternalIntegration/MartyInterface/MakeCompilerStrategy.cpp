#include "MakeCompilerStrategy.h"
#include <cstdlib>
#include <iostream>
#include <string>

void MakeCompilerStrategy::compile_run(const std::string& sourceFile, const std::string& outputBinary) {
    if (!this->check_if_compile(outputBinary)) {
        this->compile(sourceFile, outputBinary);
    }
    std::string command_run =  outputBinary + " -Q " + std::to_string(this->Q_match);
    std::cout << "------------------------------" << std::endl;
    std::cout << command_run << std::endl;
    system(command_run.c_str());
}

void MakeCompilerStrategy::compile(const std::string& sourceFile, const std::string& outputBinary) {
    std::string command = "cd " + sourceFile + " && make";
    system(command.c_str());
}