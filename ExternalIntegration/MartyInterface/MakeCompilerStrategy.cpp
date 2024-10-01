#include "MakeCompilerStrategy.h"
#include <cstdlib>
#include <iostream>

void MakeCompilerStrategy::compile_run(const std::string& sourceFile, const std::string& outputBinary) {
    // std::cout << sourceFile << "  " << outputBinary << std::endl;
    if (!this->check_if_compile(outputBinary)) {
        this->compile(sourceFile, outputBinary);
    }
    
    std::string command_run = "./" + sourceFile + outputBinary;
    system(command_run.c_str());
}

void MakeCompilerStrategy::compile(const std::string& sourceFile, const std::string& outputBinary) {
    std::string command = "cd " + sourceFile + " && make";
    system(command.c_str());
}