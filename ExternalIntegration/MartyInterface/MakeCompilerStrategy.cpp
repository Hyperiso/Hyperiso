#include "MakeCompilerStrategy.h"
#include <cstdlib>
#include <iostream>

void MakeCompilerStrategy::compile_run(const std::string& sourceFile, const std::string& outputBinary) {
    std::cout << sourceFile << "  " << outputBinary << std::endl;
    
    std::string command = "cd " + sourceFile + " && make";
    system(command.c_str());
    std::string command_run = "./" + sourceFile + outputBinary;
    system(command_run.c_str());
}