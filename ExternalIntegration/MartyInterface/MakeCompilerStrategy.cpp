#include "MakeCompilerStrategy.h"
#include <cstdlib>

void MakeCompilerStrategy::compile_run(const std::string& sourceFile, const std::string& outputBinary) {
    std::string command = "make";
    system(command.c_str());
}