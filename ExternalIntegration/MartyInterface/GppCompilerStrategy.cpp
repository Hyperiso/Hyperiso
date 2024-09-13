#include "GppCompilerStrategy.h"
#include <cstdlib>

void GppCompilerStrategy::compile_run(const std::string& sourceFile, const std::string& outputBinary) {
    std::string command_compile = "g++ -o " + outputBinary + " " + sourceFile + " -lmarty" + " -lgfortran";
    system(command_compile.c_str());
    std::string command_run = "./"+outputBinary;
    system(command_run.c_str());
}
