#include "GppCompilerStrategy.h"
#include <cstdlib>

void GppCompilerStrategy::compile(const std::string& sourceFile, const std::string& outputBinary) {
    std::string command = "g++ -o " + outputBinary + " " + sourceFile + " -lmarty" + " -lgfortran";
    system(command.c_str());
}
