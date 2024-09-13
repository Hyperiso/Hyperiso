#include "MakeCompilerStrategy.h"
#include <cstdlib>

void MakeCompilerStrategy::compile_run(const std::string& sourceFile, const std::string& outputBinary) {
    std::string directory = "libs/C7_SM";
    
    std::string command = "cd " + directory + " && make";
    system(command.c_str());
    std::string command_run = "./" + directory + "/bin/" + "example_c7_sm.x";
    system(command_run.c_str());
}