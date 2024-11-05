#include "GppCompilerStrategy.h"
#include <cstdlib>
#include "FileNameManager.h"

void GppCompilerStrategy::compile_run(const std::string& sourceFile, const std::string& outputBinary) {
    if (!this->check_if_compile(outputBinary)) {
        this->compile(sourceFile, outputBinary);
    }
    std::string command_run = "cd " + FileNameManager::getInstance(wilson, model)->getOutputDir() +" && "+outputBinary;
    system(command_run.c_str());
}

void GppCompilerStrategy::compile(const std::string& sourceFile, const std::string& outputBinary) {
    std::string root_path = project_root.data();
    std::string command_compile = "g++ -o " + outputBinary + " " + sourceFile +" -L"+
    root_path+ "/ExternalIntegration/MARTY/MARTY_INSTALL/lib -Wl,-rpath," 
    + root_path + "/ExternalIntegration/MARTY/MARTY_INSTALL/lib" + " -lmarty"  + " -lgfortran";
    system(command_compile.c_str());
}