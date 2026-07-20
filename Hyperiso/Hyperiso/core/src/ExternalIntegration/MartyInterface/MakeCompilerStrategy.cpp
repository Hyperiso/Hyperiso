#include "MakeCompilerStrategy.h"
#include "MartyRuntimeConfig.h"
#include <cstdlib>
#include <iostream>
#include <string>

void MakeCompilerStrategy::compile_run(const std::string& sourceFile, const std::string& outputBinary) {
    if (!this->check_if_compile(outputBinary)) {
        this->compile(sourceFile, outputBinary);
    }
    std::string command_run = MartyRuntimeConfig::shell_quote(outputBinary)
                            + " -Q " + std::to_string(this->Q_match);
    if (!param_file_.empty()) {
        command_run += " --param-file " + MartyRuntimeConfig::shell_quote(param_file_);
    }
    if (!output_file_.empty()) {
        command_run += " --output-file " + MartyRuntimeConfig::shell_quote(output_file_);
    }
    executeCommand(command_run);
}

void MakeCompilerStrategy::compile(const std::string& sourceFile, const std::string&) {
    std::string command = "cd " + MartyRuntimeConfig::shell_quote(sourceFile) + " && make";
    executeCommand(command);
}