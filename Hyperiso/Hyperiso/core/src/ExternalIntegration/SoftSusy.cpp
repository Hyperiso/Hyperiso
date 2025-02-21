#ifdef BUILD_WITH_SOFTSUSY
#include <iostream>
#include <fstream>
#include <filesystem>
#include <cstdlib> // For system()
#include "Logger.h"
#include "SoftSusy.h"
#include <sstream>

namespace fs = std::filesystem;

/**
 * @brief Executes the spectrum calculation using SOFTSUSY.
 * 
 * @param inputFilePath Path to the input file containing model parameters.
 * @param outputFilePath Path where the calculated spectrum will be written.
 */
void SoftsusyCalculator::calculateSpectrum(const std::string& inputFilePath, const std::string& outputFilePath) {

    std::string root_tp_file = project_tp_root.data();
    std::string root_assets_file = project_assets_root.data();
    // Example system call to SOFTSUSY - replace with actual implementation
    std::string command;
    if (inputFilePath.starts_with("/")) {
        command = root_tp_file + "SOFTSUSY/src/SOFTSUSY/softpoint.x leshouches < " + inputFilePath + " > " + outputFilePath;
    } else {
        command = root_tp_file + "SOFTSUSY/src/SOFTSUSY/softpoint.x leshouches < " + root_assets_file + inputFilePath + " > " + outputFilePath;
    }

    LOG_DEBUG("SOFTSUSY COMMAND : " + command);

    int result = system(command.c_str());
    if (result != 0) {
        LOG_ERROR("SoftwareError", "SOFTSUSY execution failed with code " + std::to_string(result));
    } else {
        LOG_INFO("SOFTSUSY execution successful.");
    }
}

#endif