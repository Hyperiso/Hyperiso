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
 * @brief Extracts the project path from the JSON configuration file.
 * 
 * This function reads the JSON file specified by `configFile` and searches for the
 * "project_root" key to extract the associated value. The absolute path of the project
 * is returned as a string.
 * 
 * @param[in] configFile The path to the configuration JSON file.
 * @return The absolute path of the project.
 */
std::string getProjectRootFromConfig() {

    std::string configFile = "../config.json";

    std::ifstream ifs(configFile);
    if (!ifs.is_open()) {
        LOG_ERROR("FileError", "Not possible to open the config file");
        return "";
    }

    std::string line;
    std::string projectRoot;
    while (std::getline(ifs, line)) {

        size_t pos = line.find("\"project_root\"");
        if (pos != std::string::npos) {

            pos = line.find_first_of('"', pos + 1); 
            size_t endPos = line.find_last_of('"'); 
            projectRoot = line.substr(pos + 1, endPos - pos - 1);
            break;
        }
    }

    // Fermer le fichier
    ifs.close();
    LOG_INFO("Project root folder is " +projectRoot);
    // Retourner la valeur du projet root
    return projectRoot;
}

/**
 * Creates a new calculator instance.
 * This method implements the Factory Method pattern, providing a way to instantiate
 * ICalculator objects. It currently defaults to creating instances of SoftsusyCalculator.
 * 
 * @return ICalculator* Pointer to the newly created ICalculator instance.
 */
// ICalculator* SoftsusyCalculatorFactory::createCalculator() {
//     return new SoftsusyCalculator();
// }



/**
 * @brief Executes the spectrum calculation using SOFTSUSY.
 * 
 * @param inputFilePath Path to the input file containing model parameters.
 * @param outputFilePath Path where the calculated spectrum will be written.
 */
void SoftsusyCalculator::calculateSpectrum(const std::string& inputFilePath, const std::string& outputFilePath) {

    std::string root_file = project_root.data();
    // Example system call to SOFTSUSY - replace with actual implementation
    std::string command;
    if (inputFilePath.starts_with("/")) {
        command = root_file + "/ExternalIntegration/SOFTSUSY/src/SOFTSUSY/softpoint.x leshouches < " + inputFilePath + " > " + outputFilePath;
    } else {
        command = root_file + "/ExternalIntegration/SOFTSUSY/src/SOFTSUSY/softpoint.x leshouches < " + root_file + "/" + inputFilePath + " > " + outputFilePath;
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