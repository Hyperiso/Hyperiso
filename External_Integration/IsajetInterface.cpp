#include "IsajetInterface.h"
#include <iostream>
#include <fstream>

IsajetInterface::IsajetInterface(const std::string& executablePath)
    : isajetExecutable(executablePath) {}

ModelOutput IsajetInterface::runIsajet(const ModelParameters& params) {
    // On écrit les paramètres dans un fichier de configuration pour ISAJET
    writeConfigFile(params);

    // On exécute ISAJET avec le fichier de configuration
    std::string command = isajetExecutable + " < config.txt";
    system(command.c_str());

    // On lit et analyse les résultats d'ISAJET
    ModelOutput output = parseOutputFile("isajet_output.txt");
    return output;
}

void IsajetInterface::writeConfigFile(const ModelParameters& params) {
    std::ofstream configFile("config.txt");
    if (configFile.is_open()) {
        //  paramètres dans le fichier de configuration
        configFile << "Paramètre1 = " << params.getParam1() << std::endl;
        configFile << "Paramètre2 = " << params.getParam2() << std::endl;
        // ... autres paramètres
        configFile.close();
    }
}

ModelOutput IsajetInterface::parseOutputFile(const std::string& filename) {
    ModelOutput output;
    std::ifstream outputFile(filename);
    if (outputFile.is_open()) {
        // Analyse du fichier de sortie
        // ...
        outputFile.close();
    }
    return output;
}
