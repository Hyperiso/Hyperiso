#include "Interpreter.h"
#include <iostream>
#include "FileNameManager.h"

Interpreter::Interpreter(const std::string& model) {
    modelDatabase = MappingDatabase::getInstance(model, FileNameManager::getInstance("C1", model)->getjsondbmodel());
    defaultDatabase = MappingDatabase::getInstance("SM", FileNameManager::getInstance("C1", "SM")->getjsondbmodel());
}

std::unordered_map<std::string, Interpreter::InterpretedParam> Interpreter::interpret(std::vector<Extractor::Parameter>& params) {
    std::unordered_map<std::string, InterpretedParam> interpretedParams;

    for (auto& param : params) {
        InterpretedParam interpreted;
        std::cout << "Processing parameter: " << param.name << std::endl;

        // Cherche d'abord dans la base de données du modèle spécifié
        auto it = modelDatabase->getParams().find(param.name);
        if (it != modelDatabase->getParams().end()) {
            // Si trouvé dans la base de données du modèle spécifié
            interpreted.block = it->second.block;
            interpreted.code = it->second.pdgCode;
        } else {
            // Sinon, on cherche dans la base de données SM
            it = defaultDatabase->getParams().find(param.name);
            if (it != defaultDatabase->getParams().end()) {
                // Si trouvé dans la base de données SM
                interpreted.block = it->second.block;
                interpreted.code = it->second.pdgCode;
            } else {
                // Si non trouvé dans aucune base, afficher une erreur
                std::cerr << "Error: Parameter " << param.name << " not found in model or SM databases." << std::endl;
                continue;
            }
        }

        // Enregistre le paramètre interprété dans le résultat final
        interpretedParams[param.name] = interpreted;
    }

    return interpretedParams;
}