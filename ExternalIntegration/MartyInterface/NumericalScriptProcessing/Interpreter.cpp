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

        auto it = modelDatabase->getParams().find(param.name);
        if (it != modelDatabase->getParams().end()) {
            interpreted.block = it->second.block;
            interpreted.code = it->second.pdgCode;
        } else {
            it = defaultDatabase->getParams().find(param.name);
            if (it != defaultDatabase->getParams().end()) {
                interpreted.block = it->second.block;
                interpreted.code = it->second.pdgCode;
            } else {
                std::cerr << "Error: Parameter " << param.name << " not found in model or SM databases." << std::endl;
                std::runtime_error("");
                continue;
            }
        }

        interpretedParams[param.name] = interpreted;
    }

    return interpretedParams;
}