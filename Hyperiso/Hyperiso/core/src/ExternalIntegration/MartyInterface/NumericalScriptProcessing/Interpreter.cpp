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
        auto it = defaultDatabase->getParams().find(param.name);
        if (it != defaultDatabase->getParams().end()) {
            interpreted.block = it->second.block;
            interpreted.code = it->second.pdgCode;
            interpreted.is_complex = param.complex;
            interpreted.is_bsm = false;
        } else {
            it = modelDatabase->getParams().find(param.name);
            if (it != modelDatabase->getParams().end()) {
                if (ModelAPI().get() == Model::SM) {
                    LOG_ERROR("LogicError", "Trying to access BSM parameter in SM calculation. Check MARTY parameter mapping files.");
                }
                interpreted.block = it->second.block;
                interpreted.code = it->second.pdgCode;
                interpreted.is_complex = param.complex;
                interpreted.is_bsm = true;
            } else {
                std::cerr << "Error: Parameter " << param.name << " not found in model or SM mapping databases." << std::endl;
                std::runtime_error("");
                continue;
            }
        }

        interpretedParams[param.name] = interpreted;
    }

    return interpretedParams;
}