#ifndef INTERPRETER_H
#define INTERPRETER_H

#include "Extractor.h"
#include "General.h"
#include <string>
#include <unordered_map>
#include <memory>
#include "ICoreAPI.h"
#include "IInterpreterPortsFactory.h" 
#include "FileNameManager.h"
#include "InterpretedParam.h"

class Interpreter {
public:
    

    Interpreter(const std::string& model,
                            std::shared_ptr<ICoreAPI<Model>> api,
                            std::shared_ptr<IInterpreterPortsFactory> ports)
    : marty_api(std::move(api))
    {
        const auto modelPath = FileNameManager::getInstance("C1", model)->getjsondbmodel();
        const auto smPath    = FileNameManager::getInstance("C1", "SM")->getjsondbmodel();
        resolver = ports->makeResolver(model, modelPath, smPath);
    }

    std::unordered_map<std::string, InterpretedParam>
    interpret(std::vector<Extractor::Parameter>& params);

    Interpreter(const Interpreter& other);
    Interpreter& operator=(const Interpreter& other);
    Interpreter(Interpreter&&) noexcept = default;
    Interpreter& operator=(Interpreter&&) noexcept = default;

private:
    std::unique_ptr<IParameterResolver> resolver;
    std::shared_ptr<ICoreAPI<Model>> marty_api;
};


#endif // INTERPRETER_H