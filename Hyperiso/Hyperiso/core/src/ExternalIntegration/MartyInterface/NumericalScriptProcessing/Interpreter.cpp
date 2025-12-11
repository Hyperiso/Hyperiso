#include "Interpreter.h"
#include <iostream>

Interpreter::Interpreter(const std::string& model,
                        std::shared_ptr<ICoreAPI<Model>> api,
                        std::shared_ptr<IInterpreterPortsFactory> ports)
: marty_api(std::move(api))
{
    const auto modelPath = FileNameManager::getInstance("C1", model)->getjsondbmodel();
    const auto smPath    = FileNameManager::getInstance("C1", "SM")->getjsondbmodel();
    resolver = ports->makeResolver(model, modelPath, smPath);
}

std::unordered_map<std::string, InterpretedParam>
Interpreter::interpret(std::vector<Extractor::Parameter>& params) {
    const bool modelIsSM = (marty_api->get() == Model::SM);

    auto resolved = resolver->resolve(params, modelIsSM);

    std::unordered_map<std::string, InterpretedParam> out;
    out.reserve(resolved.size());
    for (auto& [name, r] : resolved) {
        out[name] = { r.block, r.code, r.is_bsm, r.is_complex };
    }
    return out;
}

Interpreter::Interpreter(const Interpreter& other)
    : resolver(other.resolver ? other.resolver->clone() : nullptr),
      marty_api(other.marty_api)
{}

Interpreter& Interpreter::operator=(const Interpreter& other) {
    if (this == &other) return *this;
    std::unique_ptr<IParameterResolver> newResolver =
        other.resolver ? other.resolver->clone() : nullptr;
    marty_api = other.marty_api;
    resolver  = std::move(newResolver);
    return *this;
}