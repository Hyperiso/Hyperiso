#include "Interpreter.h"

std::unordered_map<std::string, Interpreter::InterpretedParam> Interpreter::interpret(const std::vector<Extractor::Parameter>& params) {
    std::unordered_map<std::string, InterpretedParam> interpretedParams;

    for (const auto& param : params) {
        InterpretedParam interpreted;

        if (param.name.rfind("m_", 0) == 0) {
            interpreted.block = "MASS";
            interpreted.code = getCode(param.name);
        } else if (param.name.rfind("g_", 0) == 0) {
            interpreted.block = "GAUGE";
            interpreted.code = getCode(param.name);
        } else if (param.name.rfind("e_", 0) == 0) {
            interpreted.block = "SMINPUT";
            interpreted.code = getCode(param.name);
        } else if (param.name.rfind("s_12", 0) == 0) {
            interpreted.block = "KIN";
            interpreted.code = getCode(param.name);
        } else if (param.name.rfind("theta_W", 0) == 0) {
            interpreted.block = "WEIN";
            interpreted.code = getCode(param.name);
        } else if (param.name.rfind("Finite", 0) == 0) {
            interpreted.block = "Finite";
            interpreted.code = getCode(param.name);
        }

        interpretedParams[param.name] = interpreted;
    }

    return interpretedParams;
}

int Interpreter::getCode(const std::string& name) {
    if (name.find("_s") != std::string::npos) return 3;
    if (name.find("_t") != std::string::npos) return 6;
    if (name.find("_b") != std::string::npos) return 5;
    if (name.find("em") != std::string::npos) return 2;
    if (name.find("s_12") != std::string::npos) return 4;
    if (name.find("theta_W") != std::string::npos) return 42;
    return -1; // Default case if no code is found
}
