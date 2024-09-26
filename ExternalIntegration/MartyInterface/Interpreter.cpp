#include "Interpreter.h"

std::unordered_map<std::string, Interpreter::InterpretedParam> Interpreter::interpret(std::vector<Extractor::Parameter>& params) {
    std::unordered_map<std::string, InterpretedParam> interpretedParams;

    for (auto& param : params) {
        InterpretedParam interpreted;
        
        if (param.name.rfind("m_", 0) == 0 || param.name.rfind("M_", 0) == 0) {
            interpreted.block = "MASS";
            interpreted.code = getCode(param.name);
        } else if (param.name.rfind("g_", 0) == 0) {
            interpreted.block = "GAUGE";
            interpreted.code = getCode(param.name);
        } else if (param.name.rfind("e_", 0) == 0) {
            interpreted.block = "SMINPUT";
            interpreted.code = getCode(param.name);
        } else if (param.name.rfind("G_F", 0) == 0) {
            interpreted.block = "SMINPUT";
            interpreted.code = 2;
        } else if (param.name.rfind("s_12", 0) == 0) {
            interpreted.block = "KIN";
            interpreted.code = getCode(param.name);
        } else if (param.name.rfind("V_", 0) == 0) {
            interpreted.block = "RECKM";
            interpreted.code = getCode(param.name);
            if (param.complex == true) {
                InterpretedParam interpreted_2;
                interpreted_2.block = "IMCKM";
                interpreted_2.code = interpreted.code;
                interpretedParams[param.name+"_im"] = interpreted_2;
                param.name = param.name+"_re";

            }
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
    if (name.find("_tb") != std::string::npos) return 22;
    if (name.find("_ts") != std::string::npos) return 21;
    if (name.find("us") != std::string::npos) return 01;
    if (name.find("cb") != std::string::npos) return 12;
    if (name.find("cs") != std::string::npos) return 11;
    if (name.find("ud") != std::string::npos) return 00;
    if (name.find("ub") != std::string::npos) return 02;
    if (name.find("_u") != std::string::npos) return 2;
    if (name.find("_d") != std::string::npos) return 1;
    if (name.find("_s") != std::string::npos) return 3;
    if (name.find("_c") != std::string::npos) return 4;
    if (name.find("_t") != std::string::npos) return 6;
    if (name.find("_b") != std::string::npos) return 5;
    if (name.find("_W") != std::string::npos) return 24;
    if (name.find("em") != std::string::npos) return 10;
    if (name.find("s_12") != std::string::npos) return 4;
    if (name.find("theta_W") != std::string::npos) return 42;
    return -1; // Default case if no code is found
}
