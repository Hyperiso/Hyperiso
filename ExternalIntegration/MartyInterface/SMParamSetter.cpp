#include "SMParamSetter.h"

void SMParamSetter::setParam(const std::string& name, const Interpreter::InterpretedParam& interpretedParam) {
    double value = calculateValue(name, interpretedParam);

    if (interpretedParam.block == "KIN" || interpretedParam.block == "WEIN" || interpretedParam.block == "Finite") {
        std::cout << interpretedParam.block << std::endl;
        params[name] = value;
    } else {
        std::cout << interpretedParam.block << std::endl;
        params[name] = jsonparser->getElement(interpretedParam.block, interpretedParam.code);
    }
}

double SMParamSetter::calculateValue(const std::string& name, const Interpreter::InterpretedParam& interpretedParam) {
    if (interpretedParam.block == "KIN") {
        return 18;
    }
    if (interpretedParam.block == "WEIN") {
        return 21;
    }

    return 1.0;
}
