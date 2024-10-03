#include "SMParamSetter.h"
#include <cmath>

void SMParamSetter::setParam(const std::string& name, const Interpreter::InterpretedParam& interpretedParam) {
    double value = calculateValue(name, interpretedParam);

    if (interpretedParam.block == "KIN" || interpretedParam.block == "WEIN" || interpretedParam.block == "Finite" || interpretedParam.block == "S2_THETAW") {
        params[name] = value;
    } else {
        params[name] = jsonparser->getElement(interpretedParam.block, interpretedParam.code);
    }
}

double SMParamSetter::calculateValue(const std::string& name, const Interpreter::InterpretedParam& interpretedParam) {
    if (interpretedParam.block == "KIN") {
        return (std::pow(jsonparser->getElement("MASS", 5),2.) + std::pow(jsonparser->getElement("MASS", 3), 2.))/2;
    }
    if (interpretedParam.block == "WEIN") {
        return 0.5;
    }
    if (interpretedParam.block == "S2_THETAW") {
        return pow(sin(0.5), 2);
    }

    return 1.0;
}
