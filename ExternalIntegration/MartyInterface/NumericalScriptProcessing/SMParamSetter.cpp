#include "SMParamSetter.h"
#include <cmath>
#include <set>
#include <string>

void SMParamSetter::setParam(const std::string& name, const Interpreter::InterpretedParam& interpretedParam) {
    double value = calculateValue(name, interpretedParam);

    std::set<std::string> special = {"KIN", "WEIN", "Finite", "S2_THETAW", "REGPROP"};

    if (special.find(interpretedParam.block) != special.end()) {
        params[name] = value;
    } else {
        params[name] = jsonparser->getElement(interpretedParam.block, interpretedParam.code);
    }
}

double SMParamSetter::calculateValue(const std::string& name, const Interpreter::InterpretedParam& interpretedParam) {
    if (interpretedParam.block == "KIN") {
        if (interpretedParam.code == 34) {
            return -std::pow(jsonparser->getElement("MASS", 13),2.);
        }
        return (std::pow(jsonparser->getElement("MASS", 5),2.) + std::pow(jsonparser->getElement("MASS", 3), 2.))/2;
    }
    if (interpretedParam.block == "WEIN") {
        return 0.5;
    }
    if (interpretedParam.block == "S2_THETAW") {
        return pow(sin(0.5), 2);
    }
    if (interpretedParam.block == "REGPROP") {
        return 1e-3;
    }

    return 1.0;
}
