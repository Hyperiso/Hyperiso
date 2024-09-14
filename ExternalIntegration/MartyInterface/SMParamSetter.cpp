#include "SMParamSetter.h"

void SMParamSetter::setParam(const std::string& name, const Interpreter::InterpretedParam& interpretedParam) {
    double value = calculateValue(name, interpretedParam);

    // Exemple de logique spécifique : pour certains blocs, on remplace directement par (*sm)(block, value)
    if (interpretedParam.block == "MASS" || interpretedParam.block == "SMINPUT") {
        std::cout << "(*sm)(" << interpretedParam.block << ", " << value << ");\n";
    } else {
        params[name] = value;
    }
}

double SMParamSetter::calculateValue(const std::string& name, const Interpreter::InterpretedParam& interpretedParam) {
    // Exemple de logique pour les blocs KIN et WEIN
    if (interpretedParam.block == "KIN") {
        return 18;
    }
    if (interpretedParam.block == "WEIN") {
        return 21;
    }

    // Valeur par défaut
    return 1.0; // Exemple de valeur par défaut
}
