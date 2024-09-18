#include "SMParamSetter.h"

void SMParamSetter::setParam(const std::string& name, const Interpreter::InterpretedParam& interpretedParam) {
    double value = calculateValue(name, interpretedParam);

    // Exemple de logique spécifique : pour certains blocs, on remplace directement par (*sm)(block, value)
    if (interpretedParam.block == "KIN" || interpretedParam.block == "WEIN" || interpretedParam.block == "Finite") {
        std::cout << interpretedParam.block << std::endl;
        params[name] = value;
    } else {
        std::cout << interpretedParam.block << std::endl;
        params[name] = jsonparser->getElement(interpretedParam.block, interpretedParam.code);
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
