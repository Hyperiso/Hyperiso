#include "SMParamSetter.h"

void SMParamSetter::setParam(const std::string& name, const Interpreter::InterpretedParam& interpretedParam) {
    LOG_INFO("setting parameter", name, interpretedParam.block, interpretedParam.code);
    std::set<std::string> special = {"KIN", "WEIN", "Finite", "REGPROP"};

    if (special.find(interpretedParam.block) != special.end()) {
        params[name] = calculateValue(name, interpretedParam);
    } else if (interpretedParam.block == "MASS" && (interpretedParam.code == LhaID(5) || interpretedParam.code == LhaID(6))) {
        if (interpretedParam.code == LhaID(5)) {
            params[name] = QCDHelper::mass_b_msbar();
        } else {
            params[name] = QCDHelper::mass_t_msbar();
        }
    } else {
        if (interpretedParam.is_bsm) {
            ParameterType type = ParameterTypeMapper::enum_elt(ModelMapper::str(ModelAPI().get()));
            params[name] = (*Parameters::GetInstance(type))(interpretedParam.block, interpretedParam.code);
        } else {
            params[name] = (*Parameters::GetInstance())(interpretedParam.block, interpretedParam.code);
        }
    }
}

double SMParamSetter::calculateValue(const std::string& name, const Interpreter::InterpretedParam& interpretedParam) {
    if (interpretedParam.block == "KIN") {
        if (interpretedParam.code == LhaID(34)) {
            return -std::pow((*Parameters::GetInstance())("MASS", 13),2.);
        }
        return (std::pow(QCDHelper::mass_b_msbar(),2.) + std::pow((*Parameters::GetInstance())("MASS", 3), 2.))/2;
    }
    if (interpretedParam.block == "WEIN") {
        return 0.5;
    }
    if (interpretedParam.block == "REGPROP") {
        return 1e-3;
    }

    return 1.0;
}
