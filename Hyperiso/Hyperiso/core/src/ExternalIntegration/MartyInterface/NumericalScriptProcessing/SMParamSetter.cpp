#include "SMParamSetter.h"

void SMParamSetter::setParam(const std::string& name, const Interpreter::InterpretedParam& interpretedParam) {
    LOG_INFO("setting parameter", name, interpretedParam.block, interpretedParam.code);
    std::set<std::string> special = {"KIN", "WEIN", "Finite", "REGPROP"}; //TODO : put else where, needed in wilson marty

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
            std::cout << "ohoh : " << name << std::endl;
            ParameterType type = ParameterTypeMapper::enum_elt(ModelMapper::str(ModelAPI().get()));
            std::cout << (*Parameters::GetInstance(type))(interpretedParam.block, interpretedParam.code) << std::endl;
            if (interpretedParam.is_complex) {
                params[name+ "_rel"] = (*Parameters::GetInstance(type))(interpretedParam.block, interpretedParam.code).real();
                params[name + "_img"] = (*Parameters::GetInstance(type))(interpretedParam.block, interpretedParam.code).imag();
            } else {
                params[name] = (*Parameters::GetInstance(type))(interpretedParam.block, interpretedParam.code);
            }
        } else {
            std::cout << "eheh : " << name << std::endl;
            if (interpretedParam.is_complex) {
                params[name+ "_rel"] = (*Parameters::GetInstance())(interpretedParam.block, interpretedParam.code).real();
                params[name + "_img"] = (*Parameters::GetInstance())(interpretedParam.block, interpretedParam.code).imag();
            } else {
                params[name] = (*Parameters::GetInstance())(interpretedParam.block, interpretedParam.code);
            }
            // params[name] = (*Parameters::GetInstance())(interpretedParam.block, interpretedParam.code);
        }
    }
}

scalar_t SMParamSetter::calculateValue(const std::string& name, const Interpreter::InterpretedParam& interpretedParam) {
    if (interpretedParam.block == "KIN") {
        if (interpretedParam.code == LhaID(34)) {
            return -pow((*Parameters::GetInstance())("MASS", 13),2.);
        }
        return (pow(QCDHelper::mass_b_msbar(),2.) + std::pow((*Parameters::GetInstance())("MASS", 3), 2.))/2.;
    }
    if (interpretedParam.block == "WEIN") {
        return 0.5;
    }
    if (interpretedParam.block == "REGPROP") {
        return 1e-3;
    }

    return 1.0;
}
