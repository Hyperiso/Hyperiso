#include "SMParamSetter.h"

void SMParamSetter::setParam(const std::string& name, const Interpreter::InterpretedParam& interpretedParam) {
    LOG_INFO("setting parameter", name, interpretedParam.block, interpretedParam.code);
    std::set<std::string> special = {"KIN", "WEIN", "Finite", "REGPROP"}; //TODO : put else where, needed in wilson marty
    if (special.find(interpretedParam.block) != special.end()) {
        params[name] = calculateValue(name, interpretedParam);
    } else if (interpretedParam.block == "MASS" && (interpretedParam.code == LhaID(5) || interpretedParam.code == LhaID(6))) {
        if (interpretedParam.code == LhaID(5)) {
            params[name] = sm_proxy("QCD", LhaID(5, 1));
        } else {
            params[name] = sm_proxy("QCD", 6);
        }
    } else {
        if (interpretedParam.is_bsm) {
            if (interpretedParam.is_complex) {
                params[name+ "_rel"] = bsm_proxy(interpretedParam.block, interpretedParam.code).real();
                params[name + "_img"] = bsm_proxy(interpretedParam.block, interpretedParam.code).imag();
            } else {
                params[name] = bsm_proxy(interpretedParam.block, interpretedParam.code);
            }
        } else {
            if (interpretedParam.is_complex) {
                params[name+ "_rel"] = sm_proxy(interpretedParam.block, interpretedParam.code).real();
                params[name + "_img"] = sm_proxy(interpretedParam.block, interpretedParam.code).imag();
            } else {
                params[name] = sm_proxy(interpretedParam.block, interpretedParam.code);
            }
        }
    }
}

scalar_t SMParamSetter::calculateValue(const std::string& name, const Interpreter::InterpretedParam& interpretedParam) {
    if (interpretedParam.block == "KIN") {
        if (interpretedParam.code == LhaID(34)) {
            return -pow(sm_proxy("MASS", 13),2.);
        }
        return (pow(sm_proxy("QCD", LhaID(5, 1)) ,2.) + std::pow(sm_proxy("MASS", 3), 2.))/2.;
    }
    if (interpretedParam.block == "WEIN") {
        return 0.5;
    }
    if (interpretedParam.block == "REGPROP") {
        return 1e-3;
    }

    return 1.0;
}
