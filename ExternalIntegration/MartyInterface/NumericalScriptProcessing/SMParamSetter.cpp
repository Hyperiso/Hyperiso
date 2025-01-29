#include "SMParamSetter.h"
#include <cmath>
#include <set>
#include <string>

void SMParamSetter::setParam(const std::string& name, const Interpreter::InterpretedParam& interpretedParam) {
    double value = calculateValue(name, interpretedParam);

    if (interpretedParam.complex == true) {
        std::string im_block = "IM" + interpretedParam.block.substr(2, interpretedParam.block.size());
        if (interpretedParam.is_bsm) {
            ParameterType type = ParameterTypeMapper::enum_elt(ModelMapper::str(MemoryManager::GetInstance()->getModel()));
            params[name+"_img"] = (*Parameters::GetInstance(type))(interpretedParam.block, interpretedParam.code);
            params[name+"_rel"] = (*Parameters::GetInstance(type))(im_block, interpretedParam.code);
        } else {
            params[name+"_img"] = (*Parameters::GetInstance())(interpretedParam.block, interpretedParam.code);
            params[name+"_rel"] = (*Parameters::GetInstance())(im_block, interpretedParam.code);
        }
    } else {
        std::set<std::string> special = {"KIN", "WEIN", "Finite", "S2_THETAW", "REGPROP"};

        if (special.find(interpretedParam.block) != special.end()) {
            params[name] = value;
        } else {
            if (interpretedParam.is_bsm) {
                ParameterType type = ParameterTypeMapper::enum_elt(ModelMapper::str(MemoryManager::GetInstance()->getModel()));
                params[name] = (*Parameters::GetInstance(type))(interpretedParam.block, interpretedParam.code);
            } else {
                params[name] = (*Parameters::GetInstance())(interpretedParam.block, interpretedParam.code);
            }
            // params[name] = jsonparser->getElement(interpretedParam.block, interpretedParam.code);
        }
    }
}

double SMParamSetter::calculateValue(const std::string& name, const Interpreter::InterpretedParam& interpretedParam) {
    if (interpretedParam.block == "KIN") {
        if (interpretedParam.code == 34) {
            return -std::pow((*Parameters::GetInstance())("MASS", 13),2.);
        }
        return (std::pow(QCDHelper::mass_b_msbar(),2.) + std::pow((*Parameters::GetInstance())("MASS", 3), 2.))/2;
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
