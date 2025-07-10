#ifndef SMPARAMSETTER_H
#define SMPARAMSETTER_H

#include "config.hpp"
#include "General.h"
#include "Parameters.h"
#include "ModelAPI.h"
#include "QCDHelper.h"
#include "Interpreter.h"
#include <cmath>
#include <set>
#include <string>
#include <iostream>
#include <cmath>
#include "MartyParameterProxy.h"

class SMParamSetter {
public:
    SMParamSetter(std::unordered_map<std::string, double>& params, const std::string& model, std::set<std::string> special_blocks) :  params(params), special_blocks(special_blocks) {
        std::string root_path = project_assets_root.data();

        if (model == "SM") {
            this->model_type = Model::SM;
        } else if (model == "THDM") {
            this->model_type = Model::THDM;
        } else if (model == "MSSM" || model == "NMSSM") {
            this->model_type = Model::SUSY;
        } else {
            this->model_type = Model::CUSTOM;
        }
        
        if (model_type != Model::SM) {
            bsm_proxy = MartyParameterProxy(ParameterType::BSM);
        }
    }

    void setParam(const std::string& name, const Interpreter::InterpretedParam& interpretedParam);

private:
    std::unordered_map<std::string, double>& params;

    scalar_t calculateValue(const std::string& name, const Interpreter::InterpretedParam& interpretedParam);

    Model model_type;
    std::set<std::string> special_blocks;
    MartyParameterProxy sm_proxy {ParameterType::SM};
    std::optional<MartyParameterProxy> bsm_proxy;
};

#endif // SMPARAMSETTER_H
