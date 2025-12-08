#ifndef SMPARAMSETTER_H
#define SMPARAMSETTER_H

#include "config.hpp"
#include "LhaID.h"
#include "Interpreter.h"
#include <cmath>
#include <set>
#include <string>
#include <iostream>
#include <cmath>
#include "IMartyParameterProxy.h"

class SMParamSetter {
public:
    SMParamSetter(const std::string& model, std::set<std::string> special_blocks, std::shared_ptr<IMartyParameterProxy<std::string, LhaID>> sm_proxy, std::shared_ptr<IMartyParameterProxy<std::string, LhaID>> bsm_proxy = nullptr) : special_blocks(special_blocks) {
        std::string root_path = project_assets_root.data();
        this->sm_proxy = sm_proxy;

        if (model == "SM") {
            this->model_type = Model::SM;
        } else if (model == "THDM") {
            this->model_type = Model::THDM;
        } else if (model == "MSSM" || model == "NMSSM") {
            this->model_type = Model::SUSY;
        } else {
            this->model_type = Model::MARTY;
        }
        
        if (model_type != Model::SM) {
            this->bsm_proxy = bsm_proxy;
        }
    }

    std::unordered_map<std::string, double> setParam(const std::string& name, const InterpretedParam& interpretedParam);

private:

    scalar_t calculateValue(const std::string& name, const InterpretedParam& interpretedParam);

    Model model_type;
    std::set<std::string> special_blocks;
    std::shared_ptr<IMartyParameterProxy<std::string, LhaID>> sm_proxy;
    std::shared_ptr<IMartyParameterProxy<std::string, LhaID>> bsm_proxy;
};

#endif // SMPARAMSETTER_H
