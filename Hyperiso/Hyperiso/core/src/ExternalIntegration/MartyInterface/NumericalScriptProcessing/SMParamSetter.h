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

class SMParamSetter {
public:
    SMParamSetter(std::unordered_map<std::string, double>& params, const std::string& model) : params(params) {
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
        // if (model == "MSSM" || model == "NMSSM") {
        //     jsonparser->loadFromFile(root_path + "savestate/parameters_SUSY.json");
        // } else if (model == "SM" || model == "THDM") {
        //     jsonparser->loadFromFile(root_path + "savestate/parameters_" + model + ".json");
        // } else {
        //     jsonparser->loadFromFile(root_path + "savestate/parameters_GENERAL.json");
        // }
    }

    void setParam(const std::string& name, const Interpreter::InterpretedParam& interpretedParam);

private:
    std::unordered_map<std::string, double>& params;

    double calculateValue(const std::string& name, const Interpreter::InterpretedParam& interpretedParam);

    // JSONParser *jsonparser = JSONParser::getInstance(0);
    Model model_type;
};

#endif // SMPARAMSETTER_H
