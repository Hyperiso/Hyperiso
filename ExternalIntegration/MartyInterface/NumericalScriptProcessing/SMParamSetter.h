#ifndef SMPARAMSETTER_H
#define SMPARAMSETTER_H

#include "IParamSetter.h"
#include "JsonParameters.h"
#include "config.hpp"

#include <iostream>
#include <cmath>
class SMParamSetter : public IParamSetter {
public:
    SMParamSetter(std::unordered_map<std::string, double>& params, const std::string& model) : params(params) {
        std::string root_path = project_root.data();
        this->model = model;
        if (this->model == "MSSM" || this->model == "NMSSM") {
            jsonparser->loadFromFile(root_path+ "/DataBase/Params/data_SUSY.json");
        } else if (this->model == "SM" || this->model == "THDM") {
            jsonparser->loadFromFile(root_path+ "/DataBase/Params/data_" + this->model + ".json");
        } else {
            jsonparser->loadFromFile(root_path+ "/DataBase/Params/data_GENERAL.json");
        }
    }

    void setParam(const std::string& name, const Interpreter::InterpretedParam& interpretedParam) override;

private:
    std::unordered_map<std::string, double>& params;

    double calculateValue(const std::string& name, const Interpreter::InterpretedParam& interpretedParam);

    JSONParser *jsonparser = JSONParser::getInstance(0);
    std::string model;
};

#endif // SMPARAMSETTER_H
