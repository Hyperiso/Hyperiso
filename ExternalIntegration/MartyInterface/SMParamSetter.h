#ifndef SMPARAMSETTER_H
#define SMPARAMSETTER_H

#include "IParamSetter.h"
#include "JsonParameters.h"
#include "config.hpp"

#include <iostream>

class SMParamSetter : public IParamSetter {
public:
    SMParamSetter(std::unordered_map<std::string, double>& params) : params(params) {
        std::string root_path = project_root.data();
        jsonparser->loadFromFile(root_path+ "/DataBase/Params/data_SM.json");
        }

    void setParam(const std::string& name, const Interpreter::InterpretedParam& interpretedParam) override;

private:
    std::unordered_map<std::string, double>& params;

    double calculateValue(const std::string& name, const Interpreter::InterpretedParam& interpretedParam);

    JSONParser *jsonparser = JSONParser::getInstance(0);
};

#endif // SMPARAMSETTER_H
