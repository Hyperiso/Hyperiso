#ifndef SMPARAMSETTER_H
#define SMPARAMSETTER_H

#include "IParamSetter.h"
#include <iostream>

class SMParamSetter : public IParamSetter {
public:
    SMParamSetter(std::unordered_map<std::string, double>& params) : params(params) {}

    void setParam(const std::string& name, const Interpreter::InterpretedParam& interpretedParam) override;

private:
    std::unordered_map<std::string, double>& params;

    double calculateValue(const std::string& name, const Interpreter::InterpretedParam& interpretedParam);
};

#endif // SMPARAMSETTER_H
