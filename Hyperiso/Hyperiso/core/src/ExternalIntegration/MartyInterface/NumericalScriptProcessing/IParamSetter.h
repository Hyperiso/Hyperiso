#ifndef IPARAMSETTER_H
#define IPARAMSETTER_H

#include <string>
#include <unordered_map>
#include "Interpreter.h"

class IParamSetter {
public:
    virtual void setParam(const std::string& name, const Interpreter::InterpretedParam& interpretedParam) = 0;
    virtual ~IParamSetter() = default;
};

#endif // IPARAMSETTER_H
