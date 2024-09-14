#ifndef INTERPRETER_H
#define INTERPRETER_H

#include "Extractor.h"
#include <string>
#include <unordered_map>

class Interpreter {
public:
    struct InterpretedParam {
        std::string block;
        int code;
    };

    std::unordered_map<std::string, InterpretedParam> interpret(const std::vector<Extractor::Parameter>& params);

private:
    int getCode(const std::string& name);
};

#endif // INTERPRETER_H
