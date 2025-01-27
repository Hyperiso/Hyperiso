#ifndef INTERPRETER_H
#define INTERPRETER_H

#include "Extractor.h"
#include <string>
#include <unordered_map>
#include "MappingDataBase.h"

class Interpreter {
public:
    struct InterpretedParam {
        std::string block;
        bool complex;
        int code;
    };
    Interpreter(const std::string& model = "SM");
    std::unordered_map<std::string, InterpretedParam> interpret(std::vector<Extractor::Parameter>& params);

private:
    std::shared_ptr<MappingDatabase> modelDatabase;
    std::shared_ptr<MappingDatabase> defaultDatabase;
};

#endif // INTERPRETER_H
