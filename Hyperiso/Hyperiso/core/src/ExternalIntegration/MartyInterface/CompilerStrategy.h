#ifndef COMPILER_STRATEGY_H
#define COMPILER_STRATEGY_H

#include <string>
#include <fstream>
#include <sys/stat.h>
#include <iostream>
#include <cstdio>
#include <memory>
#include <array>
#include "Include.h"

bool executeCommand(const std::string& command);

class CompilerStrategy {
public:
    CompilerStrategy(std::string model, std::string wilson) : model(model), wilson(wilson) {}
    virtual ~CompilerStrategy() = default;
    virtual void compile_run(const std::string& sourceFile, const std::string& outputBinary) = 0;
    virtual void compile(const std::string& sourceFile, const std::string& outputBinary) = 0;
    virtual bool check_if_compile(const std::string& outputBinary);

protected:
    std::string wilson{};
    std::string model{};
};

#endif
