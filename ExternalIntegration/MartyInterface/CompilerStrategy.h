#ifndef COMPILER_STRATEGY_H
#define COMPILER_STRATEGY_H

#include <string>

class CompilerStrategy {
public:
    virtual ~CompilerStrategy() = default;
    virtual void compile_run(const std::string& sourceFile, const std::string& outputBinary) = 0;
};

#endif // COMPILER_STRATEGY_H
