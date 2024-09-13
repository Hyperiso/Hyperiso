#ifndef MAKE_COMPILER_STRATEGY_H
#define MAKE_COMPILER_STRATEGY_H

#include "CompilerStrategy.h"

class MakeCompilerStrategy : public CompilerStrategy {
public:
    void compile_run(const std::string& sourceFile, const std::string& outputBinary) override;
};

#endif // GPP_COMPILER_STRATEGY_H