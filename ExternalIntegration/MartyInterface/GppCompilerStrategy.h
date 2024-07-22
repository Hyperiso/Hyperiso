#ifndef GPP_COMPILER_STRATEGY_H
#define GPP_COMPILER_STRATEGY_H

#include "CompilerStrategy.h"

class GppCompilerStrategy : public CompilerStrategy {
public:
    void compile(const std::string& sourceFile, const std::string& outputBinary) override;
};

#endif // GPP_COMPILER_STRATEGY_H
