#ifndef GPP_COMPILER_STRATEGY_H
#define GPP_COMPILER_STRATEGY_H

#include "CompilerStrategy.h"
#include "config.hpp"
class GppCompilerStrategy : public CompilerStrategy {
public:
    void compile_run(const std::string& sourceFile, const std::string& outputBinary) override;
    void compile(const std::string& sourceFile, const std::string& outputBinary) override;
    // bool check_if_compile(const std::string& outputBinary) override;
};

#endif // GPP_COMPILER_STRATEGY_H
