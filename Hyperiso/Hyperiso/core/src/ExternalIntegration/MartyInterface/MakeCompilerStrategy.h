#ifndef MAKE_COMPILER_STRATEGY_H
#define MAKE_COMPILER_STRATEGY_H

#include "CompilerStrategy.h"

class MakeCompilerStrategy : public CompilerStrategy {
public:
    MakeCompilerStrategy(std::string model, std::string wilson) : CompilerStrategy(model, wilson) {}
    void compile_run(const std::string& sourceFile, const std::string& outputBinary) override;
    void compile(const std::string& sourceFile, const std::string& outputBinary) override;

    void set_Q_match(double Q_match) {this->Q_match = Q_match;}
    double Q_match = 80.379;

};

#endif // GPP_COMPILER_STRATEGY_H