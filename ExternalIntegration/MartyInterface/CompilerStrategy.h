#ifndef COMPILER_STRATEGY_H
#define COMPILER_STRATEGY_H

#include <string>
#include <fstream>
#include <sys/stat.h>
#include <iostream>

class CompilerStrategy {
public:
    virtual ~CompilerStrategy() = default;
    virtual void compile_run(const std::string& sourceFile, const std::string& outputBinary) = 0;
    virtual void compile(const std::string& sourceFile, const std::string& outputBinary) = 0;
    virtual bool check_if_compile(const std::string& outputBinary) {
    struct stat buffer;
    if (stat(outputBinary.c_str(), &buffer) != 0) {
        return false;
    }
    if (buffer.st_size == 0) {
        return false;
    }
    std::cout << "Already compiled !" << std::endl;
    return true;
}

protected:
    std::string wilson{};
};

#endif // COMPILER_STRATEGY_H
