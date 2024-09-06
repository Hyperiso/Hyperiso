#ifndef MARTY_INTERFACE_H
#define MARTY_INTERFACE_H

#include "CodeGenerator.h"
#include "GppCompilerStrategy.h"
#include "MakeCompilerStrategy.h"
#include "SMModelModifier.h"
#include "Logger.h"
// #include "THDMModelModifier.h"

class MartyInterface {
public:
    void generate(std::string wilson, std::string model);
    void compile_run(std::string wilson, std::string model);
    void generate_numlib(std::string wilson, std::string model);
    void compile_run_libs(std::string wilson, std::string model);

private:
    bool generated = false;
    bool compiled = false;

    std::string num_file_path{};
};

#endif // MARTY_INTERFACE_H
