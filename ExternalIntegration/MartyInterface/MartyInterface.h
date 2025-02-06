#ifndef MARTY_INTERFACE_H
#define MARTY_INTERFACE_H

#include <memory>
#include <filesystem>
#include <iostream>
#include "config.hpp"
#include <algorithm>
#include "FileNameManager.h"
#include "GeneralModelModifier.h"
#include "CodeGenerator.h"
#include "GppCompilerStrategy.h"
#include "MakeCompilerStrategy.h"
#include "SMModelModifier.h"
#include "Logger.h"
#include "GeneralNumModelModifier.h"

class MartyInterface {
public:
    void generate(std::string wilson, std::string model);
    void compile_run(std::string wilson, std::string model);
    void generate_numlib(std::string wilson, std::string model, double Q_match);
    void compile_run_libs(std::string wilson, std::string model, double Q_match);

    void calculate(std::string wilson, std::string model, double Q_match, bool new_params=false);

private:

    bool already_run(std::string&& outputBinary);
    std::string output_binary_name(std::string& wilson, std::string& model);

    std::string num_file_path{};
};

#endif // MARTY_INTERFACE_H
