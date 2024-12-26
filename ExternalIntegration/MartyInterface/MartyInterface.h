#ifndef MARTY_INTERFACE_H
#define MARTY_INTERFACE_H

#include "CodeGenerator.h"
#include "GppCompilerStrategy.h"
#include "MakeCompilerStrategy.h"
#include "SMModelModifier.h"
#include "Logger.h"
#include "GeneralNumModelModifier.h"
// #include "NumModelModifier.h"

class MartyInterface {
public:
    void generate(std::string wilson, std::string model);
    void compile_run(std::string wilson, std::string model);
    void generate_numlib(std::string wilson, std::string model, double Q_match);
    void compile_run_libs(std::string wilson, std::string model, double Q_match);

    void calculate(std::string wilson, std::string model, double Q_match, bool new_params=false) {
        generate(wilson, model);
        compile_run(wilson, model);
        generate_numlib(wilson, model, Q_match);
        compile_run_libs(wilson, model, Q_match);

    }

private:

    bool already_run(std::string&& outputBinary) {
        struct stat buffer;
        if (stat(outputBinary.c_str(), &buffer) != 0) {
            return false;
        }
        if (buffer.st_size == 0) {
            return false;
        }
        std::cout << "Already run !" << std::endl;
        return true;
    }
    std::string output_binary_name(std::string& wilson, std::string& model) {
        return "generated_" + wilson+"_" + model + ".cpp";
    }

    std::string num_file_path{};
};

#endif // MARTY_INTERFACE_H
