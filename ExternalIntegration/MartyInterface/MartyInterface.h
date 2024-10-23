#ifndef MARTY_INTERFACE_H
#define MARTY_INTERFACE_H

#include "CodeGenerator.h"
#include "GppCompilerStrategy.h"
#include "MakeCompilerStrategy.h"
#include "SMModelModifier.h"
#include "Logger.h"
#include "NumModelModifier.h"

class MartyInterface {
public:
    void generate(std::string wilson, std::string model);
    void compile_run(std::string wilson, std::string model);
    void generate_numlib(std::string wilson, std::string model);
    void compile_run_libs(std::string wilson, std::string model);

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
