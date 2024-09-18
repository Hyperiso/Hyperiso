#include "MartyInterface.h"
#include "Extractor.h"
#include "Interpreter.h"

int main() {

    MartyInterface MartyInterface;
    // MartyInterface.generate("C7", "SM");
    // MartyInterface.compile_run("C7", "SM");
    MartyInterface.generate_numlib("C8", "SM");
    MartyInterface.compile_run_libs("C8", "SM");

    //  MartyInterface.generate("C8", "SM");
    // MartyInterface.compile_run("C8", "SM");

    // Extractor extractor;
    // Interpreter interpreter;

    // std::string filename = "libs/C7_SM/include/params.h";


    // std::vector<Extractor::Parameter> extractedParams = extractor.extract(filename);

    // auto interpretedParams = interpreter.interpret(extractedParams);

    // for (const auto& [name, interpreted] : interpretedParams) {
    //     std::cout << "Paramètre: " << name << "\n";
    //     std::cout << "Bloc: " << interpreted.block << "\n";
    //     std::cout << "Code: " << interpreted.code << "\n";
    //     std::cout << "-------------------\n";
    // }
    // return 0;
}