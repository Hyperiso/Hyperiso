#include "MartyInterface.h"
#include "GeneralNumModelModifier.h"
#include "Extractor.h"
#include "Interpreter.h"
#include "IncludeManager.h"
#include "LineProcessor.h"
#include "ModelWriter.h"
#include "SMParamSetter.h"
#include "FileNameManager.h"
#include "MemoryManager.h"
#include "config.hpp"

int main() {
    MemoryManager* mm = MemoryManager::GetInstance();
    std::string base_path = project_root.data();
    mm->init(base_path+"/Test/testinput_thdm.lha", Model::SM);
    // auto truc = FileNameManager::getInstance("C7", "SM");

    // std::cout << truc->getGeneratedFileName() << std::endl;
    // std::cout << truc->getExecutableFileName() << std::endl;
    // std::cout << truc->getNumGeneratedFileName() << std::endl;
    // std::cout << truc->getNumExecutableFileName() << std::endl;
    // std::cout << truc->getHelperFileName("h") << std::endl;

    
    MartyInterface MartyInterface;
    // MartyInterface.generate("C7", "SM");
    // MartyInterface.compile_run("C7", "SM");
    
    // MartyInterface.generate_numlib("C7", "SM");

    // MartyInterface.compile_run_libs("C7", "SM", 81);
    // MartyInterface.generate("C2", "SM");
    // MartyInterface.compile_run("C2", "SM");
    // MartyInterface.generate_numlib("C2", "SM");
    // MartyInterface.compile_run_libs("C2", "SM", 81);

    // MartyInterface.generate("C8", "SM");
    // MartyInterface.compile_run("C8", "SM");
    // MartyInterface.generate_numlib("C8", "SM");
    // MartyInterface.compile_run_libs("C8", "SM", 81);


    // MartyInterface.calculate("C2", "SM", 81);
    // MartyInterface.calculate("C7", "SM", 81);
    // MartyInterface.calculate("C5", "SM", 81);
    MartyInterface.calculate("C7", "THDM", 81);
    MartyInterface.calculate("C7", "SM", 160);
    MartyInterface.calculate("C7", "SM", 81);
    MartyInterface.calculate("C2", "THDM", 81);
    MartyInterface.calculate("C2", "SM", 81);
    // MartyInterface.calculate("C1", "THDM", 81);
    // MartyInterface.calculate("C2", "THDM", 81);
    // MartyInterface.calculate("C3", "THDM", 81);
    // MartyInterface.calculate("C4", "THDM", 81);
    // MartyInterface.calculate("C5", "THDM", 81);
    // MartyInterface.calculate("C6", "THDM", 81);
    // MartyInterface.calculate("C7", "THDM", 81);
    // MartyInterface.calculate("C8", "THDM", 81);
    // MartyInterface.calculate("C9", "THDM", 81);
    // MartyInterface.calculate("C10", "THDM", 81);
    // MartyInterface.generate("C7", "MSSM");
    // MartyInterface.compile_run("C7", "MSSM");

    // MartyInterface.generate("C7", "THDM");
    // MartyInterface.compile_run("C7", "THDM");
    // std::unordered_map<std::string, double> params;
    // Extractor extractor;
    // Interpreter interpreter;
    // SMParamSetter paramSetter(params);

    // std::string wilson = "C7"; 

    // ParamWriter paramWriter(params, wilson);
    // IncludeManager includeManager;
    // LineProcessor lineProcessor(paramWriter, includeManager, /*forceMode=*/false);
    // ModelWriter modelWriter(lineProcessor);

    // std::ifstream inputFile("libs/C7_SM/script/example_c7_sm.cpp");
    // std::ofstream outputFile("libs/C7_SM/script/example_c7_sm2.cpp");
    // GeneralNumModelModifier modelModifier("C7", /*force=*/false);
    // modelModifier.modify(inputFile, outputFile);

    // modelWriter.writeModel(inputFile, outputFile);
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