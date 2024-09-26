#include "SMModelModifier.h"
#include <iostream>
#include <complex>
// SMNumModelModifier::SMNumModelModifier()
//     : paramSetter(params) {  // Initialisation du paramSetter avec la référence de params
// }

void SMModelModifier::modifyLine(std::string& line) {
    if (line.find("//MODEL_SPECIFIC_CODE") != std::string::npos) {
        line.replace(line.find("//MODEL_SPECIFIC_CODE"), 20, "SM specific code");
    }
}

void SMNumModelModifier::addLine(std::ofstream& outputFile, const std::string& currentLine, bool addBefore) {
    if (currentLine.find("//42") != std::string::npos) {
        this->done = true;
    }
    if (this->done) {
        outputFile << currentLine << "\n";
        return;
    }
    if(count == 0) {
        outputFile << "//42" << "\n";
        outputFile << currentLine << "\n";
        this->count++;
        return;
    }
    if (currentLine.find("return 0;") != std::string::npos) {
        std::string filename = "libs/" + this->wilson+"_SM/include/params.h";
        auto extractedParams = extractor.extract(filename);
        auto interpretedParams = interpreter.interpret(extractedParams);
        if (addBefore) {
            outputFile << "\tparam_t param;\n";

            

            for (const auto& [name, interpreted] : interpretedParams) {
                paramSetter.setParam(name, interpreted);
            }

            for (const auto& [name, value] : params) {
                std::string real_name = name;
                if (name.find("_im") != std::string::npos) {
                    real_name = name.substr(0, name.size()-3);
                    // std::complex<double> complex_value = params[real_name+"_re"] + std::complex<double>(0,1) * value;
                    outputFile << "\tparam." << real_name << " = {" << params[real_name+"_re"] << "," << value << "};\n";
                } else if (name.find("_re") != std::string::npos) {
                    continue;
                }
                else {
                    outputFile << "\tparam." << real_name << " = " << value << ";\n";
                }
                
            }
            // outputFile << "\tsetMu(81.);\n";
            outputFile << "\tauto out = std::ofstream(\""+ this->wilson +"_SM.txt\");\n";
            outputFile << "\tout << "+ this->wilson + "(param).real() << \" \" << "+ this->wilson + "(param).imag() << std::endl;\n";
            outputFile << "\tout.close();\n";

        }
        std::cout << currentLine << std::endl;
        outputFile << currentLine << "\n";
        if (!addBefore) {
            outputFile << "    // Additional SM specific code\n";
            outputFile << "    int someVariable = 42;\n";
        }
    }
    else if (currentLine.find("using namespace") != std::string::npos){
         outputFile << "#include <fstream>\n";
         outputFile << currentLine << std::endl;
    } else {
        std::cout << currentLine << std::endl;
        outputFile << currentLine << "\n";
    }
    this->count++;
}