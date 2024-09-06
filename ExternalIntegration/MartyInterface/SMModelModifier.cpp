#include "SMModelModifier.h"
#include <iostream>

void SMModelModifier::modifyLine(std::string& line) {
    if (line.find("//MODEL_SPECIFIC_CODE") != std::string::npos) {
        line.replace(line.find("//MODEL_SPECIFIC_CODE"), 20, "SM specific code");
    }
}

void SMNumModelModifier::addLine(std::ofstream& outputFile, const std::string& currentLine, bool addBefore) {
    std::cout << "woaw" << std::endl;
    if (currentLine.find("return 0;") != std::string::npos) {
        if (addBefore) {
            outputFile << "\t double mb = 4.18;\n";
            outputFile << "\t double ms = 95e-3;\n";
            outputFile << "\t double s12 = (mb*mb+ms*ms)/2;\n";
            outputFile << "\t param_t param;\n";
            outputFile << "\t param.m_b = mb;\n";
            outputFile << "\t param.m_s = ms;\n";
            outputFile << "\t param.e_em = 0.01;\n";
            outputFile << "\t param.s_12 = s12;\n";
            outputFile << "\t param.theta_W = 2;\n";
            outputFile << "\t auto out = std::ofstream(\"C7_SM.txt\");\n";
            outputFile << "\t param.m_t = 173.2;\n";
            outputFile << "\t out << C7(param).real() << std::endl;\n";
            outputFile << "\t out.close();\n";

        }
        std::cout << currentLine << std::endl;
        // Écrire la ligne courante "return 0;"
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
}