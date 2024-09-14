#include "SMModelModifier.h"
#include <iostream>

SMNumModelModifier::SMNumModelModifier()
    : paramSetter(params) {  // Initialisation du paramSetter avec la référence de params
}

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
        std::string filename = "libs/C7_SM/include/params.h";
        auto extractedParams = extractor.extract(filename);
        auto interpretedParams = interpreter.interpret(extractedParams);
        if (addBefore) {
            // outputFile << "\t double mb = 4.18;\n";
            // outputFile << "\t double ms = 95e-3;\n";
            // outputFile << "\t double s12 = (mb*mb+ms*ms)/2;\n";
            outputFile << "\t param_t param;\n";
            // outputFile << "\t param.m_b = mb;\n";
            // outputFile << "\t param.m_s = ms;\n";
            // outputFile << "\t param.e_em = 0.01;\n";
            // outputFile << "\t param.s_12 = s12;\n";
            // outputFile << "\t param.theta_W = 2;\n";
            

            for (const auto& [name, interpreted] : interpretedParams) {
                paramSetter.setParam(name, interpreted);
            }

            for (const auto& [name, value] : params) {
                outputFile << "\tparam." << name << " = " << value << ";\n";
            }
            outputFile << "\t auto out = std::ofstream(\"C7_SM.txt\");\n";
            // outputFile << "\t param.m_t = 173.2;\n";
            outputFile << "\t out << C7(param).real() << std::endl;\n";
            outputFile << "\t out.close();\n";

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

void SMNumModelModifier::processParams(std::ofstream& outputFile) {
    // Charger et interpréter les paramètres
    std::string filename = "libs/C7_SM/include/params.h";
    auto extractedParams = extractor.extract(filename);
    auto interpretedParams = interpreter.interpret(extractedParams);

    // Utiliser SMParamSetter pour traiter chaque paramètre
    for (const auto& [name, interpreted] : interpretedParams) {
        paramSetter.setParam(name, interpreted);
    }

    // Ajouter les paramètres dans le fichier de sortie
    for (const auto& [name, value] : params) {
        outputFile << "\tparam." << name << " = " << value << ";\n";
    }

    outputFile << "\tparam.m_t = 173.2;\n";
    outputFile << "\tauto out = std::ofstream(\"C7_SM.txt\");\n";
    outputFile << "\tout << C7(param).real() << std::endl;\n";
    outputFile << "\tout.close();\n";
}