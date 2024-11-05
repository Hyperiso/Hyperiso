#ifndef NUM_MODEL_MODIFIER_H
#define NUM_MODEL_MODIFIER_H

#include <map>
#include <string>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include <sstream>
#include "ModelModifier.h"
#include "Extractor.h"
#include "Interpreter.h"
#include "SMParamSetter.h"
#include "config.hpp"

class NumModelModifier : public ModelModifier {
private:
    std::map<std::string, std::string> paramMap;
    std::unordered_map<std::string, double> params;
    bool done = false;
    bool forceMode = false;
    int count = 0; 
    std::string wilson;
    Extractor extractor;
    Interpreter interpreter;
    SMParamSetter paramSetter; 

public:
    NumModelModifier(const std::string& wilson, bool force = false)
        : wilson(wilson), forceMode(force), paramSetter(params) {
        std::string filename = "libs/" + wilson + "_SM/include/params.h";
        auto extractedParams = extractor.extract(filename);
        auto interpretedParams = interpreter.interpret(extractedParams);

        for (const auto& [name, interpreted] : interpretedParams) {
            paramSetter.setParam(name, interpreted);
        }

    }

    void addLine(std::ofstream& outputFile, const std::string& currentLine, bool addBefore) override {
        if (currentLine.find("//42") != std::string::npos && !forceMode) {
            this->done = true;
        }

        if (done && !forceMode) {
            outputFile << currentLine << "\n";
            return;
        }

        if (currentLine.find("param.") != std::string::npos) {
            std::string paramName = extractParamName(currentLine);
            std::string paramValue = extractParamValue(currentLine);

            if (params.find(paramName) != params.end()) {
                outputFile << "\tparam." << paramName << " = " << params[paramName] << ";\n";
                params.erase(paramName);
            } else {
                outputFile << currentLine << "\n";
            }
        } 
        else if (currentLine.find("return 0;") != std::string::npos) {
            if (addBefore) {
                outputFile << "std::string path = \"" << project_root.data() << "/DataBase/MartyWilson/SM_wilson.csv\";\n";
                outputFile << "\tparam_t param;\n";
                for (const auto& [name, value] : params) {
                    std::string real_name = name;
                    if (name.find("_im") != std::string::npos) {
                        real_name = name.substr(0, name.size() - 3);
                        outputFile << "\tparam." << real_name << " = {"
                                   << params[real_name + "_re"] << ", " << value << "};\n";
                    } else if (name.find("_re") != std::string::npos) {
                        continue;
                    } else {
                        outputFile << "\tparam." << real_name << " = " << value << ";\n";
                    }
                }
                outputFile << "\tsetMu(80.379);\n";
                outputFile << "\twriteWilsonCoefficients(\""+ wilson + "\", " + wilson +"(param), 200, path);\n";
                outputFile << "\tauto out = std::ofstream(\"" + wilson + "_SM.txt\");\n";
                outputFile << "\tout << " + wilson + "(param).real() << \" \" << "
                           + wilson + "(param).imag() << std::endl;\n";
                outputFile << "\tout.close();\n";
            }

            outputFile << currentLine << "\n";
        } 
        else if (currentLine.find("using namespace") != std::string::npos) {
            outputFile << "#include <fstream>\n";
            outputFile << "#include \"csv_helper.h\"\n";
            outputFile << currentLine << "\n";
        } 
        else {
            outputFile << currentLine << "\n";
        }

        count++;
    }

    void modifyLine(std::string& line) override {}
private:
    std::string extractParamName(const std::string& line) {
        size_t startPos = line.find("param.") + 6;
        size_t endPos = line.find(" =");
        return line.substr(startPos, endPos - startPos);
    }

    std::string extractParamValue(const std::string& line) {
        size_t startPos = line.find("= ") + 2;
        size_t endPos = line.find(";");
        return line.substr(startPos, endPos - startPos);
    }
};

#endif