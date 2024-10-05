#pragma once
#include <unordered_map>
#include <string>
#include <fstream>
#include "config.hpp"

class ParamWriter {
private:
    std::unordered_map<std::string, double>& params;  // Référence vers les paramètres
    std::string wilson;
    
public:
    ParamWriter(std::unordered_map<std::string, double>& params, const std::string& wilson)
        : params(params), wilson(wilson) {}

    void writeParams(std::ofstream& outputFile) {
        outputFile << "\tstd::string path = \"" << project_root.data() << "/DataBase/MartyWilson/SM_wilson.csv\";\n";
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

    }

    void writeSingleParam(std::ofstream& outputFile, const std::string& paramName, double paramValue) {
        outputFile << "\tparam." << paramName << " = " << paramValue << ";\n";
    }

    std::unordered_map<std::string, double>& getParams() {
        return params;
    }
};
