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

        add_argpars(outputFile);

        outputFile << "\tparam_t param;\n";
        
        for (const auto& [name, value] : params) {
            std::cout << name << " " << value << std::endl;
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

        outputFile << "\tsetMu(Q_match);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "\", " + wilson + "(param), Q_match, path);\n";

    }

    void add_argpars(std::ofstream& outputFile) {
        
        outputFile << "\tdouble Q_match = 80.379;\n";
        outputFile << "\tfor (int i = 1; i < argc; i++) {\n";
        outputFile << "\t\tif (std::string(argv[i]) == \"--Q_match\" || std::string(argv[i]) == \"-Q\") {\n";
        outputFile << "\t\t\tQ_match = std::stod(argv[i + 1]);\n";
        outputFile << "\t\t\ti++;\n";
        outputFile << "\t\t} else if (std::string(argv[i]) == \"--help\" || std::string(argv[i]) == \"-h\") {\n";
        outputFile << "\t\t\tstd::cout << \"Options disponibles :\" << std::endl;\n";
        outputFile << "\t\t\tstd::cout << \"--Q_match/-Q : Valeur de Q_match (par défaut 80.379)\" << std::endl;\n";
        outputFile << "\t\t\tstd::cout << \"--help/-h : Affiche ce message.\" << std::endl;\n";
        outputFile << "\t\t\treturn 0;\n";
        outputFile << "\t\t}\n";
        outputFile << "\t}\n\n";
        
    }

    void writeSingleParam(std::ofstream& outputFile, const std::string& paramName, double paramValue) {
        outputFile << "\tparam." << paramName << " = " << paramValue << ";\n";
    }

    std::unordered_map<std::string, double>& getParams() {
        return params;
    }
};
