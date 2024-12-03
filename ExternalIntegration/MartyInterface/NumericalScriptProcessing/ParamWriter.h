#pragma once
#include <unordered_map>
#include <string>
#include <fstream>
#include "config.hpp"
#include "FileNameManager.h"

class ParamWriter {
private:
    std::unordered_map<std::string, double>& params;  // Référence vers les paramètres
    std::string wilson;
    std::string model;
    bool top, bottom, charm;

public:
    ParamWriter(std::unordered_map<std::string, double>& params, const std::string& wilson, const std::string& model)
        : params(params), wilson(wilson), model(model) {}

    void writeParams(std::ofstream& outputFile) {
        outputFile << "\tstd::string path = \"" << FileNameManager::getInstance(this->wilson, this->model)->getCsvWilsonFileName() <<"\";\n";


        outputFile << "\tparam_t param;\n";
        
        for (const auto& [name, value] : params) {
            std::cout << "--..--..--..--..--..--.." << std::endl;
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
        add_argpars(outputFile);
        quark_special_case(outputFile);
        outputFile << "\tsetMu(Q_match);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "\", " + wilson + "(param), Q_match, path);\n";

    }

    void add_argpars(std::ofstream& outputFile) {
        
        outputFile << "\tdouble Q_match = 80.379;\n";
        outputFile << "\tdouble m_top = 172;\n";
        outputFile << "\tdouble m_bot = 4.18;\n";
        outputFile << "\tdouble m_char = 1.79;\n";
        outputFile << "\tfor (int i = 1; i < argc; i++) {\n";
        outputFile << "\t\tif (std::string(argv[i]) == \"--Q_match\" || std::string(argv[i]) == \"-Q\") {\n";
        outputFile << "\t\t\tQ_match = std::stod(argv[i + 1]);\n";
        outputFile << "\t\t\ti++;\n";
        outputFile << "\t\t}else if (std::string(argv[i]) == \"--m_top\" || std::string(argv[i]) == \"-mt\") {\n";
        outputFile << "\t\t\tm_top = std::stod(argv[i + 1]);\n";
        outputFile << "\t\t\ti++;\n";
        outputFile << "\t\t}else if (std::string(argv[i]) == \"--m_bot\" || std::string(argv[i]) == \"-mb\") {\n";
        outputFile << "\t\t\tm_bot = std::stod(argv[i + 1]);\n";
        outputFile << "\t\t\ti++;\n";
        outputFile << "\t\t}else if (std::string(argv[i]) == \"--m_char\" || std::string(argv[i]) == \"-mc\") {\n";
        outputFile << "\t\t\tm_char = std::stod(argv[i + 1]);\n";
        outputFile << "\t\t\ti++;\n";
        outputFile << "\t\t} else if (std::string(argv[i]) == \"--help\" || std::string(argv[i]) == \"-h\") {\n";
        outputFile << "\t\t\tstd::cout << \"Options availables :\" << std::endl;\n";
        outputFile << "\t\t\tstd::cout << \"--Q_match/-Q : Value of Q_match (default 80.379)\" << std::endl;\n";
        outputFile << "\t\t\tstd::cout << \"--m_top/-mt : Value of top quark mass (default 172 GeV).\" << std::endl;\n";
        outputFile << "\t\t\tstd::cout << \"--m_bot/-mb : Value of bottom quark mass (default 4.18 GeV).\" << std::endl;\n";
        outputFile << "\t\t\tstd::cout << \"--m_char/-mc : Value of charm quark mass (default 1.79 GeV).\" << std::endl;\n";
        outputFile << "\t\t\tstd::cout << \"--help/-h : Affiche ce message.\" << std::endl;\n";
        outputFile << "\t\t\treturn 0;\n";
        outputFile << "\t\t}\n";
        outputFile << "\t}\n\n";
        
    }

    void quark_special_case(std::ofstream& outputFile) {
        if (this->top) {
            outputFile << "\tparam.m_t = m_top;\n";
        }
        if (this->bottom) {
            outputFile << "\tparam.m_b = m_bot;\n";
        }
        if (this->charm) {
            outputFile << "\tparam.m_c = m_char;\n";
        }
    }

    void writeSingleParam(std::ofstream& outputFile, const std::string& paramName, double paramValue) {
        outputFile << "\tparam." << paramName << " = " << paramValue << ";\n";
    }

    std::unordered_map<std::string, double>& getParams() {
        return params;
    }
};
