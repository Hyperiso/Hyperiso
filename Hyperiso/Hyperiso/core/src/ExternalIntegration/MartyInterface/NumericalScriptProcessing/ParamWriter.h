#pragma once
#include <unordered_map>
#include <string>
#include <fstream>
#include "config.hpp"
#include "FileNameManager.h"

class ParamWriter {
private:
    std::unordered_map<std::string, double>& params;

public:
    ParamWriter(std::unordered_map<std::string, double>& params) : params(params) {}

    void writeParams(std::ofstream& outputFile) {
        for (const auto& [name, value] : params) {
            std::string real_name = name;
            outputFile << real_name << "," << value << "\n";
        }
    }

    std::unordered_map<std::string, double>& getParams() {
        return params;
    }
};
