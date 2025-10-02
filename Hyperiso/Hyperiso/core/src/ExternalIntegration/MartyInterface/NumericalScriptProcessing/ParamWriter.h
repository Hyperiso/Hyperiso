#pragma once
#include <unordered_map>
#include <string>
#include <fstream>
#include "config.hpp"
#include "FileNameManager.h"

class ParamWriter {
public:
    ParamWriter() =default;

    void writeParams(std::ofstream& outputFile, const std::unordered_map<std::string, double>& params) {
        for (const auto& [name, value] : params) {
            std::string real_name = name;
            outputFile << real_name << "," << value << "\n";
        }

    }

};
